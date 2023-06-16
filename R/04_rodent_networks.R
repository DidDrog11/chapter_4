source(here::here("R", "00_setup.R"))
source(here("R", "01_load_data.R"))
#source(here("R", "02_clean_data.R"))
source(here("R", "03_descriptive.R"))

# Load in the rodent data
rodents <- unique_rodents %>%
  left_join(trap_data) %>%
  rename(grid = grid_number)

landuse_visit_rodents <- rodents %>%
  st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
  group_by(landuse, visit) %>%
  group_split()
  
network_numbers <- trap_data %>%
  distinct(village, grid = grid_number, visit, landuse) %>%
  left_join(landuse_visit_rodents %>%
              bind_rows(.id = "Site") %>%
              mutate(Site = as.numeric(Site)) %>%
              tibble() %>%
              distinct(Site, village, visit, grid, landuse),
            by = c("village", "visit", "grid", "landuse")) %>%
  group_by(visit, landuse) %>%
  fill(Site, .direction = "downup") %>%
  drop_na(Site) %>%
  rename("Network" = Site,
         "Village" = village,
         "Visit" = visit,
         "Landuse" = landuse,
         "Grid" = grid) %>%
  distinct(Network, Village, Visit, Grid, Landuse)

species_n <- tibble(bind_rows(landuse_visit_rodents)) %>%
  group_by(species) %>%
  summarise(N = n())

included_species <- species_n %>%
  filter(N >= 10) %>%
  arrange(-N) %>%
  pull(species) %>%
  droplevels()

# Label rare species <10 as other
other_species <- species_n %>%
  filter(N < 10) %>%
  pull(species) %>%
  droplevels()

prepared_data <- list(rodents = rodents,
                      trap_data = trap_data %>%
                        rename(grid = grid_number),
                      network_numbers = network_numbers)

write_rds(prepared_data, here("data", "data_for_abundance.rds"))

if(file.exists(here("data", "expanded_assemblages_2023-06-13.rds"))) {
  
  expanded_assemblages <- read_rds(here("data", "expanded_assemblages_2023-06-13.rds"))
  
} else {
  
  produce_assemblages <- function(rodent_data = landuse_visit_rodents, distance = 30) {
    
    # Define distance as meters
    units(distance) <- "meters"
    
    # Produce a buffer around each individual to create their range
    individual_buffers <- rodent_data %>%
      select(rodent_uid, species, village, visit, grid) %>%
      st_buffer(dist = distance) %>%
      group_by(rodent_uid) %>%
      group_split()
    
    assemblages <- lapply(individual_buffers, function(y) {
      
      # Allocate a village
      reference_village <- unique(y$village)
      
      # Allocate a visit to a rodent
      reference_visit <- unique(y$visit)
      
      # Allocate a grid
      reference_grid <- unique(y$grid)
      
      # Allocate a species to an individual
      reference_species <- unique(y$species)
      
      # Allocate an ID to use for reference
      reference_id <- unique(y$rodent_uid)
      
      # Join the reference rodent to others that were detected in their range
      contact <- st_filter(bind_rows(rodent_data) %>%
                             filter(!rodent_uid %in% reference_id &
                                      village %in% reference_village &
                                      visit %in% reference_visit) %>%
                             rename(to_id = rodent_uid,
                                    to_species = species), 
                           y %>%
                             select(rodent_uid, geometry),
                           .predicate = st_within) %>%
        mutate(from_id = reference_id,
               from_species = reference_species) %>%
        select(any_of(x = c("from_id", "from_species", "to_id", "to_species", "village", "visit", "grid")))
      
      if(nrow(contact) == 0) {
        
        contact <- tibble(from_id = reference_id,
                          from_species = reference_species,
                          to_id = NA,
                          to_species = NA,
                          village = reference_village,
                          visit = reference_visit,
                          grid = reference_grid)
        
      } else {
        
        contact <- tibble(contact) %>%
          select(-geometry)
        
      }
      
    })
    
    # Combine the assemblages
    combined_assemblages <- bind_rows(assemblages, .id = "assemblage") %>%
      group_by(assemblage) %>%
      mutate(assemblage = as.numeric(assemblage))
    
    if(nrow(combined_assemblages) == 1) {
      
      nodes = combined_assemblages$from_id
      
    } else {
    # Nodes
    nodes = unique(c(combined_assemblages$from_id,
                     combined_assemblages %>%
                       drop_na(to_id) %>%
                       pull(to_id)))
    }
    
    vertices <- tibble(node = nodes) %>%
      left_join(tibble(bind_rows(rodent_data)) %>%
                  mutate(rodent_uid = as.character(rodent_uid)) %>%
                  select(rodent_uid,
                         species,
                         village,
                         visit,
                         grid,
                         landuse),
                by = c("node" = "rodent_uid")) %>%
      select(node, Species = species, Village = village, Visit = visit, Grid = grid, Landuse = landuse)
    
    # Edgelist
    edgelist <- combined_assemblages %>%
      ungroup() %>%
      drop_na(to_id) %>%
      select(from_id, to_id) %>%
      mutate(value = 1)
    
    species_graph <- graph_from_data_frame(d = edgelist, vertices = vertices, directed = FALSE)
    
    simple_graph <- igraph::simplify(species_graph, remove.multiple = TRUE, remove.loops = TRUE)
    
    return(list(nodelist = vertices,
                edgelist = edgelist,
                graph = simple_graph,
                full_graph = species_graph))
    
  }
  
  assemblages <- lapply(landuse_visit_rodents, function(x) {
    produce_assemblages(x, distance = 30)
  })
  
  # Add further vertex attributes
  
  for(i in 1:length(assemblages)) {
    # Species grouping
    vertex_attr(assemblages[[i]]$graph, name = "Species", index = V(assemblages[[i]]$graph)$Species %in% other_species) <- "Other"
    vertex_attr(assemblages[[i]]$full_graph, name = "Species", index = V(assemblages[[i]]$full_graph)$Species %in% other_species) <- "Other"
    
    # Distance of contact definition
    vertex_attr(assemblages[[i]]$graph, name = "Contact radius") <- "30m" # set to the distance used in producing the network
    
    # Village classification as Urban or Rural
    vertex_attr(assemblages[[i]]$graph, name = "Village locations", index = V(assemblages[[i]]$graph)$"Village" == "Lambayama") <- "Urban"
    vertex_attr(assemblages[[i]]$graph, name = "Village locations", index = V(assemblages[[i]]$graph)$"Village" != "Lambayama") <- "Rural"
    
    # Season classification as wet or dry
    vertex_attr(assemblages[[i]]$graph, name = "Season", index = V(assemblages[[i]]$graph)$"Visit" %in% c("1", "2", "5", "6", "9", "10")) <- "Dry"
    vertex_attr(assemblages[[i]]$graph, name = "Season", index = V(assemblages[[i]]$graph)$"Visit" %in% c("3", "4", "7", "8")) <- "Wet"
    
  }
  
  # Consistent colours for species
  species_palette_df <- as_tibble(bind_graphs(lapply(assemblages, function(x) x$graph))) %>%
    select(Species) %>%
    distinct(Species) %>%
    arrange(Species)
  
  species_palette <- brewer.pal(nrow(species_palette_df), "Paired")
  names(species_palette) <- c(species_palette_df$Species)
  
  nodelists <- lapply(assemblages, function(x) { x$nodelist })
  
  write_rds(nodelists, here("temp", "network_nodelist.rds"))
  
  # To plot these networks uncomment the below script
  #source(here("R", "07_plot_networks.R"))
  
  # Add unobserved individuals ----------------------------------------------
  # The number of unobserved individuals is derived from the estimating abundance script 06_estimating_abundance.R
  # This assumes a closed population with N individuals at each site at all time points
  
  estimated_abundance <- if(file.exists(here("data", "estimated_population.rds"))) {
    read_rds(here("data", "estimated_population.rds")) 
  } else {
    source(here("R", "06_estimating_abundance.R"))
    read_rds(here("data", "estimated_population.rds")) 
  }
  
  
  formatted_abundance <- estimated_abundance %>%
    pivot_longer(cols = contains("_median"), names_to = "Species", values_to = "Estimated") %>%
    mutate(Species = str_to_sentence(str_replace(str_remove_all(Species, "_median"), "_", " ")),
           Species = fct(Species, levels = c(levels(included_species), "Other")))
  
  expanded_assemblages <- list()
  expanded_assemblages$nodelist <- list()
  
  for(i in 1:length(assemblages)) {
    
    n_observed <- assemblages[[i]]$nodelist %>%
      mutate(Observed = TRUE,
             Network = i) %>%
      group_by(Species, Village) %>%
      summarise(Observed = n(), .groups = "drop")
    
    nodelist <- formatted_abundance %>%
      filter(Network == i) %>%
      left_join(n_observed, by = c("Village", "Species")) %>%
      mutate(Estimated = case_when(!is.na(Observed) ~ Estimated - Observed,
                                   Estimated - Observed < 0 ~ 0,
                                   TRUE ~ Estimated)) %>%
      select(-Observed) %>%
      uncount(Estimated) %>%
      group_by(Visit, Village) %>%
      mutate(id = 10000 - row_number()) %>%
      mutate(node = paste0(Visit, "_", str_to_upper(str_sub(Village, 1, 3)), "_", id)) %>%
      select(-id) %>%
      ungroup()
    
    expanded_assemblages$nodelist[[i]] <- bind_rows(assemblages[[i]]$nodelist %>%
                                                      mutate(Observed = TRUE,
                                                             Network = i), 
                                                    nodelist %>%
                                                      mutate(Observed = FALSE)) %>%
      select(-Grid) %>%
      fill(Landuse, .direction = "down")
    
  }
  
  expanded_assemblages$adjacency <- list()
  
  for(i in 1:length(assemblages)) {
    
    # Produce an edgelist for the observed individuals
    # Assign 0 to the trapped individuals with no observed contacts and combine with the individuals with observed contacts
    observed_edgelist <- expanded_assemblages$nodelist[[i]] %>%
      filter(!node %in% assemblages[[i]]$edgelist$from_id) %>%
      filter(!node %in% assemblages[[i]]$edgelist$to_id) %>%
      filter(Observed == TRUE) %>%
      mutate(from_id = node,
             to_id = node,
             value = 0) %>%
      select(from_id, to_id, value) %>%
      bind_rows(assemblages[[i]]$edgelist %>%
                  mutate(from_id = as.character(from_id),
                         to_id = as.character(to_id)))
    # This produces an edgelist for all observed individuals
    
    # Produce an adjacency matrix from all individuals with values for the observed individuals produced from the edgelist
    adjacency_matrix <- as_adjacency_matrix(graph_from_data_frame(observed_edgelist, vertices = expanded_assemblages$nodelist[[i]]), sparse = FALSE)
    # Remove self-loops
    diag(adjacency_matrix) <- 0
    # To indicate missing edges turn all unobserved individuals into NA based on the colname and then the rowname
    adjacency_matrix[, !colnames(adjacency_matrix) %in% c(observed_edgelist$from_id, observed_edgelist$to_id)] <- NA
    adjacency_matrix[!rownames(adjacency_matrix) %in% c(observed_edgelist$from_id, observed_edgelist$to_id), ] <- NA
    
    # Assign this to the adjacency matrix
    expanded_assemblages$adjacency[[i]] <- adjacency_matrix
  } 
  
  expanded_assemblages$network <- list()
  
  for(i in 1:length(expanded_assemblages$adjacency)) {
    
    # Remove unobserved species when no individuals of that species have been trapped for that network
    obs_species <- expanded_assemblages$nodelist[[i]] %>%
      filter(Observed == TRUE) %>%
      pull(Species)
    
    unobserved_ids <- expanded_assemblages$nodelist[[i]] %>%
      filter(!Species %in% obs_species) %>%
      pull(node)
    
    filtered_network <- expanded_assemblages$adjacency[[i]][!row.names(expanded_assemblages$adjacency[[i]]) %in% unobserved_ids,
                                                            !colnames(expanded_assemblages$adjacency[[i]]) %in% unobserved_ids]
    
    expanded_network <- network(expanded_assemblages$adjacency[[i]], directed = FALSE)
    
    if(any(c(network.vertex.names(expanded_network)) == expanded_assemblages$nodelist[[i]]$node) == FALSE) warning("Node dataframe is not in the same order as the edgelist. Needs to be checked to ensure wrong vertex attributes aren't assigned.")
    
    expanded_network%v%"Species" <- as.character(expanded_assemblages$nodelist[[i]]$Species)
    expanded_network%v%"Village" <- expanded_assemblages$nodelist[[i]]$Village
    expanded_network%v%"Visit" <- as.numeric(expanded_assemblages$nodelist[[i]]$Visit)
    expanded_network%v%"Landuse" <- as.character(expanded_assemblages$nodelist[[i]]$Landuse)
    expanded_network%v%"Observed" <- expanded_assemblages$nodelist[[i]]$Observed
    
    expanded_assemblages$network[[i]] <- expanded_network
    
  }
  
  write_rds(expanded_assemblages, here("data", "expanded_assemblages_2023-06-13.rds"))
  
}

rodent_network <- expanded_assemblages$network

write_rds(rodent_network, here("data", "rodent_network_site_visit.rds"))

# Build models ------------------------------------------------------------
source(here("R", "modified_functions.R"))

# Collapse all non-M. natalensis into Other spp
rodent_network <- lapply(rodent_network, function(x) {
  network::set.vertex.attribute(x, "Species", v = c(which(x%v%"Species" != "Mastomys natalensis")), "Other spp")
})

# Main effects only models
rodent_models <- list()

rodent_models$null_model <- lapply(rodent_network, function(x) {
  if(length(valid.eids(x)) - network.naedgecount(x) >= 1) {
    try(ergm(x ~ edges))
  } else NA
})

rodent_models_summary <- list()

rodent_models_summary$null_model <- lapply(rodent_models$null_model, function(x) {
  
  if(is.list(x)) { try(summary.ergm.david(x)) } else NA
  
})

# M. natalensis only
rodent_models$main_effects <- lapply(rodent_network, function(x) {
  if(length(valid.eids(x)) - network.naedgecount(x) >= 1) {
    try(ergm(x ~ edges + nodefactor("Species", levels = c("Mastomys natalensis"))))
  } else NA
})

rodent_models_summary$main_effects <- lapply(rodent_models$main_effects, function(x) {
  
  if(is.list(x)) { try(summary.ergm.david(x)) } else NA
  
})

# Homophily model for main effects term
rodent_models$homophily <- lapply(rodent_network, function(x) {
  if(length(valid.eids(x)) - network.naedgecount(x) >= 1) {
    try(ergm(x ~ edges + nodefactor("Species", levels = c("Mastomys natalensis")) 
             + nodematch("Species", diff = TRUE, levels = c("Mastomys natalensis"))))
  }
})

rodent_models_summary$homophily <- lapply(rodent_models$homophily, function(x) {
  
  if(is.list(x)) { try(summary.ergm.david(x)) } else NA
  
})

rodent_models$gof <- lapply(rodent_models$homophily, function(x) {
  if(is.list(x)) { try(gof(x)) } else NA
})

final_model <- list(final_model = rodent_models_summary$homophily,
                    final_model_gof = rodent_models$gof)

# Save models -------------------------------------------------------------
dir.create(here("temp"))
write_rds(rodent_network, here("temp", "rodent_networks_2023-06-13.rds"))
write_rds(rodent_models, here("temp", "rodent_models_2023-06-13.rds"))
write_rds(rodent_models_summary, here("temp", "rodent_models_summary_2023-06-13.rds"))
write_rds(final_model, here("temp", "final_model_2023-06-13.rds"))
