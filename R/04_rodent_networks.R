source(here::here("R", "00_setup.R"))
source(here("R", "01_load_data.R"))
source(here("R", "02_clean_data.R"))
source(here("R", "03_descriptive.R"))

# Load in the rodent data
rodents <- unique_rodents %>%
  drop_na(species) %>%
  drop_na(grid_number) %>%
  left_join(trap_data, by = c("trap_id", "village", "grid_number", "visit")) %>%
  group_by(rodent_id) %>%
  arrange(site_id) %>%
  slice(1) %>%
  select(rodent_id, species, village, visit, grid_number, site_id, landuse, trap_easting, trap_northing)

landuse_visit_rodents <- rodents %>%
  st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
  group_by(landuse, visit) %>%
  group_split()

trap_locations <- trap_locations %>%
  left_join(visit_season, by = "visit")

network_numbers <- trap_locations %>%
  distinct(village, grid_number, visit, landuse) %>%
  left_join(landuse_visit_rodents %>%
              bind_rows(.id = "Site") %>%
              mutate(Site = as.numeric(Site)) %>%
              tibble() %>%
              select(Site, village, visit, grid_number, landuse),
            by = c("village", "visit", "grid_number", "landuse")) %>%
  group_by(village, visit, landuse) %>%
  fill(Site, .direction = "downup") %>%
  drop_na(Site) %>%
  rename("Network" = Site,
         "Village" = village,
         "Visit" = visit,
         "Landuse" = landuse,
         "Grid" = grid_number) %>%
  distinct(Network, Village, Visit, Grid, Landuse)

# Rodent species with less than 50 individuals will be grouped as Other
species_n <- tibble(bind_rows(landuse_visit_rodents)) %>%
  group_by(species) %>%
  summarise(N = n()) %>%
  mutate(species = str_to_sentence(str_replace_all(species, "_", " ")))

other_species <- species_n$species[species_n$N <= 50]

if(file.exists(here("data", "expanded_assemblages_2023-03-02.rds"))) {
  
  expanded_assemblages <- read_rds(here("data", "expanded_assemblages_2023-03-02.rds"))
  
} else {
  
  produce_assemblages <- function(rodent_data = landuse_visit_rodents, distance = 30) {
    
    # Define distance as meters
    units(distance) <- "meters"
    
    # Produce a buffer around each individual to create their range
    individual_buffers <- lapply(rodent_data, function(x) {
      x %>%
        select(rodent_id, species, village, visit, grid_number) %>%
        st_buffer(dist = distance) %>%
        group_by(rodent_id) %>%
        group_split() })
    
    assemblages <- lapply(individual_buffers, function(y) {
      
      lapply(y, function(x) {
        
        # Allocate a village
        reference_village <- unique(x$village)
        
        # Allocate a visit to a rodent
        reference_visit <- unique(x$visit)
        
        # Allocate a grid
        reference_grid_number <- unique(x$grid_number)
        
        # Allocate a species to an individual
        reference_species <- unique(x$species)
        
        # Allocate an ID to use for reference
        reference_id <- unique(x$rodent_id)
        
        # Join the reference rodent to others that were detected in their range
        contact <- st_filter(bind_rows(rodent_data) %>%
                               filter(!rodent_id %in% reference_id &
                                        village %in% reference_village &
                                        visit %in% reference_visit &
                                        grid_number %in% reference_grid_number) %>%
                               rename(to_id = rodent_id,
                                      to_species = species), 
                             x %>%
                               select(rodent_id, geometry),
                             .predicate = st_within) %>%
          mutate(from_id = reference_id,
                 from_species = reference_species) %>%
          select(any_of(x = c("from_id", "from_species", "to_id", "to_species", "village", "visit", "grid_number")))
        
        if(nrow(contact) == 0) {
          
          contact <- tibble(from_id = reference_id,
                            from_species = reference_species,
                            to_id = NA,
                            to_species = NA,
                            village = reference_village,
                            visit = reference_visit,
                            grid_number = reference_grid_number)
          
        } else {
          
          contact <- tibble(contact) %>%
            select(-geometry)
          
        }
        
      })
    })
    
    # Combine the assemblages
    combined_assemblages <- lapply(assemblages, function(x) {bind_rows(x, .id = "assemblage") %>%
        group_by(assemblage) %>%
        mutate(assemblage = as.numeric(assemblage))
    })
    
    # Nodes
    nodes = lapply(combined_assemblages, function(x) {unique(c(x$from_id, x %>%
                                                                 drop_na(to_id) %>%
                                                                 pull(to_id)))
    })
    
    vertices <- lapply(nodes, function(x) {
      tibble(node = x) %>%
        left_join(tibble(bind_rows(rodent_data)) %>%
                    select(rodent_id,
                           species,
                           village,
                           visit,
                           grid_number,
                           landuse),
                  by = c("node" = "rodent_id")) %>%
        mutate(Species = str_to_sentence(str_replace_all(species, "_", " ")),
               Species = case_when(Species %in% other_species ~ "Other spp",
                                   TRUE ~ Species),
               Landuse = str_to_sentence(str_remove_all(landuse, "secondary_")),
               Village = str_to_sentence(village)) %>%
        rename(Visit = visit) %>%
        select(node, Species, Village, Visit, Grid_number = "grid_number", Landuse)
    })
    
    # Edgelist
    edgelist <- lapply(combined_assemblages, function(x) {x %>%
        ungroup() %>%
        drop_na(to_id) %>%
        select(from_id, to_id) %>%
        mutate(value = 1)
    })
    
    species_graph <- mapply(x = edgelist, y = vertices, function(x, y) {graph_from_data_frame(d = x, vertices = y, directed = FALSE)
    })
    
    simple_graph <- lapply(species_graph, function(x) {igraph::simplify(x, remove.multiple = TRUE, remove.loops = TRUE)})
    
    return(list(nodelist = vertices,
                edgelist = edgelist,
                graph = simple_graph,
                full_graph = species_graph))
    
  }
  
  assemblages <- produce_assemblages(landuse_visit_rodents, distance = 30)
  
  # Add further vertex attributes
  
  for(i in 1:length(assemblages$graph)) {
    # Species grouping
    vertex_attr(assemblages$graph[[i]], name = "Species", index = V(assemblages$graph[[i]])$Species %in% other_species) <- "Other"
    vertex_attr(assemblages$full_graph[[i]], name = "Species", index = V(assemblages$graph[[i]])$Species %in% other_species) <- "Other"
    
    # Distance of contact definition
    vertex_attr(assemblages$graph[[i]], name = "Contact radius") <- "30m" # set to the distance used in producing the network
    
    # Village classification as Urban or Rural
    vertex_attr(assemblages$graph[[i]], name = "Village locations", index = V(assemblages$graph[[i]])$"Village" == "Lambayama") <- "Urban"
    vertex_attr(assemblages$graph[[i]], name = "Village locations", index = V(assemblages$graph[[i]])$"Village" != "Lambayama") <- "Rural"
    
    # Season classification as wet or dry
    vertex_attr(assemblages$graph[[i]], name = "Season", index = V(assemblages$graph[[i]])$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
    vertex_attr(assemblages$graph[[i]], name = "Season", index = V(assemblages$graph[[i]])$"Visit" %in% c("3", "4", "7", "8")) <- "Wet"
    
  }
  
  # Consistent colours for species
  species_palette_df <- as_tibble(bind_graphs(assemblages$graph)) %>%
    select(Species) %>%
    distinct(Species) %>%
    arrange(Species) %>%
    mutate(colour = brewer.pal(length(.$Species), "Set1"))
  
  species_palette <- c(species_palette_df$colour)
  names(species_palette) <- c(species_palette_df$Species)
  
  write_rds(assemblages$nodelist, here("temp", "network_nodelist.rds"))
  
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
    pivot_longer(cols = contains("_mode"), names_to = "Species", values_to = "Estimated") %>%
    mutate(Species = str_to_sentence(str_replace(str_remove_all(Species, "_mode"), "_", " ")))
  
  expanded_assemblages <- list()
  expanded_assemblages$nodelist <- list()
  
  for(i in 1:length(assemblages$nodelist)) {
    
    nodelist <- formatted_abundance %>%
      filter(Network == i) %>%
      uncount(Estimated) %>%
      group_by(Visit, Village) %>%
      mutate(id = 10000 - row_number()) %>%
      mutate(node = paste0(Visit, "_", str_to_upper(str_sub(Village, 1, 3)), "_", id)) %>%
      select(-id)
    
    expanded_assemblages$nodelist[[i]] <- bind_rows(assemblages$nodelist[[i]] %>%
                                                      mutate(Observed = TRUE,
                                                             Network = i), 
                                                    nodelist %>%
                                                      mutate(Observed = FALSE)) %>%
      select(-Grid_number) %>%
      fill(Landuse, .direction = "down")
    
  }
  
  expanded_assemblages$adjacency <- list()
  
  for(i in 1:length(assemblages$edgelist)) {
    
    # Produce an edgelist for the observed individuals
    # Assign 0 to the trapped individuals with no observed contacts and combine with the individuals with observed contacts
    observed_edgelist <- expanded_assemblages$nodelist[[i]] %>%
      filter(!node %in% assemblages$edgelist[[i]]$from_id) %>%
      filter(!node %in% assemblages$edgelist[[i]]$to_id) %>%
      filter(Observed == TRUE) %>%
      mutate(from_id = node,
             to_id = node,
             value = 0) %>%
      select(from_id, to_id, value) %>%
      bind_rows(assemblages$edgelist[[i]] %>%
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
    
    expanded_network%v%"Species" <- expanded_assemblages$nodelist[[i]]$Species
    expanded_network%v%"Village" <- expanded_assemblages$nodelist[[i]]$Village
    expanded_network%v%"Visit" <- expanded_assemblages$nodelist[[i]]$Visit
    expanded_network%v%"Landuse" <- expanded_assemblages$nodelist[[i]]$Landuse
    expanded_network%v%"Observed" <- expanded_assemblages$nodelist[[i]]$Observed
    
    expanded_assemblages$network[[i]] <- expanded_network
    
  }
  
  write_rds(expanded_assemblages, here("data", "expanded_assemblages_2023-03-02.rds"))
  
}

# Describing contact networks --------------------------------------------------------
rodent_network <- expanded_assemblages$network

write_rds(rodent_network, here("data", "rodent_network_site_visit.rds"))

# 23 networks have been produced one for each landuse type and visit
# These have been expanded with non-observed individuals based on estimated abundance

rodent_network <- read_rds(here("data", "rodent_network_site_visit.rds"))

# Landuse network metrics
`%s%` <- network::`%s%`

rodent_net_descriptives <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  
  descriptives <- tibble(nodes = network.size(x),
                         observed_nodes = table(x%v%"Observed")["TRUE"],
                         unobserved_nodes = table(x%v%"Observed")["FALSE"],
                         non_missing_edges = network.edgecount(x),
                         missing_edges = network.naedgecount(x),
                         mean_degree_observed = mean(degree(observed_subgraph, gmode = "graph")),
                         sd_degree = sd(degree(observed_subgraph, gmode = "graph")),
                         landuse = unique(x%v%"Landuse"),
                         n_species = length(unique(x%v%"Species")),
                         density = network.density(x),
                         density_observed = network.density(x, na.omit = FALSE))
  
}) %>%
  bind_rows() %>%
  mutate(network_number = row_number()) %>%
  relocate(network_number, landuse, n_species) %>%
  mutate(network_number = fct_inorder(as.character(network_number)))

# Plot species level properties -------------------------------------------

species_degree <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  
  degrees <- tibble("Species" = observed_subgraph%v%"Species",
                    "Degree" = degree(observed_subgraph),
                    "Landuse" = observed_subgraph%v%"Landuse")
}) %>%
  bind_rows() %>%
  arrange(-Degree) %>%
  mutate(Species = factor(Species, levels = names(species_palette)),
         Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village")))

species_degree %>%
  ggplot() +
  geom_boxplot(aes(x = fct_rev(Species), y = Degree, fill = Species)) +
  labs(title = "Degree by species",
       x = "Species") +
  scale_fill_manual(values = species_palette) +
  facet_grid(~ Landuse) +
  coord_flip() +
  theme_bw()

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

# Save models -------------------------------------------------------------
dir.create(here("temp"))
write_rds(rodent_network, here("temp", "rodent_networks_2023-03-02.rds"))
write_rds(rodent_models, here("temp", "rodent_models_2023-03-02.rds"))
write_rds(rodent_models_summary, here("temp", "rodent_models_summary_2023-03-02.rds"))
