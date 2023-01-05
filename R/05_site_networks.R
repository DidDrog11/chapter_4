source(here::here("R", "00_setup.R"))
source(here("R", "01_load_data.R"))
source(here("R", "02_clean_data.R"))
source(here("R", "03_descriptive.R"))
estimated_abundance <- if(file.exists(here("data", "estimated_population.rds"))) {
  read_rds(here("data", "estimated_population.rds")) 
} else {
  source(here("R", "06_estimating_abundance.R"))
}

# Load in the rodent data and split by landuse and visit
rodents <- unique_rodents %>%
  drop_na(species) %>%
  drop_na(grid_number) %>%
  left_join(trap_data, by = c("trap_id", "village", "grid_number", "visit")) %>%
  group_by(rodent_id) %>%
  arrange(site_id) %>%
  slice(1) %>%
  select(rodent_id, species, village, visit, grid_number, site_id, landuse, trap_northing, trap_easting)
  
landuse_rodents <- rodents %>%
  st_as_sf(coords = c("trap_northing", "trap_easting"), crs = SL_UTM) %>%
  group_by(landuse, visit) %>%
  group_split()

trap_locations <- trap_locations %>%
  left_join(visit_season, by = "visit")

# Rodent species with less than 10 individuals will be grouped as Other
species_n <- tibble(bind_rows(landuse_rodents)) %>%
  group_by(species) %>%
  summarise(N = n()) %>%
  mutate(species = str_to_sentence(str_replace_all(species, "_", " ")))

other_species <- species_n$species[species_n$N <= 10]

produce_assemblages_site <- function(rodent_data = landuse_rodents, distance = 30) {
  
  # Define distance as meters
  units(distance) <- "meters"
  
  # Produce a buffer around each individual to create their range
  individual_buffers <- lapply(rodent_data, function(x) { x %>%
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
  edgelist = lapply(combined_assemblages, function(x) {x %>%
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

assemblages_30 <- produce_assemblages_site(landuse_rodents, distance = 30)

# Add further vertex attributes

for(i in 1:length(assemblages_30$graph)) {
  # Species grouping
  vertex_attr(assemblages_30$graph[[i]], name = "Species", index = V(assemblages_30$graph[[i]])$Species %in% other_species) <- "Other"
  vertex_attr(assemblages_30$full_graph[[i]], name = "Species", index = V(assemblages_30$graph[[i]])$Species %in% other_species) <- "Other"
  
  # Distance of contact definition
  vertex_attr(assemblages_30$graph[[i]], name = "Contact radius") <- "30m"
  
  # Village classification as Urban or Rural
  vertex_attr(assemblages_30$graph[[i]], name = "Village locations", index = V(assemblages_30$graph[[i]])$"Village" == "Lambayama") <- "Urban"
  vertex_attr(assemblages_30$graph[[i]], name = "Village locations", index = V(assemblages_30$graph[[i]])$"Village" != "Lambayama") <- "Rural"
  
  # Season classification as wet or dry
  vertex_attr(assemblages_30$graph[[i]], name = "Season", index = V(assemblages_30$graph[[i]])$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
  vertex_attr(assemblages_30$graph[[i]], name = "Season", index = V(assemblages_30$graph[[i]])$"Visit" %in% c("3", "4", "7", "8")) <- "Wet"
  
  }

# Consistent colours for species
species_palette_df <- as_tibble(bind_graphs(assemblages_30$graph)) %>%
  select(Species) %>%
  distinct(Species) %>%
  arrange(Species) %>%
  mutate(colour = brewer.pal(length(.$Species), "Dark2"))

species_palette <- c(species_palette_df$colour)
names(species_palette) <- c(species_palette_df$Species)

# To plot these networks uncomment the below script
#source(here("R", "07_plot_networks.R"))

# Add unobserved individuals ----------------------------------------------
# The number of unobserved individuals is derived from the estimating abundance script 06_estimating_abundance.R
# This assumes a closed population with N individuals at each site at all time points
# This needs to be translated to fit within the structure for these networks, i.e. multiple sites for a single landuse type at different timepoints

formatted_abundance <- estimated_abundance %>%
  mutate(Landuse = str_to_sentence(landuse),
         Village = str_to_sentence(village),
         Grid_number = grid_number) %>%
  pivot_longer(cols = contains("_mode"), names_to = "Species", values_to = "Estimated") %>%
  select(Species, Village, Grid_number, Landuse, Estimated) %>%
  mutate(Species = str_to_sentence(str_replace(str_remove_all(Species, "_mode"), "_", " ")))
  
expanded_assemblages <- list()
expanded_assemblages$nodelist <- list()

for(i in 1:length(assemblages_30$nodelist)) {
  
  visit_nodelist <- unique(assemblages_30$nodelist[[i]]$Visit)
  
  nodelist <- formatted_abundance %>%
    filter(Landuse %in% assemblages_30$nodelist[[i]]$Landuse) %>%
    filter(Species %in% assemblages_30$nodelist[[i]]$Species) %>%
    filter(Village %in% assemblages_30$nodelist[[i]]$Village) %>%
    filter(Grid_number %in% assemblages_30$nodelist[[i]]$Grid_number) %>%
    full_join(assemblages_30$nodelist[[i]] %>%
                group_by(Species, Village, Grid_number) %>%
                summarise(Observed = n())) %>%
    mutate(Observed = replace_na(Observed, 0),
           to_add = Estimated - Observed) %>%
    uncount(to_add) %>%
    select(-Estimated, -Observed) %>%
    mutate(id = 999 - row_number()) %>%
    mutate(node = paste0(visit_nodelist, "_", str_to_upper(str_sub(Village, 1, 3)), "_", id),
           Visit = as.numeric(visit_nodelist)) %>%
    select(-id)
  
  expanded_assemblages$nodelist[[i]] <- bind_rows(assemblages_30$nodelist[[i]] %>%
                                                    mutate(Observed = TRUE), 
                                                  nodelist %>%
                                                    mutate(Observed = FALSE))
  
}

expanded_assemblages$adjacency <- list()

for(i in 1:length(assemblages_30$edgelist)) {
  
  # Produce an edgelist for the observed individuals
  # Assign 0 to the trapped individuals with no observed contacts and combine with the individuals with observed contacts
  observed_edgelist <- expanded_assemblages$nodelist[[i]] %>%
    filter(!node %in% assemblages_30$edgelist[[i]]$from_id) %>%
    filter(!node %in% assemblages_30$edgelist[[i]]$to_id) %>%
    filter(Observed == TRUE) %>%
    mutate(from_id = node,
           to_id = node,
           value = 0) %>%
    select(from_id, to_id, value) %>%
    bind_rows(assemblages_30$edgelist[[i]] %>%
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
  
  expanded_network <- network(expanded_assemblages$adjacency[[i]], directed = FALSE)
  
  if(any(c(network.vertex.names(expanded_network)) == expanded_assemblages$nodelist[[i]]$node) == FALSE) warning("Node dataframe is not in the same order as the edgelist. Needs to be checked to ensure wrong vertex attributes aren't assigned.")
  
  expanded_network%v%"Species" <- expanded_assemblages$nodelist[[i]]$Species
  expanded_network%v%"Village" <- expanded_assemblages$nodelist[[i]]$Village
  expanded_network%v%"Visit" <- expanded_assemblages$nodelist[[i]]$Visit
  expanded_network%v%"Grid_number" <- expanded_assemblages$nodelist[[i]]$Grid_number
  expanded_network%v%"Landuse" <- expanded_assemblages$nodelist[[i]]$Landuse
  expanded_network%v%"Observed" <- expanded_assemblages$nodelist[[i]]$Observed
  
  expanded_assemblages$network[[i]] <- expanded_network
  
}

# Describing contact networks --------------------------------------------------------

rodent_network <- expanded_assemblages$network

write_rds(rodent_network, here("data", "rodent_network_site.rds"))

# 22 site networks have been produced
# These have been expanded with non-observed individuals based on estimated abundance

rodent_network <- read_rds(here("data", "rodent_network_site.rds"))

# Landuse network metrics

rodent_net_descriptives <- lapply(rodent_network, function(x) {
    
    observed <- x%v%"Observed"
    observed_subgraph <- x %s% which(observed == TRUE)
    
    descriptives <- tibble(nodes = length(x %v% "vertex.names"),
                           observed_nodes = table(x%v%"Observed")["TRUE"],
                           unobserved_nodes = table(x%v%"Observed")["FALSE"],
                           non_missing_edges = length(valid.eids(x)) - network.naedgecount(x),
                           missing_edges = network.naedgecount(x),
                           mean_degree_observed = mean(degree(observed_subgraph, gmode = "graph")),
                           sd_degree = sd(degree(observed_subgraph, gmode = "graph")),
                           landuse = unique(x%v%"Landuse"),
                           visit =  unique(x%v%"Visit"),
                           n_species = length(unique(x%v%"Species")),
                           density = gden(x))
    
  }) %>%
  bind_rows() %>%
  mutate(network_number = row_number()) %>%
  relocate(network_number, visit, landuse, n_species)

network_metric_df <- rodent_net_descriptives %>%
  mutate(network_number = fct_inorder(as.character(network_number)))

node_plot <- network_metric_df %>%
  ggplot() +
  geom_point(aes(x = nodes, y = network_number, colour = landuse)) +
  labs(x = "Nodes",
       colour = "Landuse",
       y = element_blank(),
       title = "Nodes") +
  scale_colour_manual(values = landuse_palette) +
  theme_bw()

edge_plot <- network_metric_df %>%
  rowwise() %>%
  mutate(edges = sum(missing_edges, non_missing_edges)) %>%
  ggplot() +
  geom_point(aes(x = edges, y = network_number, colour = landuse)) +
  labs(x = "Edges",
       colour = "Landuse",
       title = "Edges",
       y = element_blank()) +
  scale_colour_manual(values = landuse_palette) +
  theme_bw()

mean_d_plot <- network_metric_df %>%
  ggplot() +
  geom_point(aes(x = mean_degree_observed, y = network_number, colour = landuse)) +
  labs(x = "Mean degree",
       colour = "Landuse",
       title = "Degree (observed)",
       y = element_blank()) +
  scale_colour_manual(values = landuse_palette) +
  theme_bw()

density_plot <- network_metric_df %>%
  mutate(density = replace_na(density, 0)) %>%
  ggplot() +
  geom_point(aes(x = density, y = network_number, colour = landuse)) +
  scale_colour_manual(values = landuse_palette) +
  labs(x = "Density",
       y = element_blank(),
       colour = "Landuse",
       title = "Network density") +
  theme_bw()

save_plot(plot = plot_grid(plotlist = list(node_plot,
                                           edge_plot,
                                           mean_d_plot,
                                           density_plot),
                           nrow = 4),
          filename = here("output", "network_descriptives.png"),
          base_width = 8,
          base_height = 12)


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

# Main effects only models
rodent_models <- list()

rodent_models$null_model <- lapply(rodent_network, function(x) {
  if(length(valid.eids(x)) - network.naedgecount(x) >= 1) {
    try(ergm(x ~ edges))
  } else NA
})

rodent_models$main_effects <- lapply(rodent_network, function(x) {
  if(length(valid.eids(x)) - network.naedgecount(x) >= 1) {
    try(ergm(x ~ edges + nodefactor("Species")))
  } else NA
})

rodent_models_summary <- list()

rodent_models_summary$null_model <- lapply(rodent_models$null_model, function(x) {
  
  if(is.list(x)) { try(summary.ergm.david(x)) } else NA
  
})

rodent_models_summary$main_effects <- lapply(rodent_models$main_effects, function(x) {
  
  if(is.list(x)) { try(summary.ergm.david(x)) } else NA
  
})

# Homophily model for main effects term
rodent_models$homophily <- lapply(rodent_network, function(x) {
  if(length(valid.eids(x)) - network.naedgecount(x) >= 1) {
    try(ergm(x ~ edges + nodefactor("Species") + nodematch("Species", diff = TRUE)))
  }
})

rodent_models_summary$homophily <- lapply(rodent_models$homophily, function(x) {
  
  if(is.list(x)) { try(summary.ergm.david(x)) } else NA
  
})

# Adding terms, factor for season, match for village location and landuse
rodent_models$additional_terms_1 <- lapply(rodent_networks, function(x) {
  lapply(x, function(y) {
    if(length(valid.eids(y) >= 1)) {
      try(ergm(y ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE) +
         nodematch("Village.locations", levels =  "Rural", diff = TRUE) +
         nodematch("Landuse", levels = c("Agriculture", "Village"), diff = TRUE)))
      }})
})


# Current point -----------------------------------------------------------
write_rds(rodent_network, here("temp", "rodent_networks_2023-01-04.rds"))
write_rds(rodent_models, here("temp", "rodent_models_2023-01-04.rds"))
write_rds(rodent_models_summary, here("temp", "rodent_models_summary_2023-01-04.rds"))


rodent_models_summary$additional_terms_1 <- lapply(rodent_models$additional_terms_1, function(x) {
  lapply(x, function(y) {
    
    if(is.list(y)) {
      
      coeff_min <- min(y$coefficient, na.rm = TRUE)
      coeff_max <- max(y$coefficient, na.rm = TRUE)
      
      if(coeff_min == -Inf | coeff_max == Inf) NA else {
        
        try(summary.ergm.david(y)) }} else NA
  })
})


# Adding terms, match for season, village location and landuse
rodent_models$additional_terms_2 <- lapply(rodent_networks, function(x) {
  ergm(x ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE) +
         nodematch("Season", levels = "Dry") +
         nodematch("Village.locations", levels =  "Rural", diff = TRUE) +
         nodematch("Landuse", levels = c("Agriculture", "Village"), diff = TRUE))
})
rodent_models_summary$additional_terms_2 <- lapply(rodent_models$additional_terms_2, function(x) {summary.ergm.david(x)})

# Adding an additional term for isolates
rodent_models$additional_terms_3 <- lapply(rodent_networks, function(x) {
  ergm(x ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE) +
         nodematch("Season", levels = "Dry") +
         nodematch("Village.locations", levels =  "Rural", diff = TRUE) +
         nodematch("Landuse", levels = c("Agriculture", "Village"), diff = TRUE) +
         isolates())
})
rodent_models_summary$additional_terms_3 <- lapply(rodent_models$additional_terms_3, function(x) {summary.ergm.david(x)})

# Adding an additional term for triangles
rodent_models$additional_terms_4 <- lapply(rodent_networks, function(x) {
  ergm(x ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE) +
         nodematch("Season", levels = "Dry") +
         nodematch("Village.locations", levels =  "Rural", diff = TRUE) +
         nodematch("Landuse", levels = c("Agriculture", "Village"), diff = TRUE) +
         triangles())
})
rodent_models_summary$additional_terms_4 <- lapply(rodent_models$additional_terms_4, function(x) {summary.ergm.david(x)})

# Assessing Goodness-of-Fit -----------------------------------------------
# For the 15 meter contact
models_15m <- list(observed = rodent_networks$observed_15,
                   null = null_models$random_15,
                   main_effects = rodent_models$main_effects$observed_15,
                   homophily = rodent_models$homophily$observed_15,
                   additional_terms_1 = rodent_models$additional_terms_1$observed_15,
                   additional_terms_2 = rodent_models$additional_terms_2$observed_15,
                   additional_terms_3 = rodent_models$additional_terms_3$observed_15)

simulated_15m <- list(null_sim = simulate(models_15m$null, seed = 123),
                      main_sim = simulate(models_15m$main_effects, seed = 123),
                      homophily = simulate(models_15m$homophily, seed = 123),
                      additional_terms_1 = simulate(models_15m$additional_terms_1, seed = 123),
                      additional_terms_2 = simulate(models_15m$additional_terms_2, seed = 123),
                      additional_terms_3 = simulate(models_15m$additional_terms_3, seed = 123))

gof_15m <- bind_cols(model = names(models_15m),
                     rbind(summary(models_15m$observed ~ edges + degree(0:5) + triangle),
                           summary(simulated_15m$null_sim ~ edges + degree(0:5) + triangle),
                           summary(simulated_15m$main_sim ~ edges + degree(0:5) + triangle),
                           summary(simulated_15m$homophily ~ edges + degree(0:5) + triangle),
                           summary(simulated_15m$additional_terms_1 ~ edges + degree(0:5) + triangle),
                           summary(simulated_15m$additional_terms_2 ~ edges + degree(0:5) + triangle),
                           summary(simulated_15m$additional_terms_3 ~ edges + degree(0:5) + triangle)),
                     AIC = c(NA, unlist(lapply(models_15m[-1], function(x) AIC(x)[1]))),
                     BIC = c(NA, unlist(lapply(models_15m[-1], function(x) BIC(x)[1]))))

# For the 30 meter contact
models_30m <- list(observed = rodent_networks$observed_30,
                   null = null_models$random_30,
                   main_effects = rodent_models$main_effects$observed_30,
                   homophily = rodent_models$homophily$observed_30,
                   additional_terms_1 = rodent_models$additional_terms_1$observed_30,
                   additional_terms_2 = rodent_models$additional_terms_2$observed_30,
                   additional_terms_3 = rodent_models$additional_terms_3$observed_30)

simulated_30m <- list(null_sim = simulate(models_30m$null, seed = 123),
                      main_sim = simulate(models_30m$main_effects, seed = 123),
                      homophily = simulate(models_30m$homophily, seed = 123),
                      additional_terms_1 = simulate(models_30m$additional_terms_1, seed = 123),
                      additional_terms_2 = simulate(models_30m$additional_terms_2, seed = 123),
                      additional_terms_3 = simulate(models_30m$additional_terms_3, seed = 123))

gof_30m <- bind_cols(model = names(models_30m),
                     rbind(summary(models_30m$observed ~ edges + degree(0:5) + triangle),
                           summary(simulated_30m$null_sim ~ edges + degree(0:5) + triangle),
                           summary(simulated_30m$main_sim ~ edges + degree(0:5) + triangle),
                           summary(simulated_30m$homophily ~ edges + degree(0:5) + triangle),
                           summary(simulated_30m$additional_terms_1 ~ edges + degree(0:5) + triangle),
                           summary(simulated_30m$additional_terms_2 ~ edges + degree(0:5) + triangle),
                           summary(simulated_30m$additional_terms_3 ~ edges + degree(0:5) + triangle)),
                     AIC = c(NA, unlist(lapply(models_30m[-1], function(x) AIC(x)[1]))),
                     BIC = c(NA, unlist(lapply(models_30m[-1], function(x) BIC(x)[1]))))

# For the 50 meter contact
models_50m <- list(observed = rodent_networks$observed_50,
                   null = null_models$random_50,
                   main_effects = rodent_models$main_effects$observed_50,
                   homophily = rodent_models$homophily$observed_50,
                   additional_terms_1 = rodent_models$additional_terms_1$observed_50,
                   additional_terms_2 = rodent_models$additional_terms_2$observed_50,
                   additional_terms_3 = rodent_models$additional_terms_3$observed_50)

simulated_50m <- list(null_sim = simulate(models_50m$null, seed = 123),
                      main_sim = simulate(models_50m$main_effects, seed = 123),
                      homophily = simulate(models_50m$homophily, seed = 123),
                      additional_terms_1 = simulate(models_50m$additional_terms_1, seed = 123),
                      additional_terms_2 = simulate(models_50m$additional_terms_2, seed = 123),
                      additional_terms_3 = simulate(models_50m$additional_terms_3, seed = 123))

gof_50m <- bind_cols(model = names(models_50m),
                     rbind(summary(models_50m$observed ~ edges + degree(0:5) + triangle),
                           summary(simulated_50m$null_sim ~ edges + degree(0:5) + triangle),
                           summary(simulated_50m$main_sim ~ edges + degree(0:5) + triangle),
                           summary(simulated_50m$homophily ~ edges + degree(0:5) + triangle),
                           summary(simulated_50m$additional_terms_1 ~ edges + degree(0:5) + triangle),
                           summary(simulated_50m$additional_terms_2 ~ edges + degree(0:5) + triangle),
                           summary(simulated_50m$additional_terms_3 ~ edges + degree(0:5) + triangle)),
                     AIC = c(NA, unlist(lapply(models_50m[-1], function(x) AIC(x)[1]))),
                     BIC = c(NA, unlist(lapply(models_30m[-1], function(x) BIC(x)[1]))))
gof(models_50m$additional_terms_3)