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
  left_join(season, by = "visit")

produce_assemblages_site <- function(rodent_data = site_rodents, distance = 15) {
  
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
           Landuse = str_to_sentence(str_remove_all(landuse, "secondary_")),
           Village = str_to_sentence(village)) %>%
    rename(Visit = visit) %>%
    select(node, Species, Village, Visit, Grid_number = "grid_number", Landuse)
    })
  
  # Edgelist
  edgelist = lapply(combined_assemblages, function(x) {x %>%
    ungroup() %>%
    drop_na(to_id) %>%
    select(from_id, to_id) 
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

# Rodent species with less than 10 individuals will be grouped as Other
species_n <- tibble(bind_rows(landuse_rodents)) %>%
  group_by(species) %>%
  summarise(N = n()) %>%
  mutate(species = str_to_sentence(str_replace_all(species, "_", " ")))

other_species <- species_n$species[species_n$N <= 10]

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

plot_graph <- list()

plot_graph <- lapply(assemblages_30$graph, function(x) {
  as_tbl_graph(x) %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0(unique(V(x)$Landuse), ", Visit ", unique(V(x)$Visit)),
         caption = "30m") +
    theme_graph(title_size = 16, base_family = "sans")
})


# Describing contact networks --------------------------------------------------------

rodent_network <- list(observed_30 = assemblages_30$graph)

write_rds(rodent_network, here("data", "rodent_network_site.rds"))

rodent_network <- read_rds(here("data", "rodent_network_site.rds"))

# Landuse network metrics
rodent_networks <- lapply(rodent_network, function(x) lapply(x, function(y) {asNetwork(y)}))

rodent_net_descriptives <- lapply(rodent_networks, function(x) {
  
  net_descriptives <- tibble()
  
  for(i in 1:length(x)) {
    
    descriptives <- tibble(network_number = i,
                           nodes = length(x[[i]] %v% "vertex.names"),
                           edges = length(valid.eids(x[[i]])),
                           mean_degree = mean(degree(x[[i]], gmode = "graph")),
                           sd_degree = sd(degree(x[[i]], gmode = "graph")),
                           density = if(is.nan(graph.density(asIgraph(x[[i]])))) NA else graph.density(asIgraph(x[[i]])))
                           
    
    net_descriptives <- net_descriptives %>%
      bind_rows(descriptives)
    
  }
  
  return(net_descriptives)
  })

network_metric_df <- bind_rows(rodent_net_descriptives, .id = "Distance") %>%
  arrange(nodes) %>%
  mutate(network_number = fct_inorder(as.character(network_number)))

node_plot <- network_metric_df %>%
  ggplot() +
  geom_point(aes(x = nodes, y = network_number)) +
  labs(x = "Nodes") +
  theme_bw()

edge_plot <- network_metric_df %>%
  ggplot() +
  geom_point(aes(x = edges, y = network_number)) +
  labs(x = "Edges") +
  theme_bw()

mean_d_plot <- network_metric_df %>%
  ggplot() +
  geom_point(aes(x = mean_degree, y = network_number)) +
  labs(x = "Mean Degree") +
  theme_bw()

density_plot <- network_metric_df %>%
  ggplot() +
  geom_point(aes(x = density, y = network_number)) +
  labs(x = "Density") +
  theme_bw()

# Produce a random graph of similar properties for comparison
# Random graph of similar size and density

random_rodent_networks <- lapply(rodent_networks, function(x) {
  lapply(x, function(y) {
    
    as.network(rgraph(length(y %v% "vertex.names"), m = 1, tprob = graph.density(asIgraph(y)), mode = "graph"), directed = FALSE)
    
  })
})

random_net_descriptives <- lapply(random_rodent_networks, function(x) {
  
  net_descriptives <- tibble()
  
  for(i in 1:length(x)) {
    
    descriptives <- tibble(network_number = i,
                           nodes = length(x[[i]] %v% "vertex.names"),
                           edges = length(valid.eids(x[[i]])),
                           mean_degree = mean(degree(x[[i]], gmode = "graph")),
                           sd_degree = sd(degree(x[[i]], gmode = "graph")),
                           density = if(is.nan(graph.density(asIgraph(x[[i]])))) NA else graph.density(asIgraph(x[[i]])))
    
    
    net_descriptives <- net_descriptives %>%
      bind_rows(descriptives)
    
  }
  
  return(net_descriptives)
})

names(random_net_descriptives) <- c("random_30")

random_metric_df <- bind_rows(random_net_descriptives, .id = "Distance") %>%
  arrange(nodes) %>%
  mutate(network_number = fct_inorder(as.character(network_number)))

random_node_plot <- random_metric_df %>%
  ggplot() +
  geom_point(aes(x = nodes, y = network_number)) +
  labs(x = "Nodes") +
  theme_bw()

random_edge_plot <- random_metric_df %>%
  ggplot() +
  geom_point(aes(x = edges, y = network_number)) +
  labs(x = "Edges") +
  theme_bw()

random_mean_d_plot <- random_metric_df %>%
  ggplot() +
  geom_point(aes(x = mean_degree, y = network_number)) +
  labs(x = "Mean Degree") +
  theme_bw()

random_density_plot <- random_metric_df %>%
  ggplot() +
  geom_point(aes(x = density, y = network_number)) +
  labs(x = "Density") +
  theme_bw()

comparison_plot_n <- plot_grid(plotlist = list(node_plot, random_node_plot), nrow = 2)
comparison_plot_e <- plot_grid(plotlist = list(edge_plot, random_edge_plot), nrow = 2)
comparison_plot_mean_d <- plot_grid(plotlist = list(mean_d_plot, random_mean_d_plot), nrow = 2)
comparison_plot_density <- plot_grid(plotlist = list(density_plot, random_density_plot), nrow = 2)

species_degree <- lapply(rodent_networks, function(x) {
  degrees <- tibble()
  
  lapply(x, function(y) {
    degrees <- tibble("Species" = y %v% "Species",
                      "Degree" = degree(y),
                      "Distance" = y %v% "Contact.radius") })
  }) %>%
  bind_rows() %>%
  arrange(-Degree) %>%
  mutate(Species = factor(Species, levels = names(species_palette))) %>%
  ggplot() +
  geom_boxplot(aes(x = fct_rev(Species), y = Degree, fill = Species)) +
  labs(title = "Mean degree by species",
       x = "Species") +
  scale_fill_manual(values = species_palette) +
  coord_flip() +
  facet_wrap(~ Distance) +
  theme_bw()

# Build models ------------------------------------------------------------
source(here("R", "modified_functions.R"))

# Random graph models
null_models <- lapply(random_rodent_networks, function(x) {
  lapply(x, function(y) {
    if(length(valid.eids(y)) >= 1) ergm(y ~ edges) else NA
  })
})

models_run <- tibble(number = 1:length(null_models$observed_30)) %>%
  rowwise() %>%
  mutate(null_30 = case_when(is.list(null_models$observed_30[[number]]) ~ TRUE,
                             TRUE ~ FALSE))

summary_null_models <- lapply(null_models, function(x) {
  lapply(x, function(y) {
    
    if(is.list(y)) {
      
      coeff <- y$coefficient
      
      if(coeff == Inf) NA else {
        
        try(summary.ergm.david(y)) }} else NA
  })
})

models_run <- models_run %>%
  rowwise() %>%
  mutate(null_30_summary = case_when(is.list(summary_null_models$observed_30[[number]]) ~ TRUE,
                             TRUE ~ FALSE))

# Main effects only models
rodent_models <- list()

rodent_models$main_effects <- lapply(rodent_networks, function(x) {
  lapply(x, function(y) {
    if(length(valid.eids(y)) >= 1) {
      try(ergm(y ~ edges + nodefactor("Species")))
    } else NA
  })
})

rodent_models_summary <- list()

rodent_models_summary$main_effects <- lapply(rodent_models$main_effects, function(x) {
  lapply(x, function(y) {
    
    if(is.list(y)) {
      
      coeff_min <- min(y$coefficient, na.rm = TRUE)
      coeff_max <- max(y$coefficient, na.rm = TRUE)
      
      if(coeff_min == -Inf | coeff_max == Inf) NA else {
        
        try(summary.ergm.david(y)) }} else NA
  })
})

# Homophily model for main effects term
rodent_models$homophily <- lapply(rodent_networks, function(x) {
  lapply(x, function(y) {
    if(length(valid.eids(y) >= 1)) {
      try(ergm(y ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE)))
         }})
})

rodent_models_summary$homophily <- lapply(rodent_models$homophily, function(x) {
  lapply(x, function(y) {
    
    if(is.list(y)) {
      
      coeff_min <- min(y$coefficient, na.rm = TRUE)
      coeff_max <- max(y$coefficient, na.rm = TRUE)
      
      if(coeff_min == -Inf | coeff_max == Inf) NA else {
        
        try(summary.ergm.david(y)) }} else NA
  })
})


# Adding terms, factor for season, match for village location and landuse
rodent_models$additional_terms_1 <- lapply(rodent_networks, function(x) {
  ergm(x ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE) +
         nodefactor("Season", levels = "Dry") +
         nodematch("Village.locations", levels =  "Rural", diff = TRUE) +
         nodematch("Landuse", levels = c("Agriculture", "Village"), diff = TRUE))
})
rodent_models_summary$additional_terms_1 <- lapply(rodent_models$additional_terms_1, function(x) {summary.ergm.david(x)})

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