# Load in the rodent data
unique_rodents <- read_rds(here("data", "unique_rodents.rds")) %>%
  drop_na(initial_species_id) %>%
  drop_na(grid_number) %>%
  st_as_sf(crs = default_CRS) %>%
  st_transform(crs = SL_UTM)

# Calculate the area trapped for each grid in each village
trap_sites <- trap_locations %>%
  group_by(village, simple_habitat, visit, grid_number) %>%
  group_split()

site_area <- lapply(trap_sites, function(x) {
  
  st_as_sf(x, coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
    st_union() %>%
    st_convex_hull() %>%
    st_area()
  
})

produce_assemblages <- function(rodent_data = unique_rodents, distance = 50) {
  
  # Define distance as meters
  units(distance) <- "meters"
  
  # Produce a buffer around each individual to create their range
  individual_buffers <- rodent_data %>%
    select(rodent_uid, initial_species_id, visit) %>%
    st_buffer(dist = distance) %>%
    group_by(rodent_uid) %>%
    group_split()
  
  assemblages <- lapply(individual_buffers, function(x) {
    
    # Allocate a visit to a rodent
    reference_visit <- unique(x$visit)
    
    # Allocate a species to an individual
    reference_species <- unique(x$initial_species_id)
    
    # Allocate an ID to use for reference
    reference_id <- unique(x$rodent_uid)
    
    # Join the reference rodent to others that were detected in their range
    contact <- st_filter(rodent_data %>%
                           filter(!rodent_uid %in% reference_id) %>%
                           rename(to_uid = rodent_uid,
                                  to_species = initial_species_id), 
                         x %>%
                           select(rodent_uid, geometry),
                         .predicate = st_within) %>%
      mutate(from_uid = reference_id,
             from_species = reference_species) %>%
      filter(visit == reference_visit) %>%
      select(any_of(x = c("from_uid", "from_species", "to_uid", "to_species", "visit")))
    
    if(nrow(contact) == 0) {
      
      contact <- tibble(from_uid = reference_id,
                        from_species = reference_species,
                        to_uid = NA,
                        to_species = NA,
                        visit = reference_visit)
      
    } else {
      
      contact <- tibble(contact) %>%
        select(-geometry)
      
    }
    
  })
  
  # Combine the assemblages
  combined_assemblages <- bind_rows(assemblages, .id = "assemblage") %>%
    group_by(assemblage) %>%
    mutate(assemblage = as.numeric(assemblage))
  
  # Nodes
  nodes = unique(c(combined_assemblages$from_uid, combined_assemblages %>%
                     drop_na(to_uid) %>%
                     pull(to_uid)))
  
  vertices <- tibble(node = nodes) %>%
    left_join(tibble(rodent_data) %>%
                select(rodent_uid,
                       species = initial_species_id,
                       interpretation,
                       village,
                       visit,
                       simple_habitat,
                       month_set,
                       year_set,
                       sex,
                       age_group),
              by = c("node" = "rodent_uid")) %>%
    mutate(Species = str_to_sentence(str_replace_all(species, "_", " ")),
           Landuse = str_to_sentence(str_remove_all(simple_habitat, "secondary_")),
           Village = str_to_sentence(village),
           Sex = str_to_sentence(sex),
           Age = str_to_sentence(age_group)) %>%
    rename(ELISA = interpretation,
           Month = month_set,
           Year = year_set,
           Visit = visit) %>%
    select(node, Species, Village, Landuse, ELISA, Sex, Visit, Age, Month, Year)
  
  # Edgelist
  edgelist = combined_assemblages %>%
    ungroup() %>%
    drop_na(to_uid) %>%
    select(from_uid, to_uid)
  
  species_graph <- graph_from_data_frame(d = edgelist, vertices = vertices, directed = FALSE)
  
  simple_graph <- igraph::simplify(species_graph, remove.multiple = TRUE, remove.loops = TRUE)
  
  return(list(nodelist = vertices,
              edgelist = edgelist,
              graph = simple_graph,
              full_graph = species_graph))
  
}

assemblages_15 <- produce_assemblages(unique_rodents, distance = 15)
assemblages_30 <- produce_assemblages(unique_rodents, distance = 30)
assemblages_50 <- produce_assemblages(unique_rodents, distance = 50)

# Rodent species with less than 10 individuals will be grouped as Other
species_n <- tibble(species = get.vertex.attribute(assemblages_15$graph)$Species) %>%
  group_by(species) %>%
  summarise(N = n())

other_species <- species_n$species[species_n$N <= 10]

# Add further vertex attributes
# Species grouping
vertex_attr(assemblages_15$graph, name = "Species", index = V(assemblages_15$graph)$Species %in% other_species) <- "Other"
vertex_attr(assemblages_15$full_graph, name = "Species", index = V(assemblages_15$graph)$Species %in% other_species) <- "Other"
vertex_attr(assemblages_30$graph, name = "Species", index = V(assemblages_30$graph)$Species %in% other_species) <- "Other"
vertex_attr(assemblages_30$full_graph, name = "Species", index = V(assemblages_30$graph)$Species %in% other_species) <- "Other"
vertex_attr(assemblages_50$graph, name = "Species", index = V(assemblages_50$graph)$Species %in% other_species) <- "Other"
vertex_attr(assemblages_50$full_graph, name = "Species", index = V(assemblages_50$graph)$Species %in% other_species) <- "Other"

# Distance of contact definition
vertex_attr(assemblages_15$graph, name = "Contact radius") <- "15m"
vertex_attr(assemblages_30$graph, name = "Contact radius") <- "30m"
vertex_attr(assemblages_50$graph, name = "Contact radius") <- "50m"

# Village classification as Urban or Rural
vertex_attr(assemblages_15$graph, name = "Village locations", index = V(assemblages_15$graph)$"Village" == "Lambayama") <- "Urban"
vertex_attr(assemblages_15$graph, name = "Village locations", index = V(assemblages_15$graph)$"Village" != "Lambayama") <- "Rural"
vertex_attr(assemblages_30$graph, name = "Village locations", index = V(assemblages_30$graph)$"Village" == "Lambayama") <- "Urban"
vertex_attr(assemblages_30$graph, name = "Village locations", index = V(assemblages_30$graph)$"Village" != "Lambayama") <- "Rural"
vertex_attr(assemblages_50$graph, name = "Village locations", index = V(assemblages_50$graph)$"Village" == "Lambayama") <- "Urban"
vertex_attr(assemblages_50$graph, name = "Village locations", index = V(assemblages_50$graph)$"Village" != "Lambayama") <- "Rural"

# Season classification as wet or dry
vertex_attr(assemblages_15$graph, name = "Season", index = V(assemblages_15$graph)$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
vertex_attr(assemblages_15$graph, name = "Season", index = V(assemblages_15$graph)$"Visit" %in% c("3", "4")) <- "Wet"
vertex_attr(assemblages_30$graph, name = "Season", index = V(assemblages_30$graph)$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
vertex_attr(assemblages_30$graph, name = "Season", index = V(assemblages_30$graph)$"Visit" %in% c("3", "4")) <- "Wet"
vertex_attr(assemblages_50$graph, name = "Season", index = V(assemblages_50$graph)$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
vertex_attr(assemblages_50$graph, name = "Season", index = V(assemblages_50$graph)$"Visit" %in% c("3", "4")) <- "Wet"

# Set equivocal ELISA as negative
vertex_attr(assemblages_15$graph, name = "ELISA", index = V(assemblages_15$graph)$"ELISA" == "Equivocal") <- "Negative"
vertex_attr(assemblages_30$graph, name = "ELISA", index = V(assemblages_30$graph)$"ELISA" == "Equivocal") <- "Negative"
vertex_attr(assemblages_50$graph, name = "ELISA", index = V(assemblages_50$graph)$"ELISA" == "Equivocal") <- "Negative"

# Consistent colours for species
species_palette_df <- igraph::as_data_frame(assemblages_15$graph, what = "vertices") %>%
  select(Species) %>%
  distinct(Species) %>%
  arrange(Species) %>%
  mutate(colour = brewer.pal(length(.$Species), "Dark2"))

species_palette <- c(species_palette_df$colour)
names(species_palette) <- c(species_palette_df$Species)

# Combined graphs
graphs <- list(`15m` = assemblages_15$graph,
               `30m` = assemblages_30$graph,
               `50m` = assemblages_50$graph)

combined_graphs <- lapply(graphs, function(x) {
  ggraph(x, layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 2) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0("Combined landuse: ", unique(V(x)$`Contact radius`))) +
    theme_graph(title_size = 16, base_family = "sans")
  })

plot_combined_graphs <- plot_grid(plotlist = list(combined_graphs[[1]] +
                                                   theme(legend.position = "none"),
                                                  combined_graphs[[2]] +
                                                   theme(legend.position = "none"),
                                                  combined_graphs[[3]] +
                                                   theme(legend.position = "none"),
                                                 legend))

save_plot(plot = plot_combined_graphs, filename = here("output", "combined_graph.pdf"), base_height = 14, base_width = 16)

# Landuse graphs
village_graphs <- lapply(graphs, function(x) {
  subgraph(x, vids = V(x)$Landuse == "Village") %>%
  ggraph(layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Species), size = 3) +
  scale_colour_manual(values = species_palette) +
  labs(colour = "Species",
       title = paste0("Village: ", unique(V(x)$`Contact radius`))) +
  theme_graph(title_size = 16, base_family = "sans")
  })

agriculture_graphs <- lapply(graphs, function(x) {
  subgraph(x, vids = V(x)$Landuse == "Agriculture") %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0("Agriculture: ", unique(V(x)$`Contact radius`))) +
    theme_graph(title_size = 16, base_family = "sans")
})

forest_graphs <- lapply(graphs, function(x) {
  subgraph(x, vids = V(x)$Landuse == "Forest") %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0("Forest: ", unique(V(x)$`Contact radius`))) +
    theme_graph(title_size = 16, base_family = "sans")
})

legend <- get_legend(village_graphs[[1]])

plot_village_graphs <- plot_grid(plotlist = list(village_graphs[[1]] +
                                                   theme(legend.position = "none"),
                                                 village_graphs[[2]] +
                                                   theme(legend.position = "none"),
                                                 village_graphs[[3]] +
                                                   theme(legend.position = "none"),
                                                 legend))

plot_agriculture_graphs <- plot_grid(plotlist = list(agriculture_graphs[[1]] +
                                                   theme(legend.position = "none"),
                                                   agriculture_graphs[[2]] +
                                                   theme(legend.position = "none"),
                                                   agriculture_graphs[[3]] +
                                                   theme(legend.position = "none"),
                                                 legend))

plot_forest_graphs <- plot_grid(plotlist = list(forest_graphs[[1]] +
                                                   theme(legend.position = "none"),
                                                forest_graphs[[2]] +
                                                   theme(legend.position = "none"),
                                                forest_graphs[[3]] +
                                                   theme(legend.position = "none"),
                                                 legend))

save_plot(plot = plot_village_graphs, filename = here("output", "village_graph.pdf"), base_height = 14, base_width = 16)
save_plot(plot = plot_agriculture_graphs, filename = here("output", "agriculture_graph.pdf"), base_height = 14, base_width = 16)
save_plot(plot = plot_forest_graphs, filename = here("output", "forest_graph.pdf"), base_height = 14, base_width = 16)

# Describing contact networks --------------------------------------------------------

rodent_network <- graphs
names(rodent_network) <- c("observed_15", "observed_30", "observed_50")

write_rds(rodent_network, here("data", "rodent_network.rds"))

# Landuse network metrics
rodent_networks <- lapply(rodent_network, function(x) {asNetwork(x)})

options(max.print = 50)
summary(rodent_networks[[1]])
summary(rodent_networks[[2]])
summary(rodent_networks[[3]])

rodent_net_descriptives <- lapply(rodent_networks, function(x) {
  list("mean_degree" = mean(degree(x, gmode = "graph")),
       "sd_degree" = sd(degree(x, gmode = "graph")),
       "degree_table" = table(degree(x, gmode = "graph")),
       "triad_table" = triad.census(x, mode = "graph"),
       "density" = graph.density(asIgraph(x)))
})

# Produce a random graph of similar properties for comparison
# Random graph of similar size and density
random_rodent_networks <- lapply(rodent_networks, function(x) {
  as.network(rgraph(412, m = 1, tprob = graph.density(asIgraph(x)), mode = "graph"), directed = FALSE)
})
names(random_rodent_networks) <- c("random_15", "random_30", "random_50")
random_rodent_networks$random_15 %v% "Contact.radius" <- "15m"
random_rodent_networks$random_30 %v% "Contact.radius" <- "30m"
random_rodent_networks$random_50 %v% "Contact.radius" <- "50m"

# Assess the degree distribution of these random networks to the observed
degree_plots <- list()

for(i in 1:length(rodent_networks)) {
  degree_plots[[i]] <- plot_grid(plotlist = list(tibble("degree" = degree(random_rodent_networks[[i]], gmode = "graph")) %>%
                                              ggplot() + 
                                              geom_bar(aes(x = degree)) +
                                              labs(title = paste0("Degree: ",  unique(random_rodent_networks[[i]] %v% "Contact.radius")," (Random)")) +
                                              theme_bw(),
                                            tibble("degree" = degree(rodent_networks[[i]], gmode = "graph")) %>%
                                              ggplot() +
                                              geom_bar(aes(x = degree)) +
                                              labs(title = paste0("Degree: ", unique(rodent_networks[[i]] %v% "Contact.radius"), " (Rodent)")) +
                                              theme_bw()))
                                            }

# The clustering of these networks can be seen using mixing matrices
options(max.print = 200)

mixing_matrices <- lapply(rodent_networks, function(x) {
  list("Species matrix" = mixingmatrix(x, "Species"),
       "Village location matrix" = mixingmatrix(x, "Village.locations"),
       "Landuse matrix" = mixingmatrix(x, "Landuse"),
       "Serostatus matrix" = mixingmatrix(x, "ELISA"))
})

species_degree <- lapply(rodent_networks, function(x) tibble("Species" = x %v% "Species",
                                                      "Degree" = degree(x),
                                                      "Distance" = x %v% "Contact.radius")) %>%
  bind_rows() %>%
  group_by(Species, Distance) %>%
  summarise(Degree = mean(Degree)) %>%
  ungroup() %>%
  arrange(-Degree) %>%
  mutate(Species = fct_inorder(Species)) %>%
  ggplot() +
  geom_point(aes(x = Species, y = Degree, colour = Species)) +
  labs(title = "Mean degree by species") +
  scale_colour_manual(values = species_palette) +
  coord_flip() +
  facet_wrap(~ Distance) +
  theme_bw()


# Build models ------------------------------------------------------------
source(here("R", "modified_functions.R"))

# Random graph models
null_models <- lapply(random_rodent_networks, function(x) ergm(x ~ edges))
summary_null_models <- lapply(null_models, function(x) summary.ergm.david(x))

# Main effects only models
rodent_models <- list()
rodent_models$main_effects <- lapply(rodent_networks, function(x) {
  ergm(x ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")))
})
rodent_models_summary <- list()
rodent_models_summary$main_effects <- lapply(rodent_models$main_effects, function(x) {summary.ergm.david(x)})

# Homophily model for main effects term
rodent_models$homophily <- lapply(rodent_networks, function(x) {
  ergm(x ~ edges + nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
         + nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE))
})
rodent_models_summary$homophily <- lapply(rodent_models$homophily, function(x) {summary.ergm.david(x)})

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



# CEF ---------------------------------------------------------------------
# 
# cefmodel_50m <- ergm(rodent_networks$observed_50 ~ edges +
#                            nodefactor("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides")) +
#                            nodematch("Species", levels = c("Crocidura spp", "Mastomys natalensis", "Praomys spp", "Lophuromys sikapusi", "Mus musculus", "Rattus rattus", "Mus minutoides"), diff = TRUE) +
#                            nodematch("Season", levels = "Dry") +
#                            nodematch("Village.locations", levels =  "Rural", diff = TRUE) +
#                            nodematch("Landuse", levels = c("Agriculture", "Village"), diff = TRUE) +
#                            gwdegree(cutoff = 5) +
#                            gwesp(cutoff = 5) +
#                            gwdsp(cutoff = 5),
#                          eval.loglik = TRUE,
#                          verbose = TRUE)
# write_rds(cefmodel_50m, here("temp", "cef_species.rds"))
