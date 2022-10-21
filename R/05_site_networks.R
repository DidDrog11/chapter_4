# Load in the rodent data and split by site
site_rodents <- read_rds(here("data", "unique_rodents.rds")) %>%
  drop_na(initial_species_id) %>%
  drop_na(grid_number) %>%
  st_as_sf(crs = default_CRS) %>%
  st_transform(crs = SL_UTM) %>%
  mutate(grid_number = case_when(grid_number %in% c("6", "7") ~ "6",
                          TRUE ~ as.character(grid_number))) %>%
  group_by(village, visit, grid_number) %>%
  group_split()

produce_assemblages_site <- function(rodent_data = site_rodents, distance = 15) {
  
  # Define distance as meters
  units(distance) <- "meters"
  
  # Produce a buffer around each individual to create their range
  individual_buffers <- lapply(rodent_data, function(x) { x %>%
      select(rodent_uid, initial_species_id, village, visit, grid_number) %>%
      st_buffer(dist = distance) %>%
      group_by(rodent_uid) %>%
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
      reference_species <- unique(x$initial_species_id)
      
      # Allocate an ID to use for reference
      reference_id <- unique(x$rodent_uid)
      
      # Join the reference rodent to others that were detected in their range
      contact <- st_filter(bind_rows(rodent_data) %>%
                             filter(!rodent_uid %in% reference_id &
                                      village %in% reference_village &
                                      visit %in% reference_visit &
                                      grid_number %in% reference_grid_number) %>%
                             rename(to_uid = rodent_uid,
                                    to_species = initial_species_id), 
                           x %>%
                             select(rodent_uid, geometry),
                           .predicate = st_within) %>%
        mutate(from_uid = reference_id,
               from_species = reference_species) %>%
        select(any_of(x = c("from_uid", "from_species", "to_uid", "to_species", "village", "visit", "grid_number")))
      
      if(nrow(contact) == 0) {
        
        contact <- tibble(from_uid = reference_id,
                          from_species = reference_species,
                          to_uid = NA,
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
  nodes = lapply(combined_assemblages, function(x) {unique(c(x$from_uid, x %>%
                     drop_na(to_uid) %>%
                     pull(to_uid)))
    })
  
  vertices <- lapply(nodes, function(x) {
    tibble(node = x) %>%
    left_join(tibble(bind_rows(rodent_data)) %>%
                select(rodent_uid,
                       species = initial_species_id,
                       village,
                       visit,
                       grid_number,
                       interpretation,
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
    select(node, Species, Village, Visit, Grid_number = "grid_number", Landuse, ELISA, Sex, Visit, Age, Month, Year)
    })
  
  # Edgelist
  edgelist = lapply(combined_assemblages, function(x) {x %>%
    ungroup() %>%
    drop_na(to_uid) %>%
    select(from_uid, to_uid) 
    })
  
  species_graph <- mapply(x = edgelist, y = vertices, function(x, y) {graph_from_data_frame(d = x, vertices = y, directed = FALSE)
    })
  
  simple_graph <- lapply(species_graph, function(x) {igraph::simplify(x, remove.multiple = TRUE, remove.loops = TRUE)})
  
  return(list(nodelist = vertices,
              edgelist = edgelist,
              graph = simple_graph,
              full_graph = species_graph))
  
}

assemblages_15 <- produce_assemblages_site(site_rodents, distance = 15)
assemblages_30 <- produce_assemblages_site(site_rodents, distance = 30)
assemblages_50 <- produce_assemblages_site(site_rodents, distance = 50)

# Rodent species with less than 10 individuals will be grouped as Other
species_n <- bind_rows(site_rodents) %>%
  group_by(initial_species_id) %>%
  summarise(N = n()) %>%
  mutate(species = str_to_sentence(str_replace_all(initial_species_id, "_", " ")))

other_species <- species_n$species[species_n$N <= 10]

# Add further vertex attributes

for(i in 1:length(assemblages_15$graph)) {
  # Species grouping
  vertex_attr(assemblages_15$graph[[i]], name = "Species", index = V(assemblages_15$graph[[i]])$Species %in% other_species) <- "Other"
  vertex_attr(assemblages_15$full_graph[[i]], name = "Species", index = V(assemblages_15$graph[[i]])$Species %in% other_species) <- "Other"
  vertex_attr(assemblages_30$graph[[i]], name = "Species", index = V(assemblages_30$graph[[i]])$Species %in% other_species) <- "Other"
  vertex_attr(assemblages_30$full_graph[[i]], name = "Species", index = V(assemblages_30$graph[[i]])$Species %in% other_species) <- "Other"
  vertex_attr(assemblages_50$graph[[i]], name = "Species", index = V(assemblages_50$graph[[i]])$Species %in% other_species) <- "Other"
  vertex_attr(assemblages_50$full_graph[[i]], name = "Species", index = V(assemblages_50$graph[[i]])$Species %in% other_species) <- "Other"
  
  # Distance of contact definition
  vertex_attr(assemblages_15$graph[[i]], name = "Contact radius") <- "15m"
  vertex_attr(assemblages_30$graph[[i]], name = "Contact radius") <- "30m"
  vertex_attr(assemblages_50$graph[[i]], name = "Contact radius") <- "50m"
  
  # Village classification as Urban or Rural
  vertex_attr(assemblages_15$graph[[i]], name = "Village locations", index = V(assemblages_15$graph[[i]])$"Village" == "Lambayama") <- "Urban"
  vertex_attr(assemblages_15$graph[[i]], name = "Village locations", index = V(assemblages_15$graph[[i]])$"Village" != "Lambayama") <- "Rural"
  vertex_attr(assemblages_30$graph[[i]], name = "Village locations", index = V(assemblages_30$graph[[i]])$"Village" == "Lambayama") <- "Urban"
  vertex_attr(assemblages_30$graph[[i]], name = "Village locations", index = V(assemblages_30$graph[[i]])$"Village" != "Lambayama") <- "Rural"
  vertex_attr(assemblages_50$graph[[i]], name = "Village locations", index = V(assemblages_50$graph[[i]])$"Village" == "Lambayama") <- "Urban"
  vertex_attr(assemblages_50$graph[[i]], name = "Village locations", index = V(assemblages_50$graph[[i]])$"Village" != "Lambayama") <- "Rural"
  
  # Season classification as wet or dry
  vertex_attr(assemblages_15$graph[[i]], name = "Season", index = V(assemblages_15$graph[[i]])$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
  vertex_attr(assemblages_15$graph[[i]], name = "Season", index = V(assemblages_15$graph[[i]])$"Visit" %in% c("3", "4")) <- "Wet"
  vertex_attr(assemblages_30$graph[[i]], name = "Season", index = V(assemblages_30$graph[[i]])$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
  vertex_attr(assemblages_30$graph[[i]], name = "Season", index = V(assemblages_30$graph[[i]])$"Visit" %in% c("3", "4")) <- "Wet"
  vertex_attr(assemblages_50$graph[[i]], name = "Season", index = V(assemblages_50$graph[[i]])$"Visit" %in% c("1", "2", "5", "6")) <- "Dry"
  vertex_attr(assemblages_50$graph[[i]], name = "Season", index = V(assemblages_50$graph[[i]])$"Visit" %in% c("3", "4")) <- "Wet"
  
  # Set equivocal ELISA as negative
  vertex_attr(assemblages_15$graph[[i]], name = "ELISA", index = V(assemblages_15$graph[[i]])$"ELISA" == "Equivocal") <- "Negative"
  vertex_attr(assemblages_30$graph[[i]], name = "ELISA", index = V(assemblages_30$graph[[i]])$"ELISA" == "Equivocal") <- "Negative"
  vertex_attr(assemblages_50$graph[[i]], name = "ELISA", index = V(assemblages_50$graph[[i]])$"ELISA" == "Equivocal") <- "Negative"
  
}

# Consistent colours for species
species_palette_df <- as_tibble(bind_graphs(assemblages_15$graph)) %>%
  select(Species) %>%
  distinct(Species) %>%
  arrange(Species) %>%
  mutate(colour = brewer.pal(length(.$Species), "Dark2"))

species_palette <- c(species_palette_df$colour)
names(species_palette) <- c(species_palette_df$Species)

plot_graph <- list()
plot_graph$assemblages_15 <- list()

plot_graph$assemblages_15 <- lapply(assemblages_15$graph, function(x) {
  as_tbl_graph(x) %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0(unique(V(x)$Village), ", Grid ", unique(V(x)$Grid_number), ", Visit ", unique(V(x)$Visit)),
         caption = "15m") +
    theme_graph(title_size = 16, base_family = "sans")
})

plot_graph$assemblages_30 <- lapply(assemblages_30$graph, function(x) {
  as_tbl_graph(x) %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0(unique(V(x)$Village), ", Grid ", unique(V(x)$Grid_number), ", Visit ", unique(V(x)$Visit)),
         caption = "30m") +
    theme_graph(title_size = 16, base_family = "sans")
})

plot_graph$assemblages_50 <- lapply(assemblages_50$graph, function(x) {
  as_tbl_graph(x) %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette) +
    labs(colour = "Species",
         title = paste0(unique(V(x)$Village), ", Grid ", unique(V(x)$Grid_number), ", Visit ", unique(V(x)$Visit)),
         caption = "50m") +
    theme_graph(title_size = 16, base_family = "sans")
})

# Describing contact networks --------------------------------------------------------

rodent_network <- list(observed_15 = assemblages_15$graph,
                       observed_30 = assemblages_30$graph,
                       observed_50 = assemblages_50$graph)

write_rds(rodent_network, here("data", "rodent_network_site.rds"))

# Landuse network metrics
rodent_networks <- lapply(rodent_network, function(x) lapply(x, function(y) {asNetwork(y)}))

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