rodents <- unique_rodents %>%
  left_join(trap_data) %>%
  rename(grid = grid_number)

grids_polygon <- trap_data %>%
  st_as_sf(coords = c("trap_easting", "trap_northing"),
           crs = SL_UTM) %>%
  group_by(village, landuse, visit, grid_number) %>%
  summarise() %>%
  st_convex_hull() %>%
  ungroup() %>%
  mutate(area = st_area(.))
  
n_in_grid <- tibble(grids_polygon) %>%
  select(-geometry) %>%
  left_join(rodents %>%
              group_by(village, landuse, visit, grid) %>%
              summarise(n = n()),
            by = c("village", "landuse", "visit", "grid_number" = "grid")) %>%
  mutate(density_rodent = n/area)

density_rodents <- n_in_grid %>%
  drop_na(density_rodent) %>%
  group_by(landuse) %>%
  summarise(density_rodent = median(density_rodent))

detected_grids <- rodents %>%
  st_as_sf(coords = c("trap_easting", "trap_northing"),
           crs = SL_UTM) %>%
  st_buffer(dist = 30) %>%
  group_by(village, visit, grid) %>%
  group_split()

grid_networks <- function(df = detected_grid) {
  
  if(nrow(df) == 1) {
    
    edges <- tibble(rodent_uid = df$rodent_uid,
                    village = df$village,
                    visit = df$visit,
                    grid_number = df$grid,
                    landuse = df$landuse,
                    degree = 0)
    
  } else {
    
    edges <- st_overlaps(df)
    
    edges <- tibble(rodent_uid = df$rodent_uid,
                    village = df$village,
                    visit = df$visit,
                    grid_number = df$grid,
                    landuse = df$landuse,
                    degree = lapply(edges, function(x) length(x)) %>% unlist())
    
  }
}
    
    
rodent_contacts <- lapply(detected_grids, function(x) grid_networks(x)) %>%
  bind_rows()

left_join(rodent_contacts, n_in_grid) %>%
  ggplot() +
  geom_point(aes(x = density_rodent, y = degree, colour = landuse))

density_degree <- left_join(rodent_contacts, n_in_grid)

cor.test(density_degree$density_rodent, density_degree$degree)
cor.test(density_degree$n, density_degree$degree)

summarise_degree <- left_join(rodent_contacts, n_in_grid)  %>%
  group_by(landuse) %>%
  summarise(median_degree = median(degree),
            IQR_degree = IQR(degree))


# Calculate betweenness on a node level -----------------------------------
assemblages_grid <- lapply(detected_grids, function(y) {
  
  # Allocate a village
  reference_village <- unique(y$village)
  
  # Allocate a visit to a rodent
  reference_visit <- unique(y$visit)
  
  # Allocate a grid
  reference_grid <- unique(y$grid)
  
  individual_rodents <- y %>%
    group_by(rodent_uid) %>%
    group_split()
  
  lapply(individual_rodents, function(z) {
    
    reference_id <- unique(z$rodent_uid)
    reference_species <- unique(z$species)
    
    contact <- st_filter(y %>%
                filter(!rodent_uid %in% reference_id) %>%
                rename(to_id = rodent_uid,
                       to_species = species),
              z,
              .predicate = st_overlaps) 
    
    if(nrow(contact) == 0) {
    
    contact <- tibble(from_id = reference_id,
                      from_species = reference_species,
                      to_id = NA,
                      to_species = NA,
                      village = reference_village,
                      visit = reference_visit,
                      grid = reference_grid,
                      landuse = unique(z$landuse),
                      interpretation = NA)
    } else {
    
    contact <- tibble(contact) %>%
      select(-geometry) %>%
      mutate(from_id = reference_id,
             from_species = reference_species) %>%
      tibble() %>%
      select(from_id, from_species, to_id, to_species, village, visit, grid, landuse, interpretation)
    
  }
  }) %>%
    bind_rows()
  
  })

# Combine the assemblages
combined_assemblages_grid <- bind_rows(assemblages_grid, .id = "assemblage") %>%
  group_by(assemblage) %>%
  mutate(assemblage = as.numeric(assemblage)) %>%
  group_split()

grid_graphs <- lapply(combined_assemblages_grid, function(x) {
  
  if(nrow(x) == 1) {
    
    nodes = x$from_id
    
  } else {
    # Nodes
    nodes = unique(c(x$from_id,
                     x %>%
                       drop_na(to_id) %>%
                       pull(to_id)))
  }
  
  vertices <- tibble(node = nodes) %>%
    left_join(x,
              by = c("node" = "from_id")) %>%
    distinct(node, Species = from_species, Village = village, Visit = visit, Grid = grid, Landuse = landuse, assemblage)
  
  # Edgelist
  edgelist <- x %>%
    ungroup() %>%
    drop_na(to_id) %>%
    select(from_id, to_id) %>%
    mutate(value = 1)
  
  species_graph <- graph_from_data_frame(d = edgelist, vertices = vertices, directed = FALSE)
  
  x %>%
    distinct(to_id, interpretation)
  
  V(species_graph)$ELISA_status <- tibble(vertex_name = V(species_graph)$name) %>%
    left_join(x %>%
                distinct(to_id, interpretation),
              by = c("vertex_name" = "to_id")) %>%
    pull(interpretation)
  
  grid_descriptives <- tibble(village = unique(x$village),
                              visit = unique(x$visit),
                              grid = unique(x$grid),
                              landuse = unique(x$landuse),
                              rodent_id = V(species_graph)$name,
                              species = V(species_graph)$Species,
                              degree = igraph::degree(species_graph),
                              contact_positive = V(species_graph)$ELISA_status,
                              betweenness = estimate_betweenness(species_graph, cutoff = -1),
                              density = igraph::edge_density(species_graph))
  
  return(grid_descriptives)
  
}) %>%
  bind_rows()

grid_graphs %>% group_by(landuse) %>% summarise(richness = unique(length(unique(species))),
                                                n = n(),
                                                mean_betweenness = mean(betweenness),
                                                mean_density = mean(density, na.rm = TRUE),
                                                mean_degree = mean(degree, na.rm = TRUE),
                                                sd_degree = sd(degree, na.rm = TRUE),
                                                min_degree = min(degree, na.rm = TRUE),
                                                max_degree = max(degree, na.rm = TRUE))

grid_graphs %>% group_by(species) %>% summarise(n = n(),
                                                mean_degree = mean(degree, na.rm = TRUE),
                                                sd_degree = sd(degree, na.rm = TRUE),
                                                min_degree = min(degree, na.rm = TRUE),
                                                max_degree = max(degree, na.rm = TRUE)) %>% 
  arrange(-n)
