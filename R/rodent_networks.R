# Load in the rodent data
unique_rodents <- read_rds(here("data", "unique_rodents.rds")) %>%
  drop_na(initial_species_id) %>%
  st_as_sf(crs = default_CRS) %>%
  st_transform(crs = SL_UTM)

produce_assemblages <- function(rodent_data = unique_rodents, distance = 50) {
  
  units(distance) <- "meters"
  
  individual_buffers <- rodent_data %>%
    mutate(rodent_id = row_number()) %>%
    st_buffer(dist = distance) %>%
    group_by(rodent_id) %>%
    group_split()
  
  comparator_rodents <- rodent_data %>%
    mutate(rodent_id = row_number())
  
  assemblages <- lapply(individual_buffers, function(x) {
    
    reference_visit <- unique(x$visit)
    
    reference_species <- unique(x$clean_names)
    
    reference_id <- unique(x$rodent_id)
    
    assemblage <- st_join(comparator_rodents %>%
                            filter(rodent_id != reference_id) %>%
                            rename(co_occurring_id = rodent_id), 
                          x %>%
                            select(rodent_id, geometry) %>%
                            rename(reference_id = rodent_id),
                          left = FALSE, join = st_within) %>%
      bind_rows(comparator_rodents %>%
                  filter(rodent_id == reference_id)) %>%
      mutate(reference_rodent = case_when(is.na(rodent_id) ~ FALSE,
                                          TRUE ~ TRUE),
             same_visit = case_when(visit == x$visit[!is.na(x$rodent_id)] ~ TRUE,
                                    TRUE ~ FALSE),
             rodent_id = coalesce(rodent_id, co_occurring_id)) %>%
      select(-co_occurring_id) %>%
      arrange(rodent_id)
    
  })
  
  combined_assemblages <- bind_rows(assemblages, .id = "assemblage") %>%
    group_by(assemblage) %>%
    mutate(assemblage = as.numeric(assemblage))
  
  summarise_assemblages <- combined_assemblages %>%
    tibble() %>%
    group_by(assemblage, initial_species_id, reference_rodent, same_visit) %>%
    summarise(n_individuals  = n()) %>%
    filter(!is.na(initial_species_id))
  
  same_visit <- summarise_assemblages %>%
    filter(same_visit == TRUE)
  
  summarise_assemblages_same_visit_plot <- same_visit %>%
    group_by(assemblage) %>%
    mutate(n_other_species = n()-1,
           clean_names = str_to_sentence(str_replace_all(initial_species_id, "_", " "))) %>%
    filter(reference_rodent == TRUE) %>%
    ggplot(aes(x = n_other_species)) + 
    geom_bar() +
    facet_wrap(~ clean_names) +
    theme_bw() +
    labs(x = "Co-located species (N)",
         y = element_blank(),
         title = "Same visit")
  
  all_visits <- summarise_assemblages
  
  summarise_assemblages_all_visits_plot <- all_visits %>%
    group_by(assemblage) %>%
    mutate(n_other_species = n()-1,
           clean_names = str_to_sentence(str_replace_all(initial_species_id, "_", " "))) %>%
    filter(reference_rodent == TRUE) %>%
    ggplot(aes(x = n_other_species)) + 
    geom_bar() +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    facet_wrap(~ clean_names) +
    theme_bw() +
    labs(x = "Co-located species (N)",
         y = element_blank(),
         title = "All visits")
  
  colocated_species_same_visit <- same_visit %>%
    group_by(assemblage) %>%
    mutate(reference_species = case_when(reference_rodent == TRUE ~ initial_species_id,
                                         TRUE ~ as.character(NA)),
           colocated_species = case_when(reference_rodent == FALSE ~ initial_species_id,
                                         TRUE ~ as.character(NA)))
  
  from_same_visit <- colocated_species_same_visit %>%
    filter(reference_rodent == TRUE) %>%
    select(-colocated_species) %>%
    rename(from_species = reference_species)
  
  to_same_visit <- colocated_species_same_visit %>%
    filter(reference_rodent == FALSE) %>%
    select(-reference_species) %>%
    rename(to_species = colocated_species)
  
  edgelist_same_visit <- left_join(to_same_visit, from_same_visit %>%
                                     select(assemblage, from_species), by = c("assemblage")) %>%
    select(assemblage, from_species, to_species, n_individuals)
  
  prop_same_visit <- edgelist_same_visit %>%
    group_by(from_species, to_species) %>%
    summarise(prop = sum(n_individuals)) %>%
    group_by(from_species) %>%
    mutate(total = sum(prop)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(prop = round(prop/total, 2))
  
  order_species <- c(prop_same_visit %>% 
                       arrange(-total) %>%
                       distinct(from_species) %>%
                       pull(from_species),
                     "gerbilliscus_spp")
  
  order_species_N <- tibble(from_species = order_species) %>%
    left_join(prop_same_visit %>%
                distinct(from_species, total)) %>%
    mutate(total = replace_na(total, replace = 0)) %>%
    mutate(from_species = paste0(str_replace_all(str_to_sentence(from_species), "_", "\n"), " (", total, ")")) %>%
    pull(from_species)
  
  prop_same_visit <- prop_same_visit %>%
    mutate(from_species = paste0(str_replace_all(str_to_sentence(from_species), "_", "\n"), " (", total, ")"),
           from_species = factor(from_species, levels = order_species_N),
           to_species = factor(to_species, levels = order_species, labels = str_replace_all(str_to_sentence(order_species), "_", "\n")))
  
  colocated_species_same_visit_plot <- ggplot(prop_same_visit) +
    geom_tile(aes(y = fct_rev(from_species), x = to_species, fill = prop, width = 0.95, height = 0.95)) +
    geom_label(aes(y = fct_rev(from_species), x = to_species, label = prop)) +
    scale_fill_viridis_c() +
    theme_bw() +
    labs(y = "Reference species (Edges)",
         x = "Co-located species",
         fill = "Proportion",
         title = "Same visit")
  
  
  # Graph based interpretation
  species_graph <- graph_from_data_frame(edgelist_same_visit %>%
                                           group_by(from_species, to_species) %>%
                                           summarise(prop = sum(n_individuals)) %>%
                                           group_by(from_species) %>%
                                           mutate(total = sum(prop)) %>%
                                           ungroup() %>%
                                           rowwise() %>%
                                           mutate(prop = round(prop/total, 2)),
                                         directed = TRUE)
  
  landuse_for_edgelist <- combined_assemblages %>%
    tibble() %>%
    filter(same_visit == TRUE) %>%
    distinct(assemblage, initial_species_id, reference_rodent, simple_habitat)
  
  edgelist_landuse <- edgelist_same_visit %>%
    left_join(landuse_for_edgelist %>%
                filter(reference_rodent == TRUE) %>%
                rename(from_species = initial_species_id,
                       from_landuse = simple_habitat) %>%
                distinct(assemblage, from_species, from_landuse),
              by = c("assemblage", "from_species")) %>%
    left_join(landuse_for_edgelist %>%
                filter(reference_rodent == FALSE) %>%
                rename(to_species = initial_species_id,
                       to_landuse = simple_habitat) %>%
                distinct(assemblage, to_species, to_landuse),
              by = c("assemblage", "to_species")) %>%
    filter(from_landuse == to_landuse)
  
  species_graph_village <- graph_from_data_frame(edgelist_landuse %>%
                                                   filter(from_landuse == "village") %>%
                                                   group_by(from_species, to_species) %>%
                                                   summarise(prop = sum(n_individuals)) %>%
                                                   group_by(from_species) %>%
                                                   mutate(total = sum(prop)) %>%
                                                   ungroup() %>%
                                                   rowwise() %>%
                                                   mutate(prop = round(prop/total, 2)),
                                                 directed = TRUE)
  
  species_graph_agriculture <- graph_from_data_frame(edgelist_landuse %>%
                                                       filter(from_landuse == "agriculture") %>%
                                                       group_by(from_species, to_species) %>%
                                                       summarise(prop = sum(n_individuals)) %>%
                                                       group_by(from_species) %>%
                                                       mutate(total = sum(prop)) %>%
                                                       ungroup() %>%
                                                       rowwise() %>%
                                                       mutate(prop = round(prop/total, 2)),
                                                     directed = TRUE)
  
  species_graph_forest <- graph_from_data_frame(edgelist_landuse %>%
                                                  filter(from_landuse == "secondary_forest") %>%
                                                  group_by(from_species, to_species) %>%
                                                  summarise(prop = sum(n_individuals)) %>%
                                                  group_by(from_species) %>%
                                                  mutate(total = sum(prop)) %>%
                                                  ungroup() %>%
                                                  rowwise() %>%
                                                  mutate(prop = round(prop/total, 2)),
                                                directed = TRUE)
  
  
  return(list(assemblages = combined_assemblages,
              summary_plots = list(same_visit = summarise_assemblages_same_visit_plot,
                                   all_visits = summarise_assemblages_all_visits_plot),
              co_location_plots = list(same_visit = colocated_species_same_visit_plot),
              graphs = list(same_visit_graph = species_graph,
                            species_graph_village = species_graph_village,
                            species_graph_agriculture = species_graph_agriculture,
                            species_graph_forest = species_graph_forest)))
  
}

assemblages <- produce_assemblages(unique_rodents, distance = 50)


ggraph(assemblages$graphs$species_graph_village, layout = "kk") +
  geom_edge_fan0(aes(colour = prop, width = prop)) +
  geom_edge_loop(aes(colour = prop,  width = prop, span = 90), alpha = 0.4) +
  geom_node_point() +
  geom_node_label(aes(label = name)) +
  scale_edge_color_viridis(option = "D", direction = -1, limits = c(0, 1), guide = "none") +
  scale_edge_width(range = c(0, 3)) +
  theme_graph(title_size = 16, base_family = "sans") +
  labs(title = "Village")

ggraph(assemblages$graphs$species_graph_agriculture, layout = "kk") +
  geom_edge_fan0(aes(colour = prop, width = prop)) +
  geom_edge_loop(aes(colour = prop,  width = prop, span = 90), alpha = 0.4) +
  geom_node_point() +
  geom_node_label(aes(label = name)) +
  scale_edge_color_viridis(option = "D", direction = -1, limits = c(0, 1), guide = "none") +
  scale_edge_width(range = c(0, 3)) +
  theme_graph(title_size = 16, base_family = "sans") +
  labs(title = "Agriculture")

ggraph(assemblages$graphs$species_graph_forest, layout = "kk") +
  geom_edge_fan0(aes(colour = prop, width = prop)) +
  geom_edge_loop(aes(colour = prop,  width = prop, span = 90), alpha = 0.4) +
  geom_node_point() +
  geom_node_label(aes(label = name)) +
  scale_edge_color_viridis(option = "D", direction = -1, limits = c(0, 1), guide = "none") +
  scale_edge_width(range = c(0, 3)) +
  theme_graph(title_size = 16, base_family = "sans") +
  labs(title = "Forest")
