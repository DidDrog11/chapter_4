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

produce_assemblages <- function(rodent_data = unique_rodents, distance = 15) {
  
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

assemblages <- produce_assemblages(unique_rodents, distance = 15)

# Consistent colours for species
species_palette_df <- assemblages$assemblages %>%
  tibble() %>%
  select(initial_species_id) %>%
  mutate(species = str_to_sentence(str_replace_all(initial_species_id, "_", " "))) %>%
  distinct(species) %>%
  mutate(colour = colorRampPalette(brewer.pal(9, 
                                              "Set1"))(13))
species_palette <- c(species_palette_df$colour)
names(species_palette) <- c(species_palette_df$species)

# Landuse graphs

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

# Contact networks --------------------------------------------------------

network <- assemblages$assemblages %>%
  tibble() %>%
  filter(same_visit == TRUE) %>%
  select(village, visit, grid_number, simple_habitat, interpretation, reference_id, rodent_id, initial_species_id, reference_rodent) %>%
  mutate(reference_id = coalesce(reference_id, rodent_id)) %>%
  group_by(reference_id) %>%
  mutate(reference_name = case_when(reference_rodent == TRUE ~ initial_species_id,
                                    reference_rodent == FALSE ~ initial_species_id[reference_rodent == TRUE]),
         interpretation = case_when(reference_rodent == TRUE ~ interpretation,
                                    reference_rodent == FALSE ~ interpretation[reference_rodent == TRUE])) %>%
  group_by(simple_habitat) %>%
  group_split()

names(network) <- lapply(network, function(x) {unique(x$simple_habitat)})

edgelist <- lapply(network, function(x) {x %>%
    select(from_id = reference_id, to_id = rodent_id, from_species = reference_name, to_species = initial_species_id, edge_sero = interpretation, )})

nodelist <- lapply(network, function(x) {
  
  from_rodents <- x %>%
    select(node_id = reference_id, species = reference_name, village, visit, grid_number, simple_habitat, interpretation)
  
  to_rodents <- x %>%
    filter(!rodent_id %in% from_rodents$node_id) %>%
    select(node_id = reference_id, species = initial_species_id, village, visit, grid_number, simple_habitat, interpretation)
  
  nodelist = bind_rows(from_rodents, to_rodents) %>%
    distinct()
  
  return(nodelist)
})

graph_landuse <- mapply(X = edgelist,
                        Y = nodelist,
                        function(X, Y) 
                        { g_landuse <- graph_from_data_frame(d = X, directed = FALSE)
                        
                        v_species <- Y %>% 
                          filter(node_id %in% V(g_landuse)$name)
                        
                        g_landuse <- set_vertex_attr(g_landuse, "Species", value = str_to_sentence(str_replace_all(v_species$species, "_", " ")))
                        g_landuse <- set_vertex_attr(g_landuse, "Village", value = str_to_sentence(v_species$village))
                        g_landuse <- set_vertex_attr(g_landuse, "Serostatus", value = str_to_sentence(v_species$interpretation))
                        g_landuse <- set_vertex_attr(g_landuse, "Landuse", value = str_to_sentence(v_species$simple_habitat))
                        g_landuse <- set_vertex_attr(g_landuse, "Visit", value = str_to_sentence(v_species$visit))
                        g_landuse <- set_vertex_attr(g_landuse, "Grid", value = str_to_sentence(v_species$grid_number))
                        
                        return(g_landuse)
                        })

plot_graphs <- list()

for(i in 1:length(graph_landuse)) {
  
  plot_graphs[[i]] <- ggraph(graph_landuse[[i]], layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    # geom_mark_hull(aes(x = x, y = y, fill = Village)) +
    theme_graph(base_family = 'Helvetica', title_size = 16) +
    scale_colour_manual(values = species_palette, drop = FALSE) +
    labs(title = str_to_sentence(str_replace_all(names(graph_landuse[i]), "_", " ")))
  
}

save_plot(plot = plot_graphs[[1]], filename = here("output", "agriculture_graph.pdf"), base_height = 12, base_width = 16)
save_plot(plot = plot_graphs[[2]], filename = here("output", "forest_graph.pdf"), base_height = 12, base_width = 16)
save_plot(plot = plot_graphs[[3]], filename = here("output", "village_graph.pdf"), base_height = 12, base_width = 16)

# Landuse network metrics

graph_metrics <- bind_graphs(as_tbl_graph(graph_landuse$agriculture),
                        as_tbl_graph(graph_landuse$secondary_forest),
                        as_tbl_graph(graph_landuse$village)) %>%
  activate(nodes) %>%
  mutate(n_neighbours = local_size(mindist = 0),
         degree = centrality_degree(),
         betweenness = centrality_betweenness(directed = FALSE),
         effective_network = node_effective_network_size(),
         closeness = centrality_closeness_harmonic()) %>%
  as_tibble()

graph_metrics %>%
  ggplot() +
  geom_bar(aes(x = degree)) +
  facet_wrap(~ Landuse) +
  theme_bw()

graph_metrics %>%
  ggplot() +
  geom_histogram(aes(x = betweenness)) +
  facet_wrap(~ Landuse) +
  theme_bw()

graph_metrics %>%
  ggplot() +
  geom_histogram(aes(x = effective_network)) +
  facet_wrap(~ Landuse) +
  theme_bw()

graph_metrics %>%
  ggplot() +
  geom_histogram(aes(x = closeness)) +
  facet_wrap(~ Landuse) +
  theme_bw()

# Networks for each visit and grid

# Agriculture
agricultural_graphs <- as_tbl_graph(graph_landuse$agriculture) %>%
  activate(nodes) %>%
  mutate(Species = factor(Species)) %>%
  to_split(Village, Grid, Visit, split_by = "nodes")

agricultural_plots <- lapply(agricultural_graphs, function(x) {
  if(nrow(x %>%
          activate(edges) %>%
          as_tibble()) > 
     nrow(x %>%
          activate(nodes) %>%
          as_tibble())) {
    plot_name <- paste0("Grid: ", unique(as_tibble(x)$Grid), ", Visit: ", unique(as_tibble(x)$Visit))
    
    plot <- ggraph(x, layout = "fr") +
      geom_edge_link() +
      geom_node_point(aes(colour = Species), size = 3) +
      scale_colour_manual(values = species_palette, drop = FALSE) +
      labs(title = plot_name) +
      theme_graph(title_size = 10,
                  foreground = "black")
    
    return(plot)
  } else {
    plot_name <- paste0("Grid: ", unique(as_tibble(x)$Grid), ", Visit: ", unique(as_tibble(x)$Visit))
    
    plot <- ggraph(x, layout = "fr") +
      geom_node_point(aes(colour = Species), size = 3) +
      scale_colour_manual(values = species_palette, drop = FALSE) +
      labs(title = plot_name) +
      theme_graph(title_size = 10,
                  foreground = "black")
    
    return(plot)
  }
})

baiama_agricultural <- (agricultural_plots[[1]] + agricultural_plots[[2]]) /
  (agricultural_plots[[3]] + agricultural_plots[[4]]) /
  (agricultural_plots[[5]] + agricultural_plots[[6]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Baiama - Agricultural")

write_rds(x = baiama_agricultural, file = here("report", "plots", "baiama_agriculture.rds"))

write_rds(x = baiama_agricultural, file = here("report", "plots", "baiama_agriculture.rds"))

lalehun_agricultural_1 <- (agricultural_plots[[10]] + agricultural_plots[[11]]  + agricultural_plots[[12]] + agricultural_plots[[13]] + agricultural_plots[[14]])  + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lalehun - Agricultural 1")
lalehun_agricultural_2 <- (agricultural_plots[[15]] + agricultural_plots[[16]] + agricultural_plots[[17]] + agricultural_plots[[18]] + agricultural_plots[[19]] + agricultural_plots[[20]])   + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lalehun - Agricultural 2")
lalehun_agricultural_3 <- (agricultural_plots[[21]] + agricultural_plots[[22]] + agricultural_plots[[23]] + agricultural_plots[[24]] + agricultural_plots[[25]] + agricultural_plots[[26]])   + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lalehun - Agricultural 3")
lalehun_agricultural_4 <- (agricultural_plots[[27]]) /
  (agricultural_plots[[28]] + agricultural_plots[[29]] + agricultural_plots[[30]] + agricultural_plots[[31]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lalehun - Agricultural 4")

lambayama_agricultural <- (agricultural_plots[[32]] + agricultural_plots[[33]] + agricultural_plots[[34]]) /
  agricultural_plots[[35]] /
  (agricultural_plots[[36]] + agricultural_plots[[37]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lambayama - Agricultural")

seilama_agricultural_1 <- (agricultural_plots[[38]] + agricultural_plots[[39]]  + agricultural_plots[[40]]) /
  (agricultural_plots[[41]] + agricultural_plots[[42]] + agricultural_plots[[43]])  + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Seilama - Agricultural 1")
seilama_agricultural_2 <- (agricultural_plots[[44]] + agricultural_plots[[45]] + agricultural_plots[[46]]) /
  (agricultural_plots[[47]] + agricultural_plots[[48]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Seilama - Agricultural 2")
seilama_agricultural_3 <- (agricultural_plots[[49]] + agricultural_plots[[50]] + agricultural_plots[[51]]) /
  (agricultural_plots[[52]] + agricultural_plots[[53]] + agricultural_plots[[54]])   + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Seilama - Agricultural 3")
seilama_agricultural_4 <- (agricultural_plots[[55]] + agricultural_plots[[56]] + agricultural_plots[[57]]) /
  (agricultural_plots[[58]] + agricultural_plots[[59]] + agricultural_plots[[60]])   + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Seilama - Agricultural 4")

# Forest
forest_graphs <- as_tbl_graph(graph_landuse$secondary_forest) %>%
  activate(nodes) %>%
  mutate(Species = factor(Species)) %>%
  to_split(Village, Grid, Visit, split_by = "nodes")

forest_plots <- lapply(forest_graphs, function(x) {
  if(nrow(x %>%
          activate(edges) %>%
          as_tibble()) > 
     nrow(x %>%
          activate(nodes) %>%
          as_tibble())) {
    plot_name <- paste0("Grid: ", unique(as_tibble(x)$Grid), ", Visit: ", unique(as_tibble(x)$Visit))
    
    plot <- ggraph(x, layout = "fr") +
      geom_edge_link() +
      geom_node_point(aes(colour = Species), size = 3) +
      scale_colour_manual(values = species_palette, drop = FALSE) +
      labs(title = plot_name) +
      theme_graph(title_size = 10,
                  foreground = "black")
    
    return(plot)
  } else {
    plot_name <- paste0("Grid: ", unique(as_tibble(x)$Grid), ", Visit: ", unique(as_tibble(x)$Visit))
    
    plot <- ggraph(x, layout = "fr") +
      geom_node_point(aes(colour = Species), size = 3) +
      scale_colour_manual(values = species_palette, drop = FALSE) +
      labs(title = plot_name) +
      theme_graph(title_size = 10,
                  foreground = "black")
    
    return(plot)
  }
})

baiama_forest <- forest_plots[[1]] +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Baiama - Forest")

lalehun_forest <- forest_plots[[3]] +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lalehun - Forest")

seilama_forest <- (forest_plots[[4]] + forest_plots[[5]] + forest_plots[[6]])/
  (forest_plots[[7]] + forest_plots[[8]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Seilama - Forest")

# Village

village_graphs <- as_tbl_graph(graph_landuse$village) %>%
  activate(nodes) %>%
  mutate(Species = factor(Species)) %>%
  to_split(Village, Grid, Visit, split_by = "nodes")

village_plots <- lapply(village_graphs, function(x) {
  if(nrow(x %>%
          activate(edges) %>%
          as_tibble()) > 
     nrow(x %>%
          activate(nodes) %>%
          as_tibble())) {
    plot_name <- paste0("Grid: ", unique(as_tibble(x)$Grid), ", Visit: ", unique(as_tibble(x)$Visit))
    
    plot <- ggraph(x, layout = "fr") +
      geom_edge_link() +
      geom_node_point(aes(colour = Species), size = 3) +
      scale_colour_manual(values = species_palette, drop = FALSE) +
      labs(title = plot_name) +
      theme_graph(title_size = 10,
                  foreground = "black")
    
    return(plot)
  } else {
    plot_name <- paste0("Grid: ", unique(as_tibble(x)$Grid), ", Visit: ", unique(as_tibble(x)$Visit))
    
    plot <- ggraph(x, layout = "fr") +
      geom_node_point(aes(colour = Species), size = 3) +
      scale_colour_manual(values = species_palette, drop = FALSE) +
      labs(title = plot_name) +
      theme_graph(title_size = 10,
                  foreground = "black")
    
    return(plot)
  }
})

baiama_village <- (village_plots[[1]] + village_plots[[2]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Baiama - Village")

lalehun_village <- (village_plots[[3]] + village_plots[[4]] + village_plots[[5]]) /
  (village_plots[[6]] + village_plots[[7]] + village_plots[[8]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lalehun - Village")

lambayama_village <- (village_plots[[9]] + village_plots[[10]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Lambayama - Village")

seilama_village <- (village_plots[[11]] + village_plots[[12]] + village_plots[[13]]) /
  (village_plots[[14]] + village_plots[[15]]) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Seilama - Village")

# Describing networks by habitat type
network_size <- bind_graphs(as_tbl_graph(graph_landuse$agriculture),
                            as_tbl_graph(graph_landuse$secondary_forest),
                            as_tbl_graph(graph_landuse$village)) %>%
  activate(nodes)%>%
  filter(Village != "Bambawo") %>%
  mutate(Species = factor(Species)) %>%
  to_split(Village, Grid, Visit, split_by = "nodes") %>%
  lapply(., function(x) {
  x %>% 
    activate(nodes) %>% 
    mutate(Individuals = graph_order(),
           Size = graph_size(),
           Contacts = unlist(graph_mutual_count()),
           Clique_size = graph_clique_num(),
           Number_components = graph_component_count(),
           Connectedness = graph_adhesion()) %>%
    as_tibble() %>%
    distinct(Landuse, Village, Grid, Visit, Individuals, Size, Contacts, Clique_size, Number_components, Connectedness)
}) %>%
  do.call(rbind.data.frame, .)

network_summary <- network_size %>%
  group_by(Landuse, Visit) %>%
  summarise(mean_individuals = mean(Individuals),
            sd_individuals = sd(Individuals),
            mean_contacts = mean(Contacts),
            sd_contacts = sd(Contacts))

ggplot(network_size) +
  geom_boxplot(aes(x = Visit, y = Individuals, fill = Landuse))

ggplot(network_size) +
  geom_boxplot(aes(x = Visit, y = Contacts, fill = Landuse))

# Calculating connectedness







as_tbl_graph(graph_landuse$agriculture) %>%
  activate(nodes) %>%
  filter(Village != "Bambawo") %>%
  group_by(Village, Visit, Grid) %>%
  mutate(n_individuals = n(),
         n_contacts = graph_size())

# Networks for positive rodents
positive_networks <- lapply(graph_landuse, function(x) {
  
  g2 <- as_tbl_graph(x) %>%
    activate(edges) %>%
    filter(edge_sero == "Positive") %>%
    activate(nodes) %>%
    filter(!node_is_isolated())})

plot_positive <- list()

for(i in 1:length(positive_networks)) {
  
  plot_positive[[i]] <- ggraph(positive_networks[[i]], layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Serostatus), size = 2) +
    # geom_mark_hull(aes(x = x, y = y, fill = Village)) +
    theme_graph(base_family = 'Helvetica', title_size = 16) +
    scale_colour_discrete() +
    labs(title = str_to_sentence(str_replace_all(names(positive_networks[i]), "_", " ")))
  
}
