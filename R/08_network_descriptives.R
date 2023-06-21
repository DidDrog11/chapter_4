source(here::here("R", "00_setup.R"))

# Describing contact networks --------------------------------------------------------
# 25 networks have been produced one for each landuse type and visit
# These have been expanded with non-observed individuals based on estimated abundance
# Removing the network using rm(rodent_network) will speed up processing of code once it is no longer needed

rodent_network <- read_rds(here("data", "rodent_network_site_visit.rds"))

# Landuse network metrics
`%s%` <- network::`%s%`

rodent_net_descriptives <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  
  network_descriptives <- tibble(nodes = network.size(x),
                                 observed_nodes = table(x%v%"Observed")["TRUE"],
                                 unobserved_nodes = table(x%v%"Observed")["FALSE"],
                                 non_missing_edges = network.edgecount(x),
                                 missing_edges = network.naedgecount(x),
                                 median_degree_observed = median(igraph::degree(observed_igraph, mode = "out")),
                                 lower_quartile = quantile(igraph::degree(observed_igraph, mode = "out"), probs = 0.25),
                                 upper_quartile = quantile(igraph::degree(observed_igraph, mode = "out"), probs = 0.75),
                                 IQR = IQR(igraph::degree(observed_igraph, mode = "out")),
                                 mean_degree = mean(igraph::degree(observed_igraph, mode = "out")),
                                 sd_degree = sd(igraph::degree(observed_igraph, mode = "out")),
                                 mean_betweenness = mean(estimate_betweenness(observed_igraph, cutoff = -1)),
                                 sd_betweenness = sd(estimate_betweenness(observed_igraph, cutoff = -1)),
                                 landuse = unique(x%v%"Landuse"),
                                 n_species = length(unique(observed_subgraph%v%"Species")),
                                 density = edge_density(observed_igraph),
                                 density_observed = network.density(x, na.omit = FALSE),
                                 visit = unique(x%v%"Visit"))
  
  node_level_descriptives <- tibble(rodent_uid = V(observed_igraph)$vertex.names,
                                    species =  V(observed_igraph)$Species,
                                    village = V(observed_igraph)$Village,
                                    landuse = V(observed_igraph)$Landuse,
                                    visit = V(observed_igraph)$Visit,
                                    degree = igraph::degree(observed_igraph, mode = "out"),
                                    betweenness = igraph::estimate_betweenness(observed_igraph, cutoff = -1))
  
  return(list(network_level = network_descriptives,
              node_level = node_level_descriptives))
  
})

network_level_descriptives <- lapply(rodent_net_descriptives, function(x) x$network_level) %>%
  bind_rows() %>%
  mutate(network_number = row_number()) %>%
  relocate(network_number, landuse, visit, n_species) %>%
  mutate(network_number = fct_inorder(as.character(network_number)),
         landuse = fct(landuse, levels = c("Forest", "Agriculture", "Village")))

network_level_descriptives %>%
  group_by(landuse) %>%
  summarise(richness = max(n_species),
            nodes = sum(observed_nodes),
            edges = sum(non_missing_edges),
            max_degree = max(mean_degree),
            max_betweenness = max(mean_betweenness),
            max_density = max(density_observed))


# Network descriptive figures ---------------------------------------------

network_node_plot <- network_level_descriptives %>%
  select(network_number, landuse, visit, observed_nodes, unobserved_nodes) %>%
  pivot_longer(cols = c("observed_nodes", "unobserved_nodes"), names_to = "node", values_to = "n") %>%
  mutate(node = dt_case_when(str_detect(node, "unobserved") ~ "Not-observed",
                             str_detect(node, "observed") ~ "Observed"),
         network = fct_rev(fct_inorder(paste0(network_number, ": ", landuse, ", visit ", visit)))) %>%
  ggplot() +
  geom_col(aes(x = n, y = network, alpha = node, fill = landuse), position = position_stack(reverse = FALSE)) +
  scale_alpha_manual(values = c(0.5, 1),
                     guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = landuse_palette) +
  theme_bw() +
  labs(y = element_blank(),
       x = "Nodes (N)",
       fill = element_blank(),
       alpha = element_blank())


node_level_descriptives <- lapply(rodent_net_descriptives, function(x) x$node_level) %>%
  bind_rows() %>%
  left_join(network_numbers %>%
              select(-Grid) %>%
              rename("network_number" = Network) %>%
              distinct(),
            by = c("village" = "Village",
                   "visit" = "Visit",
                   "landuse" = "Landuse")) %>%
  relocate(network_number, rodent_uid, species, village, landuse, visit, degree) %>%
  mutate(network_number = fct_inorder(as.character(network_number)),
         landuse = fct(landuse, levels = c("Forest", "Agriculture", "Village")))


# Node level descriptives -------------------------------------------------


node_level_descriptives %>% 
  group_by(network_number, landuse, visit) %>%
  summarise(max_degree = max(degree)) %>%
  arrange(-max_degree)
  
node_level_descriptives %>% 
  group_by(landuse) %>%
  summarise(mean_degree = mean(degree)) %>%
  arrange(-mean_degree)

node_level_descriptives %>% 
  group_by(landuse) %>%
  summarise(mean_betweenness = mean(betweenness)) %>%
  arrange(-mean_betweenness)

node_level_descriptives %>% 
  group_by(species, landuse) %>%
  summarise(n = n(),
            max_degree = max(degree),
            mean_degree = mean(degree),
            sd_degree = sd(degree, na.rm = TRUE),
            median_degree = median(degree),
            max_betweenness = max(betweenness),
            mean_betweenness = mean(betweenness, na.rm = TRUE),
            sd_betweenness = sd(betweenness),
            median_betweenness = median(betweenness)) %>%
  arrange(-max_degree)


# Figure 2 ----------------------------------------------------------------

figure_2_updated <- node_level_descriptives %>%
  mutate(species = fct(species, levels = species_order_plots)) %>%
  ungroup() %>%
  ggplot(aes(x = species, y = degree, fill = landuse)) +
  geom_boxplot(position = position_dodge(0.9, preserve = "total"), alpha = 1, outlier.shape = 21) +
  coord_flip() +
  scale_fill_manual(values = landuse_palette) +
  scale_colour_manual(values = landuse_palette) +
  facet_wrap(~ species, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(fill = "Land use",
       colour = "Land use",
       x = "Species",
       y = "Degree")

save_plot(plot = figure_2_updated, filename = here("output", "Figure_2_updated.png"), base_height = 9, base_width = 6)

# Species interactions ----------------------------------------------------

species_contact_graphs <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  observed_tblgraph <- as_tbl_graph(observed_igraph)
  
})

species_contact_graph <- species_contact_graphs %>%
  bind_graphs()

observed_individuals <- species_contact_graph %>%
  activate(nodes) %>%
  data.frame() %>%
  rownames_to_column("rowid") %>%
  mutate(rowid = as.integer(rowid))

observed_contacts <- bind_rows(species_contact_graph %>%
                                 activate(edges) %>%
                                 data.frame() %>%
                                 tibble(),
                               species_contact_graph %>%
                                 activate(edges) %>%
                                 data.frame() %>%
                                 tibble() %>%
                                 rename("to" = "from",
                                        "from" = "to")) %>%
  distinct()

contacts <- observed_contacts %>%
  left_join(observed_individuals %>%
              select(-na, -Observed), by = c("from" = "rowid")) %>%
  rename(Species_from = Species,
         Vertex_from = vertex.names,
         Village_from = Village,
         Visit_from = Visit,
         Landuse_from = Landuse) %>%
  select(-from) %>%
  left_join(observed_individuals %>%
              select(-na, -Observed), by = c("to" = "rowid")) %>%
  rename(Species_to = Species,
         Vertex_to = vertex.names,
         Village_to = Village,
         Visit_to = Visit,
         Landuse_to = Landuse) %>%
  select(-to)

summarise_contacts <- contacts %>%
  group_by(Species_from, Species_to, Landuse_from, Landuse_to) %>%
  summarise(n = n()) %>%
  mutate(Species_from = factor(Species_from, levels = species_order_plots, ordered = TRUE),
         Species_to = factor(Species_to, levels = species_order_plots, ordered = TRUE),
         Landuse_from = fct(Landuse_from, levels = c("Forest", "Agriculture", "Village")))

contact_df <- summarise_contacts %>%
  rowwise() %>%
  rename(Landuse = Landuse_from) %>%
  group_by(Species_from, Landuse) %>%
  mutate(n_contacts = sum(n)) %>%
  group_by(Species_from, Species_to, n_contacts, Landuse) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(perc_contacts = round((n/n_contacts) * 100, 0),
         label_perc_contacts = case_when(perc_contacts == 0 ~ "<1%",
                                         TRUE ~ paste0(perc_contacts, "%")))

n_species_contacts <- contact_df %>%
  group_by(Species_from) %>%
  summarise(n_distinct_species = length(unique(Species_to))) %>%
  left_join(rodent_detail %>%
              group_by(species) %>%
              summarise(n = n()) %>%
              mutate(species = factor(species, levels = species_order_plots, ordered = TRUE)),
            by = c("Species_from" = "species"))

cor.test(n_species_contacts$n_distinct_species, n_species_contacts$n, method = "pearson")
  
forest_contact_plot <- contact_df %>%
  filter(Landuse == "Forest") %>%
  ggplot() +
  geom_tile(aes(x = Species_from, y = Species_to, fill = perc_contacts)) +
  geom_label(aes(x = Species_from, y = Species_to, label = label_perc_contacts)) +
  scale_x_discrete(drop = FALSE,  labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() +
  labs(x = "Contact to",
       y = "Contact from",
       title = "Forest",
       fill = "Percentage of contacts (%)")

ag_contact_plot <- contact_df %>%
  filter(Landuse == "Agriculture") %>%
  ggplot() +
  geom_tile(aes(x = Species_from, y = Species_to, fill = perc_contacts)) +
  geom_label(aes(x = Species_from, y = Species_to, label = label_perc_contacts)) +
  scale_x_discrete(drop = FALSE,  labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() +
  labs(x = "Contact to",
       y = "Contact from",
       title = "Agriculture",
       fill = "Percentage of contacts (%)")

vil_contact_plot <- contact_df %>%
  filter(Landuse == "Village") %>%
  ggplot() +
  geom_tile(aes(x = Species_from, y = Species_to, fill = perc_contacts)) +
  geom_label(aes(x = Species_from, y = Species_to, label = label_perc_contacts)) +
  scale_x_discrete(drop = FALSE,  labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() +
  labs(x = "Contact to",
       y = "Contact from",
       title = "Village",
       fill = "Percentage of contacts (%)")

save_plot(plot = forest_contact_plot +
            theme(legend.position = "bottom"), filename = here("output", "Supplementary_Figure_4a.png"), base_height = 12, base_width = 14)
save_plot(plot = ag_contact_plot +
            theme(legend.position = "bottom"), filename = here("output", "Figure_3_updated.png"), base_height = 10, base_width = 14)
save_plot(plot = vil_contact_plot +
            theme(legend.position = "bottom"), filename = here("output", "Supplementary_Figure_4b.png"), base_height = 12, base_width = 14)

# Node level descriptive figures ------------------------------------------

node_level_descriptives <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  
  descriptives <- tibble(species = observed_subgraph%v%"Species",
                         landuse = observed_subgraph%v%"Landuse",
                         visit = observed_subgraph%v%"Visit",
                         village = observed_subgraph%v%"Village",
                         degree = igraph::degree(observed_igraph, mode = "out"))
}) %>%
  bind_rows(.id = "network_number") %>%
  relocate(network_number, landuse, visit, village, species) %>%
  mutate(network_number = fct(network_number),
         species = fct(species, levels = species_order_plots),
         landuse = fct(landuse, levels = c("Forest", "Agriculture", "Village")))


plot_landuse_degree <- node_level_descriptives %>% 
  mutate(species = fct_rev(species),
         landuse = fct_rev(landuse)) %>%
  ggplot() + 
  geom_dotplot(aes(y = degree, x = landuse, fill = landuse, colour = landuse),
               binaxis = "y",
               stackdir = "center",
               position = position_jitter(width = 0.02, height = 0.08),
               binwidth = 0.2) +
  scale_fill_manual(values = landuse_palette,
                    guide = "none") +
  scale_colour_manual(values = landuse_palette,
                    guide = "none") +
  labs(x = element_blank(),
       y = "Degree",
       fill = "Landuse") +
  theme_bw() +
  coord_flip()

plot_species_degree <- node_level_descriptives %>% 
  mutate(species = fct_rev(species),
         landuse = fct_rev(landuse)) %>%
  ggplot() + 
  geom_boxplot(aes(y = species, x = degree, fill = landuse), position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = landuse_palette,
                    guide = guide_legend(reverse = TRUE)) +
  labs(y = element_blank(),
       x = "Degree",
       fill = "Landuse") +
  theme_bw()

save_plot(plot_grid(plot_landuse_degree, plot_species_degree, labels = c("A", "B"), rel_widths = c(0.6, 1)),
          filename = here("output", "Figure_2.png"), base_width = 10, base_height = 9)
  
# Plot networks -----------------------------------------------------------

network_plots <- function(network, fig_label) {

  a <- as.edgelist(network)
  vnames <- tibble(node = attr(a, "vnames"))
  obs_edges <- tibble(from = attr(a, "vnames")[a[, 1]],
                  to =  attr(a, "vnames")[a[, 2]])
  gg_graph <- tbl_graph(nodes = vnames, edges = obs_edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(Landuse = network%v%"Landuse",
           Observed = network%v%"Observed",
           Species = network%v%"Species",
           Village = network%v%"Village",
           Visit = network%v%"Visit")
    
  fig <- ggraph(gg_graph, layout = "kk")  +
    geom_edge_fan() +
    geom_node_point(aes(colour = Species, alpha = Observed)) +
    scale_alpha_manual(values = c(0.1, 1)) +
    scale_colour_manual(values = species_palette) +
    theme_graph() +
    guides(alpha = "none") +
    labs(title = paste0(unique(gg_graph %>%
                                 pull(Landuse)), ", visit ",
                        unique(gg_graph %>%
                                 pull(Visit))))
  
  save_plot(filename = here("output", paste0("Supplementary_Figure_3", fig_label, ".png")), plot = fig, base_width = 8, base_height = 6)

  }

for(i in 1:length(rodent_network)) {
  
  network_plots(network = rodent_network[[i]], fig_label = LETTERS[i])
  
}


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
         Landuse = fct_rev(factor(Landuse, levels = c("Forest", "Agriculture", "Village"))))

species_degree %>%
  ggplot() +
  geom_point(aes(x = Degree, y = Landuse, colour = Landuse), position = position_jitter(width = 0)) +
  geom_violin(aes(x = Degree, y = Landuse, fill = Landuse), alpha = 0.5) +
  labs(title = "Degree by species and landuse",
       x = "Degree") +
  scale_colour_manual(values = landuse_palette) +
  scale_fill_manual(values = landuse_palette) +
  facet_wrap(~ Species, ncol = 1) +
  theme_bw()

# Association of degree with antibody status ------------------------------

elisa_status <- combined_data$ELISA_enriched %>%
  group_by(rodent_uid) %>%
  arrange(interpretation) %>%
  filter(rodent_uid %in% rodents$rodent_uid) %>%
  select(rodent_uid, interpretation, clean_names) %>%
  left_join(node_level_descriptives,
            by = "rodent_uid") %>%
  drop_na(interpretation)

plot(elisa_status$interpretation, elisa_status$degree)
plot(elisa_status$interpretation, elisa_status$betweenness)

networks_containing_positive <- elisa_status %>%
  filter(interpretation == "Positive") %>%
  distinct(network_number)

positive_m_nat_landuse <- elisa_status %>%
  filter(interpretation == "Positive") %>%
  filter(species == "Mastomys natalensis")
