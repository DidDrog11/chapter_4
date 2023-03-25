source(here::here("R", "00_setup.R"))

# Describing contact networks --------------------------------------------------------
# 25 networks have been produced one for each landuse type and visit
# These have been expanded with non-observed individuals based on estimated abundance

rodent_network <- read_rds(here("data", "rodent_network_site_visit.rds"))

# Landuse network metrics
`%s%` <- network::`%s%`

rodent_net_descriptives <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  
  descriptives <- tibble(nodes = network.size(x),
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
                         landuse = unique(x%v%"Landuse"),
                         n_species = length(unique(x%v%"Species")),
                         density = edge_density(observed_igraph),
                         density_observed = network.density(x, na.omit = FALSE),
                         visit = unique(x%v%"Visit"))
  
}) %>%
  bind_rows() %>%
  mutate(network_number = row_number()) %>%
  relocate(network_number, landuse, visit, n_species) %>%
  mutate(network_number = fct_inorder(as.character(network_number)),
         landuse = fct(landuse, levels = c("Forest", "Agriculture", "Village")))

landuse_level_degree <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  
}) %>%
  bind_graphs() %>%
  mutate(degree = centrality_degree())

landuse_level_degree %>%
  group_by(Landuse) %>%
  mutate(mean = mean(degree),
         sd = sd(degree)) %>%
  distinct(Landuse, mean, sd)

# Network descriptive figures ---------------------------------------------

network_node_plot <- rodent_net_descriptives %>%
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
       alpha = "Opacity")

observed_edges_plot <- rodent_net_descriptives %>%
  select(network_number, landuse, visit, non_missing_edges) %>%
  mutate(network = fct_rev(fct_inorder(paste0(network_number, ": ", landuse, ", visit ", visit))),
         network_number = fct_rev(network_number)) %>%
  ggplot() +
  geom_col(aes(y = network_number, x = non_missing_edges, fill = landuse)) +
  scale_fill_manual(values = landuse_palette) +
  theme_bw() +
  theme() +
  guides(fill = "none") +
  labs(y = element_blank(),
       x = "Observed edges (N)",
       fill = element_blank())

save_plot(plot = plot_grid(network_node_plot, observed_edges_plot, rel_widths = c(1, 0.5), labels = c("A", "B")),
          filename = here("output", "Figure_1.png"), base_width = 10, base_height = 8)
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
         species = dt_case_when(str_detect(species, "Lemnis") ~ "Lemniscomys striatus",
                                str_detect(species, "Mala") ~ "Malacomys edwardsi",
                                TRUE ~ species),
         species = fct(species, levels = species_order_all),
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
          filename = here("output", "Figure_2.png"), base_width = 10, base_height = 7)
  
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
