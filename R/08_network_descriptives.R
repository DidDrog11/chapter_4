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
                         mean_betweenness = mean(estimate_betweenness(observed_igraph, cutoff = -1)),
                         sd_betweenness = sd(estimate_betweenness(observed_igraph, cutoff = -1)),
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

rodent_net_descriptives %>%
  group_by(landuse) %>%
  summarise(median_edges = median(non_missing_edges),
            IQR_edges = IQR(non_missing_edges))

landuse_level_degree <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  
  return(observed_igraph)
  
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
       alpha = element_blank())

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


# Species interactions ----------------------------------------------------

species_contact_graph <- lapply(rodent_network, function(x) {
  
  observed <- x%v%"Observed"
  observed_subgraph <- x %s% which(observed == TRUE)
  observed_igraph <- asIgraph(observed_subgraph)
  observed_tblgraph <- as_tbl_graph(observed_igraph)
  
}) %>%
  bind_graphs()

observed_individuals <- species_contact_graph %>%
  activate(nodes) %>%
  data.frame() %>%
  rownames_to_column("rowid") %>%
  mutate(rowid = as.integer(rowid))

observed_contacts <- species_contact_graph %>%
  activate(edges) %>%
  data.frame()

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

triangle_contact <- summarise_contacts %>%
  rowwise() %>%
  mutate(species_1 = min(Species_from, Species_to),
         species_2 = max(Species_from, Species_to)) %>%
  rename(Landuse = Landuse_from) %>%
  group_by(species_1, species_2, Landuse) %>%
  summarise(n = sum(n)) %>%
  ungroup()

fill_triangle <- triangle_contact %>%
  expand(species_1, species_2, Landuse) %>%
  mutate(n_0 = dt_case_when(species_1 <= species_2 ~ as.integer(0),
                            TRUE ~ as.integer(NA)))

contact_df <- fill_triangle %>%
  left_join(triangle_contact,
            by = c("species_1", "species_2", "Landuse")) %>%
  mutate(n = coalesce(n, n_0))

forest_contact <- ggplot(contact_df %>%
                                filter(Landuse == "Forest") %>%
                                mutate(alpha = dt_case_when(n == 0 ~ 0,
                                                            TRUE ~ 1),
                                       n_discrete = cut(n, breaks = c(0, 1, 5, 10, 20, 50, 100, 200), include.lowest = TRUE, right = FALSE))) +
  geom_tile(aes(x = species_1, y = species_2, fill = n_discrete, alpha = alpha)) +
  scale_x_discrete(drop = FALSE, labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis(direction = -1, discrete = TRUE, na.translate = FALSE) +
  scale_alpha(range = c(0.1, 1)) +
  guides(alpha = "none") +
  facet_wrap(~ Landuse) +
  theme_bw() +
  labs(x = element_blank(),
       y = element_blank(),
       alpha = element_blank(),
       fill = "Contacts (N)")

agriculture_contact <- ggplot(contact_df %>%
                                filter(Landuse == "Agriculture") %>%
                                mutate(alpha = dt_case_when(n == 0 ~ 0,
                                                            TRUE ~ 1),
                                       n_discrete = cut(n, breaks = c(0, 1, 5, 10, 20, 50, 100, 200), include.lowest = TRUE, right = FALSE))) +
  geom_tile(aes(x = species_1, y = species_2, fill = n_discrete, alpha = alpha)) +
  scale_x_discrete(drop = FALSE, labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis(direction = -1, discrete = TRUE, na.translate = FALSE) +
  scale_alpha(range = c(0.1, 1)) +
  guides(alpha = "none") +
  facet_wrap(~ Landuse) +
  theme_bw() +
  labs(x = element_blank(),
       y = element_blank(),
       alpha = element_blank(),
       fill = "Contacts (N)")

village_contact <- ggplot(contact_df %>%
                            filter(Landuse == "Village") %>%
                            mutate(alpha = dt_case_when(n == 0 ~ 0,
                                                        TRUE ~ 1),
                                   n_discrete = cut(n, breaks = c(0, 1, 5, 10, 20, 50, 100, 200), include.lowest = TRUE, right = FALSE))) +
  geom_tile(aes(x = species_1, y = species_2, fill = n_discrete, alpha = alpha)) +
  scale_x_discrete(drop = FALSE, labels = function(x) str_wrap(x, width = 10)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_viridis(direction = -1, discrete = TRUE, na.translate = FALSE) +
  scale_alpha(range = c(0.1, 1)) +
  guides(alpha = "none") +
  facet_wrap(~ Landuse) +
  theme_bw() +
  labs(x = element_blank(),
       y = element_blank(),
       alpha = element_blank(),
       fill = "Contacts (N)")

combined_contact <- plot_grid(plotlist = list(forest_contact, agriculture_contact, village_contact), ncol = 1, labels = c("A", "B", "C"))

save_plot(plot = combined_contact, filename = here("output", "Figure_3.png"), base_height = 10, base_width = 12)

n_contacts_n_observed <- bind_rows(contacts,
                                   contacts %>%
                                     rename(Species_to = Species_from,
                                            Species_from = Species_to)) %>%
  distinct(Species_from, Species_to) %>%
  filter(Species_from != Species_to) %>%
  group_by(Species_from) %>%
  summarise(n_species = n()) %>%
  left_join(observed_individuals %>%
              group_by(Species) %>%
              summarise(n_individuals = n()),
            by = c("Species_from" = "Species")) %>%
  arrange(-n_species)

correlation_n_species <- cor.test(n_contacts_n_observed$n_species, n_contacts_n_observed$n_individuals, method = "pearson")


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
