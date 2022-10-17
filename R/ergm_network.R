graph_landuse <- read_rds(here("data", "graphs_landuse.rds"))

rodent_graph <- bind_graphs(graph_landuse$agriculture, graph_landuse$secondary_forest, graph_landuse$village)

# Species with more than 10 nodes
retain_rodents <- names(table(V(rodent_graph)$Species))[table(V(rodent_graph)$Species) > 10]

rodent_graph <- rodent_graph %>%
  mutate(Season = case_when(Visit %in% c("1", "2", "5", "6") ~ "Dry",
                            TRUE ~ "Wet"),
         Species = case_when(!Species %in% retain_rodents ~ "Other",
                             TRUE ~ Species),
         Village_loc = case_when(Village == "Lambayama" ~ "Urban",
                                 TRUE ~ "Rural"))

rodent_connected <- rodent_graph %>%
  activate(edges) %>%
  filter(from != to) %>%
  activate(nodes) %>%
  filter(!node_is_isolated())

# Convert to network and plot

rodent_graph <- igraph::simplify(as.igraph(rodent_graph))

rodent_net <- asNetwork(rodent_graph)

summary(rodent_net)

# Descriptive plots

ggraph(rodent_net, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Species)) +
  scale_colour_brewer(palette = "Dark2") +
  theme_graph() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
ggraph(rodent_net, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Village_loc)) +
  scale_colour_brewer(palette = "Set1") +
  labs(colour = "Village location") +
  theme_graph() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
ggraph(rodent_net, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Landuse)) +
  scale_colour_brewer(palette = "Set2") +
  theme_graph() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
ggraph(rodent_net, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Serostatus)) +
  scale_colour_brewer(palette = "Set2") +
  theme_graph() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# Network descriptives
rodent_net

rodent_net_descriptives <- list("mean_degree" = mean(degree(rodent_net, gmode = "graph")),
                                "sd_degree" = sd(degree(rodent_net, gmode = "graph")),
                                "degree_table" = table(degree(rodent_net, gmode = "graph")),
                                "triad_table" = triad.census(rodent_net, mode = "graph"),
                                "density" = graph.density(asIgraph(rodent_net)))

# Random graph of similar size and density
randomg <- rgraph(412, m = 1, tprob = rodent_net_descriptives$density, mode = "graph")
random_net <- as.network(randomg, directed = FALSE)

# Comparing distribution of degree between rodent network and random graph
plot_grid(plotlist = list(random = tibble("degree" = degree(random_net, gmode = "graph")) %>%
                            ggplot() + 
                            geom_bar(aes(x = degree)) +
                            labs(title = "Distribution of degree (Random)") +
                            theme_bw(),
                          rodent = tibble("degree" = degree(rodent_net, gmode = "graph")) %>%
                            ggplot() +
                            geom_bar(aes(x = degree)) +
                            labs(title = "Distribution of degree (Rodent)") +
                            theme_bw()))

summary(rodent_net ~ edges + nodematch("Species", diff = T) + nodematch("Landuse", diff = T))

rodent.01 <- ergm(rodent_net ~ edges + nodematch("Species", diff = T, levels = retain_rodents) + nodematch("Landuse", diff = T)) # Estimate the model 

summary(rodent.01)

# Removing those without evidence of contact

rodent_connected <- igraph::simplify(as.igraph(rodent_connected))

rodent_net_con <- asNetwork(rodent_connected)

rodent_net_con <- read_rds(here("rodent_net.rds"))

ggraph(rodent_net_con, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Species)) +
  facet_nodes(~ Landuse,  switch = "x", scales = "free")

rodent.02 <- ergm(rodent_net_con ~ edges + nodefactor("Species") + nodefactor("Landuse")) # Estimate the model 

summary(rodent.02)

rodent.03 <- ergm(rodent_net_con ~ edges + nodefactor("Species") + nodematch("Species", diff = T)) # Estimate the model 

summary(rodent.03)

rodent.04 <- ergm(rodent_net_con ~ edges + nodefactor("Species") + nodematch("Species", diff = T) + nodefactor("Landuse") + nodematch("Landuse", diff = T)) # Estimate the model 

summary(rodent.04)

rodent.04_gof <- gof(rodent.04 ~ degree + espartners + distance + model)

r4_prediction <- predict(rodent.04, conditional = FALSE, output = "data.frame", nsim = 1000, type = "response")

r4_data <- tibble(node = 1:length((rodent_net_con %v% "Species")),
                  species = (rodent_net_con %v% "Species"),
                  landuse = (rodent_net_con %v% "Landuse"))

test <- r4_prediction %>%
  left_join(r4_data %>%
              rename("species_tail" = "species",
                     "landuse_tail" = "landuse"), by = c("tail" = "node")) %>%
  left_join(r4_data %>%
              rename("species_head" = "species",
                     "landuse_head" = "landuse"), by = c("head" = "node")) %>%
  group_by(species_tail, species_head, landuse_head, landuse_tail)

# Probability of minutoides-minutoides link
plogis(-23.5 + 0.44 + 0.58 + 1.91 + 15.8)
# Probability of minutoides-mus link
plogis(-23.5 + 0.44 + 1.91 + 15.8)

# Probability of a mastomys forming a tie, with another mastomys, in village settings
plogis(-23.47 + -0.53 + 2.46 + 1.15 + 17.61)
# Probability of a mastomys forming a tie, with a mus, in agricultural settings
plogis(-23.47 + -0.53 + 1.15 + 17.61)

# Probability of a mus forming a tie, with another mus, in village settings
plogis(-23.47 + -1.91 + 5.15 + 1.15 + 17.6)
# Probability of a mus forming a tie, with a non-mus, in village settings
plogis(-23.47 + -1.91 + 1.15 + 17.6)
