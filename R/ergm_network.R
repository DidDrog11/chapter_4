graph_landuse <- read_rds(here("data", "graphs_landuse.rds"))

rodent_graph <- bind_graphs(graph_landuse$agriculture, graph_landuse$secondary_forest, graph_landuse$village)

rodent_graph <- rodent_graph %>%
  mutate(Season = case_when(Visit %in% c("1", "2", "5", "6") ~ "Dry",
                            TRUE ~ "Wet"))

rodent_connected <- rodent_graph %>%
  activate(edges) %>%
  filter(from != to) %>%
  activate(nodes) %>%
  filter(!node_is_isolated())

# Convert to network and plot

rodent_graph <- igraph::simplify(as.igraph(rodent_graph))

rodent_net <- asNetwork(rodent_graph)

ggraph(rodent_net, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = Species)) +
  facet_nodes(~ Landuse + Season,  switch = "x", scales = "free")

# All nodes and species with more than 10 nodes
retain_rodents <- names(table(V(rodent_graph)$Species))[table(V(rodent_graph)$Species) > 10]

summary(rodent_net ~ edges + nodematch("Species", diff = T, levels = retain_rodents) + nodematch("Landuse", diff = T))

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

summary(rodent_net_con ~ edges + nodefactor("Species") + nodematch("Species", diff = T) + nodefactor("Landuse") + nodematch("Landuse", diff = T))

rodent.02 <- ergm(rodent_net_con ~ edges + nodefactor("Species") + nodematch("Species", diff = T) + nodefactor("Landuse") + nodematch("Landuse", diff = T)) # Estimate the model 

summary(rodent.02)

# Probability of a mastomys forming a tie, with another mastomys, in agricultural settings
plogis(-4.59 + -0.57 + 2.42)
