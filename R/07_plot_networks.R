
plot_graph <- list()

plot_graph <- lapply(assemblages, function(x) {
  as_tbl_graph(x$graph) %>%
    mutate(Species = factor(Species, levels = names(species_palette))) %>%
    ggraph(layout = "fr") +
    geom_edge_link() +
    geom_node_point(aes(colour = Species), size = 3) +
    scale_colour_manual(values = species_palette, drop = FALSE) +
    labs(colour = "Species",
         title = paste0("Visit ", unique(V(x$graph)$Visit), ": ", unique(V(x$graph)$Landuse)),
         caption = "30m",
         x = element_blank(),
         y = element_blank()) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
})

assemblages_forest <- plot_grid(plotlist = list(plot_graph[[1]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[2]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[3]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[4]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[5]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[6]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[7]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[8]] +
                                                  theme(legend.position = "none"),
                                                plot_graph[[9]] +
                                                  theme(legend.position = "none")),
                                ncol = 3)

assemblages_agriculture <- plot_grid(plotlist = list(plot_graph[[10]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[11]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[12]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[13]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[14]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[15]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[16]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[17]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[18]] +
                                                       theme(legend.position = "none"),
                                                     plot_graph[[19]] +
                                                       theme(legend.position = "none")),
                                     ncol = 4)

assemblages_village <- plot_grid(plotlist = list(plot_graph[[20]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[21]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[22]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[23]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[24]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[25]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[26]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[27]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[28]] +
                                                   theme(legend.position = "none"),
                                                 plot_graph[[29]] +
                                                   theme(legend.position = "none")),
                                 ncol = 4)

save_plot(plot = assemblages_forest, here("output", "forest_landuse_graph.png"), base_width = 10, base_height = 12)
save_plot(plot = assemblages_agriculture, here("output", "agriculture_landuse_graph.png"), base_width = 16, base_height = 12)
save_plot(plot = assemblages_village, here("output", "village_landuse_graph.png"), base_width = 16, base_height = 12)
