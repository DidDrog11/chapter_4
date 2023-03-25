library(tidyterra)
library(terra)

rodents <- rodent_trapping_data$rodent_data %>%
  select(rodent_uid, species, field_id) %>%
  mutate(genus = str_to_sentence(str_split(coalesce(species, field_id), "_", simplify = TRUE)[, 1]))

a <- rodent_trapping_data$trap_data %>%
  select(village, visit, grid_number, trap_number, habitat_group, trap_uid, rodent_uid) %>%
  mutate(alpha = case_when(is.na(rodent_uid) ~ 0.1,
                           TRUE ~ 1),
         size = case_when(is.na(rodent_uid) ~ "trap",
                          TRUE ~ "rodent")) %>%
  filter(village != "bambawo") %>%
  left_join(rodents %>%
              select(rodent_uid, genus)) %>%
  distinct(rodent_uid, geometry, .keep_all = TRUE) %>%
  group_by(village) %>%
  group_split()

b <- lapply(a, vect)

pal <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
names(pal) <- levels(fct(rodents$genus))

c <- lapply(b, function(x) {
  
    ggplot() +
    geom_spatvector(data = x[is.na(x$genus)], aes(colour = fct(genus), alpha = alpha, size = size)) +
    geom_spatvector(data = x[!is.na(x$genus)], aes(colour = fct(genus), alpha = alpha, size = size)) +
    scale_colour_manual(values = pal,
                        drop = FALSE) +
    scale_size_manual(values = c(3, 1)) +
    theme_bw() +
    labs(colour = "Genus",
         title = str_to_sentence(unique(x$village))) +
    guides(alpha = "none",
           size = "none")
  
})

save_plot(plot = c[[1]], filename = here("output", "presentation_rodent_locations_1.png"), base_width = 10, base_height = 12)
save_plot(plot = c[[2]], filename = here("output", "presentation_rodent_locations_2.png"), base_width = 5, base_height = 7)
save_plot(plot = c[[3]], filename = here("output", "presentation_rodent_locations_3.png"), base_width = 8, base_height = 5)
save_plot(plot = c[[4]], filename = here("output", "presentation_rodent_locations_4.png"), base_width = 8, base_height = 5)
