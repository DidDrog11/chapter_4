# Not currently used as ELISA has only performed on a subset
# unique_rodents <- combined_data$ELISA_enriched %>%
#   group_by(rodent_uid) %>%
#   arrange(interpretation) %>%
#   slice(1) %>%
#   ungroup() %>%
#   mutate(village = case_when(is.na(village) & str_detect(rodent_uid, "LAL") ~ "lalehun",
#                              is.na(village) & str_detect(rodent_uid, "SEI") ~ "seilama",
#                              is.na(village) & str_detect(rodent_uid, "LAM") ~ "lambayama",
#                              is.na(village) & str_detect(rodent_uid, "BAI") ~ "baiama",
#                              TRUE ~ village),
#          visit = case_when(is.na(visit) ~ as.numeric(str_split(rodent_uid, "_", simplify = TRUE)[, 1]),
#                            TRUE ~ as.numeric(visit)),
#          simple_habitat = case_when(str_detect(site_habitat, "agriculture") ~ "agriculture",
#                                     str_detect(site_habitat, "village") ~ "village",
#                                     str_detect(site_habitat, "forest") ~ "secondary_forest",
#                                     str_detect(site_habitat, "fallow") ~ "agriculture",
#                                     is.na(site_habitat) & str_detect(habitat_group, "agriculture") ~ "agriculture",
#                                     is.na(site_habitat) & str_detect(rodent_uid, "3_LAM_001") ~ "agriculture",
#                                     is.na(site_habitat) & str_detect(rodent_uid, "LAL") ~ "agriculture",
#                                     is.na(site_habitat) & str_detect(rodent_uid, "SEI") ~ "agriculture"))
# 
# write_rds(unique_rodents, here("data", "unique_rodents.rds"))

# trap_locations <- trap_data$sites %>%
#   mutate(simple_habitat = case_when(str_detect(site_habitat, "agriculture") ~ "agriculture",
#                                     str_detect(site_habitat, "village") ~ "village",
#                                     str_detect(site_habitat, "forest") ~ "secondary_forest",
#                                     str_detect(site_habitat, "fallow") ~ "agriculture",
#                                     str_detect(habitat_group, "agriculture") ~ "agriculture"))


# Using rodent trapping extract -------------------------------------------------
add_grid_number <- rodent_trapping_data$rodent_data %>%
  left_join(tibble(rodent_trapping_data$trap_data) %>%
              select(trap_uid, grid_number), by = "trap_uid")


unique_rodents <- add_grid_number %>%
  filter(!str_detect(rodent_uid, "BAM")) %>%
  mutate(species = dt_case_when( 
    # using PCR allocated species
    !is.na(species) & str_detect(species, "crocidura") ~ "crocidura_spp",
    !is.na(species) & str_detect(species, "gerbilliscus") ~ "gerbilliscus_guineae",
    !is.na(species) & str_detect(species, "hybomys") ~ "hybomys_planifrons",
    !is.na(species) & str_detect(species, "hylomyscus") ~ "hylomyscus_simus",
    !is.na(species) & str_detect(species, "lemniscomys") ~ "lemniscomys_striatus",
    !is.na(species) & str_detect(species, "lophuromys") ~ "lophuromys_sikapusi",
    !is.na(species) & str_detect(species, "malacomys") ~ "malacomys_edwardsi",
    !is.na(species) & str_detect(species, "mastomys") ~ "mastomys_natalensis",
    !is.na(species) & str_detect(species, "musculus") ~ "mus_musculus",
    !is.na(species) & str_detect(species, "setulosus") ~ "mus_setulosus",
    !is.na(species) & str_detect(species, "praomys") ~ "praomys_rostratus",
    !is.na(species) & str_detect(species, "rattus") ~ "rattus_rattus",
    # using field id and location
    str_detect(rodent_uid, "LAM") & str_detect(field_id, "mus_") ~ "mus_musculus",
    grid_number != "7" & str_detect(field_id, "mus_") ~ "mus_setulosus",
    str_detect(field_id, "mus_") ~ "mus_musculus",
    # using field id only
    str_detect(field_id, "crocidura") ~ "crocidura_spp",
    str_detect(field_id, "gerbilliscus") ~ "gerbilliscus_guineae",
    str_detect(field_id, "hybomys") ~ "hybomys_planifrons",
    str_detect(field_id, "hylomyscus") ~ "hylomyscus_simus",
    str_detect(field_id, "lemniscomys") ~ "lemniscomys_striatus",
    str_detect(field_id, "lophuromys") ~ "lophuromys_sikapusi",
    str_detect(field_id, "malacomys") ~ "malacomys_edwardsi",
    str_detect(field_id, "mastomys") ~ "mastomys_natalensis",
    str_detect(field_id, "praomys") ~ "praomys_rostratus",
    str_detect(field_id, "rattus") ~ "rattus_rattus",
    
    TRUE ~ field_id))

trap_locations <- tibble(rodent_trapping_data$trap_data %>%
                           filter(!str_detect(village, "bambawo")) %>%
                           mutate(landuse = fct(dt_case_when(str_detect(village, "seilama") & grid_number %in% c(1, 2, 3, 4) ~ "agriculture",
                                                             str_detect(village, "seilama") & grid_number %in% c(5) ~ "forest",
                                                             str_detect(village, "seilama") & grid_number %in% c(6, 7) ~ "village",
                                                             str_detect(village, "lalehun") & grid_number %in% c(1, 2, 3, 5) ~ "agriculture",
                                                             str_detect(village, "lalehun") & grid_number %in% c(4) ~ "forest",
                                                             str_detect(village, "lalehun") & grid_number %in% c(6, 7) ~ "village",
                                                             str_detect(village, "lambayama") & grid_number %in% c(1, 2, 3, 4) ~ "agriculture",
                                                             str_detect(village, "lambayama") & grid_number %in% c(6, 7) ~ "village",
                                                             str_detect(village, "baiama") & grid_number %in% c(2, 3, 4) ~ "agriculture",
                                                             str_detect(village, "baiama") & grid_number %in% c(1) ~ "forest",
                                                             str_detect(village, "baiama") & grid_number %in% c(6, 7) ~ "village"),
                                                levels = c("forest", "agriculture", "village")),
                                  lon = st_coordinates(geometry)[, 1],
                                  lat = st_coordinates(geometry)[, 2])) %>%
  left_join(visit_season %>%
              mutate(visit = factor(visit)), by = "visit") %>%
  rename(grid = grid_number) %>%
  select(-geometry)

write_rds(unique_rodents, here("data", "chapter_4_rodents.rds"))
write_rds(trap_locations, here("data", "chapter_4_traps.rds"))
