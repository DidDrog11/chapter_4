unique_rodents <- combined_data$ELISA_enriched %>%
  group_by(rodent_uid) %>%
  arrange(interpretation) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(village = case_when(is.na(village) & str_detect(rodent_uid, "LAL") ~ "lalehun",
                             is.na(village) & str_detect(rodent_uid, "SEI") ~ "seilama",
                             is.na(village) & str_detect(rodent_uid, "LAM") ~ "lambayama",
                             TRUE ~ village),
         month_set = month(date_set),
         year_set = year(date_set),
         visit = factor(case_when(year_set == 2020 & month_set %in% c(11, 12) ~ 1,
                           year_set == 2021 & month_set %in% c(4) ~ 2,
                           year_set == 2021 & month_set %in% c(6, 7) ~ 3,
                           year_set == 2021 & month_set %in% c(10) ~ 4,
                           year_set == 2022 & month_set %in% c(1, 2) ~ 5,
                           year_set == 2022 & month_set %in% c(4) ~ 6,
                           is.na(year_set) & str_detect(rodent_uid, "3_LAM") ~ 5,
                           is.na(year_set) & str_detect(rodent_uid, "4_SEI|4_LAL") ~ 4)),
         simple_habitat = case_when(str_detect(site_habitat, "agriculture") ~ "agriculture",
                                    str_detect(site_habitat, "village") ~ "village",
                                    str_detect(site_habitat, "forest") ~ "secondary_forest",
                                    str_detect(site_habitat, "fallow") ~ "agriculture",
                                    is.na(site_habitat) & str_detect(habitat_group, "agriculture") ~ "agriculture",
                                    is.na(site_habitat) & str_detect(rodent_uid, "3_LAM_001") ~ "agriculture",
                                    is.na(site_habitat) & str_detect(rodent_uid, "LAL") ~ "agriculture",
                                    is.na(site_habitat) & str_detect(rodent_uid, "SEI") ~ "agriculture"))

write_rds(unique_rodents, here("data", "unique_rodents.rds"))

trap_locations <- trap_data$sites %>%
  mutate(simple_habitat = case_when(str_detect(site_habitat, "agriculture") ~ "agriculture",
                                    str_detect(site_habitat, "village") ~ "village",
                                    str_detect(site_habitat, "forest") ~ "secondary_forest",
                                    str_detect(site_habitat, "fallow") ~ "agriculture",
                                    str_detect(habitat_group, "agriculture") ~ "agriculture"))


# Using chapter 3 extract -------------------------------------------------

unique_rodents <- chapter_3_extract$detections %>%
  mutate(species = case_when(village == "lambayama" & str_detect(species, "mus_") ~ "mus_musculus",
                                        grid_number != "7" & str_detect(species, "mus_") ~ "mus_minutoides",
                                        str_detect(species, "mus_") ~ "mus_musculus",
                                        str_detect(species, "mastomys_") ~ "mastomys_natalensis",
                                        str_detect(species, "rattus_") ~ "rattus_rattus",
                                        str_detect(species, "lophuromys_") ~ "lophuromys_sikapusi",
                                        TRUE ~ species))

trap_locations <- chapter_3_extract$sites

write_rds(unique_rodents, here("data", "chapter_4_rodents.rds"))
write_rds(trap_locations, here("data", "chapter_4_traps.rds"))
