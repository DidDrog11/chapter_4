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
