combined_data <- readRDS(gzcon(url("https://github.com/DidDrog11/SL_lassa_ELISA/raw/main/output/ELISA_output.rds")))

# Classify rodent species for now based on locations
combined_data$ELISA_enriched <- combined_data$ELISA_enriched %>%
  mutate(initial_species_id = case_when(village == "lambayama" & str_detect(initial_species_id, "mus_") ~ "mus_musculus",
                                        grid_number != "7" & str_detect(initial_species_id, "mus_") ~ "mus_minutoides",
                                        str_detect(initial_species_id, "mus_") ~ "mus_musculus",
                                        str_detect(initial_species_id, "mastomys_") ~ "mastomys_natalensis",
                                        str_detect(initial_species_id, "rattus_") ~ "rattus_rattus",
                                        str_detect(initial_species_id, "lophuromys_") ~ "lophuromys_sikapusi",
                                        TRUE ~ initial_species_id))

trap_data <- readRDS(gzcon(url("https://github.com/DidDrog11/chapter_3/raw/observed_data/data/data_for_export/chapter_4_extract.rds")))

write_rds(combined_data, here("input", "ELISA_data.rds"))
write_rds(trap_data, here("input", "trap_data.rds"))

# Use observed data, no ELISA ---------------------------------------------

rodent_trapping_data <- readRDS(gzcon(url("https://github.com/DidDrog11/rodent_trapping/raw/main/data/data_for_export/combined_data.rds")))

