
unique_rodents <- read_rds(here("input", "ELISA_data.rds"))$ELISA_enriched %>%
  drop_na(trap_uid) %>%
  mutate(village = factor(village, levels = village_order, labels = str_to_sentence(village_order)),
         species = factor(clean_names, levels = str_replace_all(str_to_lower(all_species_order), " ", "_"), labels = all_species_order),
         # convert trap uid to first night as will be grouped
         trap_uid = paste0(str_split(trap_uid, pattern = "_", simplify = TRUE)[, 1], "_",
                           str_split(trap_uid, pattern = "_", simplify = TRUE)[, 2], "_",
                           1, "_",
                           str_split(trap_uid, pattern = "_", simplify = TRUE)[, 4], "_",
                           str_split(trap_uid, pattern = "_", simplify = TRUE)[, 5]),
         trap_uid = factor(trap_uid)) %>%
  select(rodent_uid, village, visit, trap_uid, species, genus, interpretation, sex, age_group, weight, head_body, hind_foot, ear, length_skull) %>%
  arrange(village, visit, rodent_uid) %>%
  mutate(rodent_uid = fct_inorder(rodent_uid))
  
trap_data <- read_rds(here("input", "trap_data.rds")) %>%
  tibble() %>%
  select(-geometry) %>%
  mutate(village = factor(village, levels = village_order, labels = str_to_sentence(village_order)),
         landuse = factor(landuse, levels = str_to_lower(names(landuse_palette)), labels = str_to_sentence(names(landuse_palette)))) %>%
  select(date_set, village, trap_uid, visit, grid_number, trap_number, landuse, trap_easting, trap_northing) %>%
  arrange(village, visit, grid_number, trap_number) %>%
  mutate(trap_uid = fct_inorder(trap_uid))

# Number of rodents trapped
nrow(unique_rodents)

# Number of trapnights
43266

# Number of species
unique_rodents %>%
  mutate(total = n()) %>%
  group_by(species) %>%
  mutate(n = n(),
         prop = n/total) %>%
  distinct(species, n, prop) %>%
  arrange(-prop)


# Describing communities --------------------------------------------------

left_join(unique_rodents, trap_data) %>%
  group_by(landuse, visit) %>%
  summarise(richness = length(unique(species)),
            n = n()) %>%
  group_by(landuse) %>%
  summarise(median_richness = median(richness),
            IQR_richness = IQR(richness),
            median_n = median(n),
            IQR_n = IQR(n))


# For running when using ELISA --------------------------------------------

# Number of antibody positive rodents
unique_rodents %>%
  filter(interpretation == "Positive") %>%
  distinct(rodent_uid) %>%
  nrow()

# Positive rodents by species
table_1_df <- unique_rodents %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0),
         total_positive = sum(n_positive)) %>%
  group_by(species) %>%
  summarise(n_individuals = n(),
            n_positive = sum(n_positive),
            total_positive = unique(total_positive)) %>%
  mutate(perc_positive = round((n_positive/n_individuals) * 100, 1),
         perc_all_positive = round((n_positive/total_positive) * 100, 1)) %>%
  arrange(-n_positive, -n_individuals) %>%
  mutate(Species = species,
         `Indviduals (N)` = n_individuals,
         `LASV Antibody detected (%)` = paste0(n_positive, " (", perc_positive, "%)"),
         `Percentage of all positive individuals` = paste0(perc_all_positive, "%")) %>%
  select(Species, `Indviduals (N)`, `LASV Antibody detected (%)`, `Percentage of all positive individuals`)

table_1_ft <- flextable(table_1_df) %>%
  italic(j = 1, part = "body")

write_rds(table_1_ft, here("output", "Table_1.rds"))

# Positive rodents by location and visit

rodent_detail <- unique_rodents %>%
  select(rodent_uid, trap_uid, species, interpretation) %>%
  left_join(trap_data %>%
              select(date_set, village, visit, trap_uid, landuse),
            by = c("trap_uid"))

rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                         TRUE ~ 0)) %>%
  group_by(village) %>%
  summarise(tested_village = n(),
         n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_village = round(n_positive/tested_village * 100, 1))

rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0)) %>%
  group_by(landuse) %>%
  summarise(tested_landuse = n(),
            n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_landuse = round(n_positive/tested_landuse * 100, 1))

rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0)) %>%
  left_join(visit_season, by = c("visit")) %>%
  group_by(season) %>%
  summarise(tested_season = n(),
            n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_season = round(n_positive/tested_season * 100, 1))

## Needs sorting when more ELISA data available
# rodent_detail %>%
#   group_by(visit) %>%
#   mutate(tested_visit = n(),
#          month = min(month(date_set), na.rm = TRUE)) %>%
#   group_by(visit, month, interpretation) %>%
#   summarise(n = n()) %>%
#   group_by(visit) %>%
#   mutate(prop = n/sum(n))
# 
# rodent_detail %>%
#   group_by(visit, interpretation) %>%
#   summarise(tested_visit = n(),
#             month = min(month(date_set), na.rm = TRUE)) %>%
#   mutate(interpretation = as.character(interpretation),
#          interpretation = replace_na(interpretation, "Not tested")) %>%
#   arrange(month, visit, interpretation)
