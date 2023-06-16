source(here::here("R", "00_setup.R"))
source(here("R", "01_load_data.R"))
source(here("R", "02_clean_data.R"))
source(here("R", "03_descriptive.R"))

rodent_species <- unique_rodents %>%
  filter(!str_detect(species, "spp")) %>%
  distinct(species) %>%
  mutate(species = str_to_sentence(str_replace_all(species, "_", " "))) %>%
  pull(species)

rodent_genera <- unique_rodents %>%
  distinct(species) %>%
  filter(!species %in% c("rattus_rattus", "mus_musculus", "mus_setulosus", "lemniscomys_striatus", "mastomys_natalensis")) %>%
  mutate(genus = str_to_sentence(str_split(species, "_", simplify = TRUE)[, 1])) %>%
  pull(genus)

# Load HomeRange data
HomeRangeData <- GetHomeRangeData()

rodent_ranges <- HomeRangeData %>%
  mutate(genus = str_split(Species, " ", simplify = TRUE)[, 1]) %>%
  mutate(match = case_when(str_detect(Species, paste(rodent_species, collapse = "|")) ~ "Species",
                           genus %in% rodent_genera ~ "Genus",
                           TRUE ~ as.character(NA))) %>%
  filter(!is.na(match)) %>%
  select(Species, Home_Range_km2, HR_Spread_km2, Spread_Units, No_Individuals, Latitude, Longitude, Country, match) %>%
  mutate(Home_Range_m2 = Home_Range_km2 * 1000000,
         Home_Range_radius = sqrt(Home_Range_m2/pi),
         individuals = as.numeric(str_remove_all(No_Individuals, "[^0-9.]")),
         individuals = replace_na(individuals, 1)) %>%
  uncount(individuals)

n_species_hr <- rodent_ranges %>%
  filter(match == "Species") %>%
  distinct(Species)
# 4 species have home ranges, only two of these include data from Africa, for m natalensis all data is originated from Tanzania

s_1 <- rodent_ranges %>%
  ggplot() +
  stat_ecdf(aes(x = Home_Range_radius,
                colour = Species),
            lwd = 0.4) +
  geom_vline(aes(xintercept = 30), colour = "black", lwd = 0.2, linetype = "dashed") +
  facet_wrap(~ match) +
  theme_bw() +
  labs(y = "ECDF",
       x = "Home range radius (m)") +
  coord_cartesian(xlim = c(0, 100))

save_plot(filename = here("output", "Supplementary_Figure_1.png"), s_1)

prop_30 <- rodent_ranges %>% 
  group_by(Species) %>% 
  mutate(ecdfFun = list(ecdf(Home_Range_radius))) %>%
  distinct(Species, ecdfFun) %>%
  mutate(prop_30 = round(unlist(lapply(ecdfFun, function(x) {x(30)})) * 100, 4))
