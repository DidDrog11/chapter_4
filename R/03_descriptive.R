# General descriptives
unique_rodents <- read_rds(here("data", "unique_rodents.rds"))


# Using chapter 3 extract -------------------------------------------------
unique_rodents <- read_rds(here("data", "chapter_4_rodents.rds"))
trap_data <- read_rds(here("data", "chapter_4_traps.rds"))

# Number of rodents trapped
nrow(unique_rodents)

# Number of trapnights
nrow(trap_data) * 4

# Number of species
unique_rodents %>%
  mutate(total = n()) %>%
  group_by(species) %>%
  mutate(n = n(),
         prop = n/total) %>%
  distinct(species, n, prop)


# For running when using ELISA --------------------------------------------

# # Number of antibody positive rodents
# unique_rodents %>%
#   filter(interpretation == "Positive") %>%
#   distinct(rodent_uid) %>%
#   nrow()
# 
# # Positive rodents by species
# unique_rodents %>%
#   group_by(initial_species_id) %>%
#   mutate(tested = n()) %>%
#   filter(interpretation == "Positive") %>%
#   distinct(rodent_uid, initial_species_id, tested) %>%
#   mutate(n = n(),
#          total_positive = nrow(.)) %>%
#   distinct(initial_species_id, n, total_positive, tested) %>%
#   mutate(prop_of_positive = n/total_positive,
#          rate_of_positive = n/tested)
# 
# # Positive rodents by location and visit
# unique_rodents %>%
#   group_by(village) %>%
#   mutate(tested_village = n()) %>%
#   filter(interpretation == "Positive") %>%
#   group_by(village, tested_village) %>%
#   summarise(n_positive = n()) %>%
#   ungroup() %>%
#   mutate(total_positive = sum(n_positive),
#          prop_positive = n_positive/total_positive,
#          prop_village = n_positive/tested_village)
# 
# unique_rodents %>%
#   group_by(visit) %>%
#   mutate(tested_visit = n(),
#          month = min(month_set, na.rm = TRUE)) %>%
#   group_by(visit, month, interpretation) %>%
#   summarise(n = n()) %>%
#   group_by(visit) %>%
#   mutate(prop = n/sum(n))
# 
# unique_rodents %>%
#   group_by(visit) %>%
#   mutate(tested_visit = n(),
#          month = min(month_set, na.rm = TRUE)) %>%
#   group_by(visit, month, interpretation) %>%
#   summarise(n = n()) %>%
#   group_by(month, interpretation) %>%
#   mutate(n = sum(n)) %>%
#   pivot_wider(names_from = interpretation,
#               values_from = n) %>%
#   rowwise() %>%
#   mutate(tested = sum(Positive, Negative, Equivocal, na.rm = TRUE)) %>%
#   group_by(month) %>%
#   mutate(positive = sum(Positive),
#          tested = sum(tested)) %>%
#   distinct(month, tested, positive) %>%
#   rowwise() %>%
#   mutate(prop_positive = positive/tested) %>%
#   ggplot() +
#   geom_col(aes(x = month, y = positive)) +
#   geom_point(aes(x = month, y = prop_positive*100)) +
#   geom_line(aes(x = month, y = prop_positive*100)) +
#   scale_y_continuous(limits = c(0, 15), "Number positive",
#                      sec.axis = sec_axis(~ . / 1.2, name = "Percentage positive")) +
#   labs(x = "Month") +
#   theme_bw()
#   
# # Positive rodents by habitat type
# unique_rodents %>%
#   group_by(simple_habitat) %>%
#   mutate(n_tested = n()) %>%
#   filter(interpretation == "Positive") %>%
#   summarise(n_positive = n(),
#             prop_positive = n()/n_tested) %>%
#   distinct()
