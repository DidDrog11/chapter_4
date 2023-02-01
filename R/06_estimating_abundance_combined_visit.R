source(here("R", "00_setup.R"))

unique_rodents <- read_rds(here("data", "chapter_4_rodents.rds"))
trap_data <- read_rds(here("data", "chapter_4_traps.rds"))

# The site will be landuse type, the visit will be the replicate

# Rodent species with less than 10 individuals will be grouped as Other
species_n <- unique_rodents %>%
  group_by(species) %>%
  summarise(N = n())

other_species <- species_n$species[species_n$N <= 10]

rodents <- unique_rodents %>%
  select(rodent_id, village, visit, trap_id, species) %>%
  mutate(species = case_when(species %in% other_species ~ "other_spp",
                             TRUE ~ species)) %>%
  left_join(trap_data %>%
              select(trap_id, site_id, landuse),
            by = c("trap_id")) %>%
  group_by(rodent_id) %>% 
  arrange(site_id) %>%
  slice(1) %>%
  select(-trap_id)

sites <- trap_data %>%
  distinct(site_id, trap_id, village, visit, grid_number, landuse) %>%
  group_by(landuse, village, visit) %>%
  mutate(tn = 4 * n()) %>%
  group_by(site_id, village, visit, grid_number, landuse) %>%
  select(-trap_id) %>%
  group_by(village, landuse) %>%
  mutate(site_number = cur_group_id()) %>%
  ungroup() %>%
  distinct(site_number, visit, tn) %>%
  arrange(visit) %>%
  pivot_wider(names_from = visit, names_prefix = "visit_", values_from = tn) %>%
  ungroup() %>%
  relocate(site_number) %>%
  arrange(site_number)

site_numbers <- trap_data %>%
  distinct(site_id, trap_id, village, visit, grid_number, landuse) %>%
  group_by(landuse, village, visit) %>%
  mutate(tn = 4 * n()) %>%
  group_by(site_id, village, visit, grid_number, landuse) %>%
  select(-trap_id) %>%
  group_by(village, landuse) %>%
  mutate(site_number = cur_group_id()) %>%
  ungroup() %>%
  select(site_number, site_id) %>%
  distinct()

surveyed_site <- trap_data %>%
  distinct(site_id, visit) %>%
  left_join(site_numbers) %>%
  select(site_number, visit) %>%
  arrange(site_number, visit) %>%
  distinct()

site_covs <- sites %>%
  left_join(site_numbers %>%
              mutate(village = str_split(site_id, "_", simplify = TRUE)[, 1]) %>%
              distinct(site_number, village), by = c("site_number"))  %>%
  mutate(location = case_when(str_detect(village, "lambayama") ~ "peri-urban",
                              TRUE ~ "rural")) %>%
  left_join(trap_data %>%
              select(site_id, landuse) %>%
              left_join(site_numbers) %>%
              distinct(landuse, site_number)) %>%
  select(site_number, location, landuse) %>%
  arrange(site_number) %>%
  distinct() %>%
  select(-site_number) %>%
  as.data.frame()

final_site_df <- site_numbers %>%
  left_join(trap_data %>%
              distinct(site_id, landuse, village, grid_number)) %>%
  select(site_number, village, landuse) %>%
  arrange(site_number) %>%
  distinct()

obs_covs <- list(tn = sites %>%
                   select(-site_number) %>%
                   replace(., . == 0, NA) %>%
                   as.matrix(),
                 season = trap_data %>%
                   left_join(visit_season, by = "visit") %>%
                   distinct(site_id, visit, season) %>%
                   left_join(site_numbers) %>%
                   arrange(visit) %>%
                   distinct(visit, season, site_number) %>%
                   pivot_wider(names_from = visit, names_prefix = "visit_", values_from = season) %>%
                   ungroup() %>%
                   arrange(site_number) %>%
                   select(-site_number) %>%
                   as.matrix())

rodents_ab <- rodents %>%
  left_join(site_numbers, by = "site_id") %>%
  ungroup() %>%
  arrange(species) %>%
  group_by(site_number, visit, species) %>%
  summarise(n = n()) %>%
  group_by(species) %>%
  group_split()

names(rodents_ab) <- unique(sort(rodents$species))

rodents_ab_binary <- lapply(rodents_ab, function(x) {
  
  left_join(surveyed_site, x) %>%
    mutate(species = case_when(is.na(species) ~ 0,
                               TRUE ~ 1)) %>%
    mutate(species = coalesce(n, species)) %>%
    group_by(site_number, visit) %>%
    mutate(n = sum(species)) %>%
    arrange(visit) %>%
    distinct(site_number, visit, n) %>%
    pivot_wider(names_from = visit, names_prefix = "visit_", values_from = n) %>%
    arrange(site_number)
  
})

names(rodents_ab_binary) <- names(rodents_ab)

# This doesn't make much sense to do for other species, probably best to run for all individually and then group

# Create a table to manage the K input into the N-mixture model
k_table = tibble(species = names(rodents_ab),
                 poisson_k = c(120, 95, 1200, 95, 120, 160, 300, 500, 1200),
                 nb_k = c(60, 40, 700, 75, 800, 300, 500, 350, 950),
                 zip_k = c(155, 85, 700, 40, 130, 230, 290, 220, 800))

# Abundance function ------------------------------------------------------
produce_abundance <- function(binary_list = rodents_ab_binary, species_name = "mus_musculus") {
  
  ab <- binary_list[[species_name]] %>%
    ungroup() %>%
    select(-site_number) %>%
    as.matrix()
  
  umf <- unmarkedFramePCount(y = ab, siteCovs = site_covs, obsCovs = obs_covs)
  
  ab_p <- pcount(formula = ~ tn + season ~ location + landuse, umf, mixture = "P", K = k_table$poisson_k[k_table$species == species_name])
  ab_nb <- pcount(formula = ~ tn + season ~ location + landuse, umf, mixture = "NB", K = k_table$nb_k[k_table$species == species_name])
  ab_zip <- pcount(formula = ~ tn + season ~ location + landuse, umf, mixture = "ZIP", K = k_table$zip_k[k_table$species == species_name])
  
  ab_list <- list("poisson" = ab_p,
                  "neg_binomial" = ab_nb,
                  "zero_inflated_poisson" = ab_zip)
  
  # Identify which model has the lowest AIC
  aic <- which.min(c(ab_p@AIC, ab_nb@AIC, ab_zip@AIC))
  
  # Select this model
  ab_selected <- ab_list[[aic]]
  
  lowest_aic <- names(ab_list[aic])
  
  # Estimate the posterior distributions of the random variables
  ab_ranef <- ranef(ab_selected)
  
  observed_n <- sum(ab, na.rm = TRUE)
  estimated_n <- sum(bup(ab_ranef))
  
  ab_df <- show(ab_ranef) %>%
    as.data.frame() %>%
    mutate(species = species_name,
           site_number = row_number()) %>%
    relocate(site_number, species) %>%
    left_join(final_site_df) %>%
    select(site_number, village, landuse, Mean, Mode) %>%
    tibble()
  
  colnames(ab_df) <- c(colnames(ab_df[1:3]), paste0(species_name, "_", "mean"), paste0(species_name, "_", "mode"))
  
  return(list(model_selected = lowest_aic,
              model = ab_selected,
              comparison = tibble("observed" = observed_n,
                                  "estimated" = estimated_n),
              output = ab_df))
  
}

crocidura <- produce_abundance(species_name = "crocidura_spp")
lophuromys <- produce_abundance(species_name = "lophuromys_sikapusi")
mastomys_natalensis <- produce_abundance(species_name = "mastomys_natalensis")
mus_minutoides <- produce_abundance(species_name = "mus_minutoides")
mus_musculus <- produce_abundance(species_name = "mus_musculus")
other_spp <- produce_abundance(species_name = "other_spp")
praomys_spp <- produce_abundance(species_name = "praomys_spp")
rattus_rattus <- produce_abundance(species_name = "rattus_rattus")

combined_abundance <- left_join(crocidura$output, lophuromys$output, by = c("site_number", "village", "landuse")) %>%
  left_join(mastomys_natalensis$output, by = c("site_number", "village", "landuse")) %>%
  left_join(mus_minutoides$output, by = c("site_number", "village", "landuse")) %>%
  left_join(mus_musculus$output, by = c("site_number", "village", "landuse")) %>%
  left_join(other_spp$output, by = c("site_number", "village", "landuse")) %>%
  left_join(praomys_spp$output, by = c("site_number", "village", "landuse")) %>%
  left_join(rattus_rattus$output, by = c("site_number", "village", "landuse")) %>%
  select(-contains("_mean"))

ab_df <- combined_abundance %>%
  pivot_longer(cols = contains("mode"), names_to = "species", values_to = "predicted_n") %>%
  mutate(species = factor(str_remove(species, "_mode")),
         landuse = factor(str_to_sentence(landuse), levels = c("Forest", "Agriculture", "Village")),
         village = str_to_sentence(village))

ggplot(ab_df) +
  geom_col(aes(y = predicted_n, x = species, fill = landuse, group = landuse), position = position_dodge()) +
  geom_segment(data = rodents %>%
                 group_by(species, landuse) %>%
                 mutate(observed = n(),
                        landuse = factor(str_to_sentence(landuse), levels = c("Forest", "Agriculture", "Village")),
                        species = factor(species, levels = levels(ab_df$species))) %>%
                 distinct(species, landuse, observed),
               aes(y = observed, yend = observed, x = as.numeric(species) - 0.45, xend = as.numeric(species) + 0.45, linetype = landuse), position = position_dodge()) +
  theme_bw()

write_rds(combined_abundance, here("data", "estimated_population.rds"))

