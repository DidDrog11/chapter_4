source(here::here("R", "00_setup.R"))

set.seed(123)

df <- read_rds(here("data", "data_for_abundance.rds"))

rodents <- df$rodents

trap_data <- df$trap_data

network_numbers <- df$network_numbers

# Rodent species with less than 20 individuals will be grouped as Other
species_n <- rodents %>%
  group_by(species) %>%
  summarise(N = n())

other_species <- species_n$species[species_n$N < 20]

sites <- trap_data %>%
  distinct(trap_uid, village, visit, grid, landuse) %>%
  mutate(tn = 1) %>%
  left_join(network_numbers, by = c("village" = "Village", "visit" = "Visit", "grid" = "Grid", "landuse" = "Landuse")) %>%
  group_by(Network, visit, village, grid, landuse) %>%
  summarise(tn = sum(tn)) %>%
  drop_na(Network) %>%
  distinct(Network, visit, village, grid, landuse, tn) %>%
  left_join(visit_season %>%
              mutate(visit = factor(visit))) %>%
  arrange(Network, village, grid, visit) %>%
  mutate(location = dt_case_when(str_detect(village, "lambayama") ~ "peri-urban",
                                 TRUE ~ "rural")) %>%
  group_by(village, grid) %>%
  mutate(site_id = cur_group_id())

site_id <- sites %>%
  distinct(village, grid, landuse, site_id)

surveyed_site <- trap_data %>%
  distinct(trap_uid, village, visit, grid, landuse) %>%
  left_join(network_numbers, by = c("village" = "Village", "visit" = "Visit", "grid" = "Grid", "landuse" = "Landuse")) %>%
  distinct(village, visit, grid, landuse, Network) %>%
  arrange(visit, village, grid) 

site_covs <- sites %>%
  select(-Network) %>%
  arrange(site_id, visit, village) %>%
  as.data.frame() 

obs_df <- sites %>%
  group_by(site_id, visit, village) %>%
  summarise(tn = sum(tn),
            season = unique(season),
            village = unique(village)) %>%
  arrange(visit, site_id)

obs_covs <- list(tn = obs_df %>%
                   select(-season) %>%
                   pivot_wider(names_from = visit, names_prefix = "visit_", values_from = tn) %>%
                   arrange(site_id) %>%
                   ungroup() %>%
                   select(-site_id, -village) %>%
                   as.matrix(),
                 season = obs_df %>%
                   select(-tn) %>%
                   pivot_wider(names_from = visit, names_prefix = "visit_", values_from = season) %>%
                   arrange(site_id) %>%
                   ungroup() %>%
                   select(-site_id, -village) %>%
                   as.matrix())

rodents_ab <- rodents %>%
  left_join(site_id, by = c("village", "grid", "landuse")) %>%
  mutate(species = dt_case_when(species %in% other_species ~ "other_spp",
                                TRUE ~ species)) %>%
  ungroup() %>%
  arrange(species) %>%
  select(site_id, visit, species) %>%
  group_by(species) %>%
  group_split()

names(rodents_ab) <- bind_rows(rodents_ab) %>%
  arrange(species) %>%
  distinct(species) %>%
  pull(species)

rodents_ab_binary <- lapply(rodents_ab, function(x) {
  
  left_join(site_covs, x) %>%
    mutate(species = case_when(is.na(species) ~ 0,
                               TRUE ~ 1)) %>%
    group_by(site_id, visit) %>%
    mutate(n = sum(species)) %>%
    arrange(visit) %>%
    distinct(site_id, visit, n) %>%
    pivot_wider(names_from = visit, names_prefix = "visit_", values_from = n) %>%
    arrange(site_id)
  
})

names(rodents_ab_binary) <- names(rodents_ab)

# Abundance function ------------------------------------------------------
produce_abundance <- function(binary_list = rodents_ab_binary, species_name = "crocidura_spp", site_covs_df = site_covs, obs_covs_df = obs_covs, k = k) {
  
  ab <- binary_list[[species_name]] %>%
    ungroup() %>%
    select(-site_id) %>%
    as.matrix()
  
  site_covs_summarise = site_covs_df %>%
    distinct(site_id, location, landuse) %>%
    select(-site_id) %>%
    as.data.frame()
  
  umf <- unmarkedFramePCount(y = ab, siteCovs = site_covs_summarise, obsCovs = obs_covs_df)
  
  ab_p <- pcount(formula = ~ tn + season ~ location + landuse, umf, mixture = "P", K = k)
  ab_nb <- pcount(formula = ~ tn + season ~ location + landuse, umf, mixture = "NB", K = k)
  ab_zip <- pcount(formula = ~ tn + season ~ location + landuse, umf, mixture = "ZIP", K = k)
  
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
  ppd <- posteriorSamples(ab_ranef, nsims = 1000)
  
  median_ab <- vector()
  
  for(i in 1:dim(ppd@samples)[1]) {
  
    median_ab[i] <- median(ppd@samples[i, ,])
    
  }
  
  observed_n <- sum(ab, na.rm = TRUE)
  estimated_n <- sum(bup(ab_ranef))
  
  ab_df <- show(ab_ranef) %>%
    as.data.frame() %>%
    mutate(species = species_name,
           site_id = row_number()) %>%
    relocate(site_id, species) %>%
    left_join(site_id) %>%
    mutate(Median = median_ab) %>%
    select(site_id, grid, village, landuse, Mean, Mode, Median) %>%
    mutate(observed = rowSums(binary_list[[species_name]] %>%
                                ungroup() %>%
                                select(-site_id), na.rm = TRUE)) %>%
  tibble()
  
  colnames(ab_df) <- c(colnames(ab_df[1:4]), paste0(species_name, "_", "mean"),
                       paste0(species_name, "_", "mode"),
                       paste0(species_name, "_", "median"),
                       paste0(species_name, "_", "observed"))
  
  return(list(model_selected = lowest_aic,
              model = ab_selected,
              comparison = tibble("observed" = observed_n,
                                  "estimated" = estimated_n),
              output = ab_df))
  
}

crocidura <- produce_abundance(species_name = "crocidura_spp", k = 60)
lophuromys <- produce_abundance(species_name = "lophuromys_sikapusi", k = 60)
mastomys_natalensis <- produce_abundance(species_name = "mastomys_natalensis", k = 80)
mus_musculus <- produce_abundance(species_name = "mus_musculus", k = 120)
mus_setulosus <- produce_abundance(species_name = "mus_setulosus", k = 50)
other_spp <- produce_abundance(species_name = "other_spp", k = 60)
praomys_spp <- produce_abundance(species_name = "praomys_rostratus", k = 120)
rattus_rattus <- produce_abundance(species_name = "rattus_rattus", k = 60)

combined_abundance <- left_join(crocidura$output, mastomys_natalensis$output, by = c("site_id", "grid", "village", "landuse")) %>%
  left_join(lophuromys$output, by = c("site_id", "grid", "village", "landuse")) %>%
  left_join(mus_musculus$output, by = c("site_id", "grid", "village", "landuse")) %>%
  left_join(mus_setulosus$output, by = c("site_id", "grid", "village", "landuse")) %>%
  left_join(other_spp$output, by = c("site_id", "grid", "village", "landuse")) %>%
  left_join(praomys_spp$output, by = c("site_id", "grid", "village", "landuse")) %>%
  left_join(rattus_rattus$output, by = c("site_id", "grid", "village", "landuse")) %>%
  select(-matches("_mean|_mode|_observed")) %>%
  mutate(across(.cols = matches("_median"), \(x) round(x, digits = 0)))

estimated_population <- network_nodelist %>%
  left_join(combined_abundance %>%
              mutate(Grid = grid,
                     Village = str_to_sentence(village),
                     Landuse = str_to_sentence(landuse)) %>%
              select(-site_id, -grid, -village, -landuse)) %>%
  group_by(Network, Visit, Village) %>% 
  summarise(across(.cols = matches("_median"), \(x) sum(x, na.rm = FALSE)))

plot_abundance <- function(species, fig_label = "A", label_y_position = 250, label_x_shift = 5) {
  
  ppd <- posteriorSamples(ranef(species$model), nsims = 1000)
  
  site_landuse_number <- site_covs %>%
    distinct(site_id, landuse)
  
  site_village_number <- site_id %>% 
    ungroup() %>%
    select(village, site_id) %>% 
    distinct()
  
  site_label <- site_landuse_number %>%
    left_join(site_village_number, by = c("site_id")) %>%
    mutate(label = str_to_title(paste0(site_id, "-", village, ": ", landuse)),
           landuse = str_to_title(landuse)) %>%
    select(site_id, label, landuse)
  
  sample_list <- list()
  
  for(i in 1:dim(ppd@samples)[1]) {
    
    s <- ppd@samples[i, , ]
    
    sample_list[[i]] <- tibble(abundance = s) %>%
      bind_cols(site_label[i, ])
    
  }
  
  sample_df <- sample_list %>%
    bind_rows() %>%
    arrange(site_id) %>%
    mutate(site_id = factor(site_id, labels = unique(label)))
  
  median_df <- sample_df %>%
    group_by(site_id) %>%
    mutate(median_abundance = round(median(abundance), 0)) %>%
    distinct(site_id, median_abundance)
  
  species_name = str_to_sentence(str_replace(str_remove_all(names(species$output)[str_detect(names(species$output), "_median")], "_median"), "_", " "))
  
  species_plot <- ggplot() +
    geom_bar(data = sample_df, aes(x = abundance, fill = landuse)) +
    geom_vline(data = median_df, aes(xintercept = median_abundance), linetype = "dashed") +
    geom_text(data = median_df, aes(x = median_abundance + label_x_shift, label = median_abundance), y = label_y_position, size = 4) +
    facet_wrap(~ site_id) +
    theme_bw() +
    scale_fill_manual(values = landuse_palette) +
    labs(title = paste0("S2", fig_label, ") ", species_name, " estimated abundance"),
         x = "Estimated abundance",
         y = "Posterior samples",
         fill = "Landuse")
  
  plot_name <- paste0("Supplementary_Figure_2", fig_label, ".png")
  
  save_plot(here("output", plot_name), plot = species_plot, base_height = 8)
  
}


plot_abundance(species = mastomys_natalensis, fig_label = "A", label_y_position = 250, label_x_shift = 5)
plot_abundance(species = mus_musculus, fig_label = "B", label_y_position = 400, label_x_shift = 5)
plot_abundance(species = rattus_rattus, fig_label = "C", label_y_position = 250, label_x_shift = 5)
plot_abundance(species = crocidura, fig_label = "D", label_y_position = 100, label_x_shift = 5)
plot_abundance(species = lophuromys, fig_label = "E", label_y_position = 100, label_x_shift = 5)
plot_abundance(species = mus_setulosus, fig_label = "F", label_y_position = 100, label_x_shift = 5)
plot_abundance(species = praomys_spp, fig_label = "G", label_y_position = 250, label_x_shift = 5)
plot_abundance(species = other_spp, fig_label = "H", label_y_position = 250, label_x_shift = 5)

write_rds(estimated_population, here("data", "estimated_population.rds"))

