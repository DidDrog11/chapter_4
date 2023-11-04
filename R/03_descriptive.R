
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
         `Individuals (N)` = n_individuals,
         `LASV Antibody detected (%)` = paste0(n_positive, " (", perc_positive, "%)"),
         `Percentage of all positive individuals` = paste0(perc_all_positive, "%")) %>%
  select(Species, `Individuals (N)`, `LASV Antibody detected (%)`, `Percentage of all positive individuals`)

# Prevalence Odds Ratio (POR) by species 
por_df <- unique_rodents %>%
  mutate(seropositive = case_when(interpretation == "Positive" ~ 1,
                                  TRUE ~ 0)) %>%
  select(species, seropositive) %>%
  group_by(species) %>%
  filter(n() >= 10) %>%
  ungroup()

# Set "Mastomys natalensis" as the reference category
por_df$species <- relevel(por_df$species, ref = "Mastomys natalensis")
por_df$species <- droplevels(por_df$species)

# Define the prior probability levels
prior_mas = 0.3       # Prior probability for Mastomys natalensis > 10%
prior_non_mas = 0.05  # Prior probability for non-Mastomys, non-Crocidura species (1-10%)
prio_croc = NULL      # Non-informative prior for Crocidura species

# Define the Bayesian logistic regression model
model <- brm(
  seropositive ~ species,
  data = por_df,
  family = bernoulli(link = "logit"),
  prior = set_prior("normal(0, 1)", class = "b"),
  sample_prior = TRUE,
  iter = 2000,
  warmup = 1000,
  chains = 4
)

# Extract the posterior samples of the coefficients
post_samples <- as_draws_df(model, variable = "^b_species", regex = TRUE) %>%
  gather(Species, Estimate, starts_with("b_species"))

species_mapping <- c(
  "Crocidurabuettikoferi" = "Crocidura buettikoferi",
  "Crociduragrandiceps" = "Crocidura grandiceps",
  "Crociduraolivieri" = "Crocidura olivieri",
  "Lemniscomysstriatus" = "Lemniscomys striatus",
  "Lophuromyssikapusi" = "Lophuromys sikapusi",
  "Malacomysedwardsi" = "Malacomys edwardsi",
  "Musmusculus" = "Mus musculus",
  "Mussetulosus" = "Mus setulosus",
  "Praomysrostratus" = "Praomys rostratus",
  "Rattusrattus" = "Rattus rattus"
)

# Remove the "b_species" prefix and rename the species
summary_df <- post_samples %>%
  mutate(
    Species = gsub("^b_species", "", Species),
    Species = species_mapping[Species],
    Species = fct_rev(factor(Species, levels = species_order_n))
  )

# Calculate the mean and 95% credible interval
summary_stats <- summary_df %>%
  group_by(Species) %>%
  summarize(
    Mean = mean(Estimate),
    CI_low = quantile(Estimate, 0.025),
    CI_high = quantile(Estimate, 0.975)
  )

# Create a forest plot using ggplot
ggplot(summary_df, aes(x = Estimate, fill = Species)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ Species, ncol = 1) +
  labs(
    x = "Difference in Probability (vs. Mastomys natalensis)",
    y = "Density",
    title = "Density Plot of Difference in Probability of Seropositivity by Species (vs. Mastomys natalensis)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")

odds_seropositivity <- ggplot() +
  geom_density_ridges(data = summary_df, aes(x = Estimate, y = Species), scale = 1, fill = "lightblue", alpha = 0.6) +
  geom_point(data = summary_stats, aes(x = Mean, y = Species), size = 3, color = "black") +
  geom_errorbarh(data = summary_stats, aes(xmin = CI_low, xmax = CI_high, y = Species), height = 0.2, color = "black", size = 1) +
  labs(
    x = expression('Difference in Log-Odds vs.'~italic('Mastomys natalensis')~'(Reference)'),
    y = "Species",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "italic"),
  ) +
  scale_x_continuous(breaks = seq(-2, 2, 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(clip = 'off')

save_plot(filename = here("output", "seropositivity_OR_output.png"), plot = odds_seropositivity, base_width = 8)

# Add ORs to table
table_1_updated <- table_1_df %>%
  left_join(summary_stats %>%
              mutate(`OR (95% CrI)` = paste0(round(Mean, 2), " (", round(CI_low, 2), "-", round(CI_high, 2), ")")) %>%
              select(Species, `OR (95% CrI)`),
            by = "Species")
# Add "ref" to table
table_1_updated[1, 5] <- "ref"

write_rds(table_1_updated, here("output", "unformatted_Table_1.rds"))

table_1_ft <- flextable(table_1_updated) %>%
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

season_tab <- rodent_detail %>%
  mutate(n_positive = case_when(interpretation == "Positive" ~ 1,
                                TRUE ~ 0)) %>%
  left_join(visit_season, by = c("visit")) %>%
  group_by(season) %>%
  summarise(tested_season = n(),
            n_positive = sum(n_positive)) %>%
  mutate(total_positive = sum(n_positive),
         prop_positive = round(n_positive/total_positive * 100, 1),
         prop_season = round(n_positive/tested_season * 100, 1))

season_tab_chi <- as.table(rbind(c(16, 444-16),
                                 c(23, 240-23)))
dimnames(season_tab_chi) <- list(season = c("Dry", "Rainy"),
                                 serostatus = c("Positive", "Negative"))

(seasonal <- chisq.test(season_tab_chi))

# Figure 1 ----------------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

africa <- world %>%
  filter(str_detect(continent, "Africa")) %>%
  select(admin) %>%
  mutate(fill = case_when(str_detect(admin, "Sierra Leone") ~ "black",
                          TRUE ~ "grey")) %>%
  ms_filter_islands(min_area = 1e10)

sl_bbox <- africa %>%
  filter(str_detect(admin, "Sierra Leone")) %>%
  st_bbox() %>%
  st_as_sfc()

africa_map <- ggplot() +
  geom_sf(data = africa, aes(fill = fill)) +
  geom_sf(data = sl_bbox, fill = NA, colour = "black", size = 4) +
  scale_fill_manual(values = c("grey", "white")) +
  guides(fill = "none") +
  theme_void()

sle_sf <- geodata::gadm(country = "SLE", level = 2, path = here("data", "geodata")) %>%
  st_as_sf()

fig_1_palette <- c(village_palette, "#FFFFFF")
names(fig_1_palette) <- c(names(village_palette)[1:4], "poi")

poi <- tibble(name = c("Kenema", "Baiama", "Lalehun", "Lambayama", "Seilama"),
              cat = c("poi", "Baiama", "Lalehun", "Lambayama", "Seilama"),
              lat = c(7.876161956810467, 7.837529372181356, 8.197392257077409, 7.850593096948891, 8.12230048178563),
              lon = c(-11.190811585001954, -11.268407665149846, -11.08032958100431, -11.196939025872055, -11.1935976318381)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = default_CRS)


# Focus map on Kenema
kenema_map <- get_googlemap(st_coordinates(poi %>% filter(name == "Kenema"))[1, ], zoom = 9, maptype = "satellite", scale = 2) %>%
  rast() %>%
  mask(vect(sle_sf %>%
              filter(str_detect(NAME_2, "Kenema"))))

zoom_kenema <- sle_sf %>%
  filter(str_detect(NAME_2, "Kenema")) %>% st_bbox()

study_map <- ggplot() + 
  geom_sf(data = sle_sf, fill = "grey", colour = "black") +
  geom_spatraster_rgb(data = kenema_map) +
  geom_sf(data = sle_sf %>%
            filter(str_detect(NAME_2, "Kenema")), fill = NA, colour = "white") +
  geom_sf(data = poi, size = 0.8) +
  ggrepel::geom_label_repel(data = poi, aes(label = name, geometry = geometry, fill = cat), stat = "sf_coordinates", min.segment.length = 0) +
  labs(x = element_blank(),
       y = element_blank()) +
  scale_fill_manual(values = fig_1_palette) +
  theme_bw() +
  annotation_north_arrow(style = north_arrow_minimal()) +
  annotation_scale(location = "br") +
  coord_sf(xlim = c(zoom_kenema[1], zoom_kenema[3]), ylim = c(zoom_kenema[2], zoom_kenema[4])) +
  guides(fill = "none")

sl_inset_map <- ggdraw() +
  draw_plot(study_map) +
  draw_plot(africa_map, x = 0.7, y = 0.7, width = 0.3, height = 0.3)

ggsave2(plot = sl_inset_map, filename = here("output", "Figure_1.png"), dpi = 300, width = 8, height = 6)
