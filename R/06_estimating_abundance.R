install.packages("unmarked")
library(unmarked)

rodents <- read_rds(here("data", "unique_rodents.rds")) %>%
  filter(village != "bambawo") %>%
  drop_na(initial_species_id) %>%
  drop_na(grid_number) %>%
  mutate(grid_number = case_when(grid_number %in% c("6", "7") ~ "6",
                                 TRUE ~ as.character(grid_number)))

sites <- tibble(village = c(rep("baiama", 19), rep("lalehun", 36), rep("lambayama", 19), rep("seilama", 36)),
                grid_number = as.character(c(rep(1:4, each = 1), rep(rep(c(1:4, 6)), times = 3), rep(1:6, times = 6), rep(1:4, each = 1), rep(rep(c(1:4, 6)), times = 3), rep(1:6, times = 6))),
                site = paste0(village, "_", grid_number),
                visit = c(rep(3, each = 4), rep(4:6, each = 5), rep(1:6, each = 6), rep(3, each = 4), rep(4:6, each = 5), rep(1:6, each = 6)))

site_numbers <- tibble(site = unique(sites$site)) %>%
                         mutate(site_number = row_number())

site_covs <- trap_locations %>%
  mutate(grid_number = case_when(grid_number %in% c("6", "7") ~ "6",
                                 TRUE ~ as.character(grid_number))) %>%
  select(date_set, village, grid_number, visit, simple_habitat) %>%
  distinct() %>%
  mutate(habitat = factor(case_when(simple_habitat == "secondary_forest" ~ "forest",
                                    is.na(simple_habitat) ~ "agriculture",
                                    TRUE ~ simple_habitat)),
         location = factor(case_when(village == "lambayama" ~ "urban",
                                     TRUE ~ "rural")),
         site = paste0(village, "_", grid_number)) %>%
  left_join(visit_season, by = "visit") %>%
  select(site, habitat, location, season) %>%
  distinct() %>%
  left_join(site_numbers) %>%
  select(site_number, location, habitat, season) %>%
  as.data.frame()

long_sites <- sites %>%
  left_join(site_numbers) %>%
  select(site_number, visit) %>%
  filter(!(site_number == 9 & visit == 4)) %>%
  filter(!(site_number == 10 & visit == 4)) %>%
  filter(!(site_number == 21 & visit == 4))

trap_nights <- trap_locations %>%
  filter(visit != 7) %>%
  mutate(grid_number = case_when(grid_number %in% c("6", "7") ~ "6",
                                 TRUE ~ as.character(grid_number)),
         visit = case_when(village == "baiama" & visit == 1 ~ 3,
                           village == "baiama" & visit == 2 ~ 4,
                           village == "baiama" & visit == 3 ~ 5,
                           village == "baiama" & visit == 4 ~ 6,
                           village == "lambayama" & visit == 1 ~ 3,
                           village == "lambayama" & visit == 2 ~ 4,
                           village == "lambayama" & visit == 3 ~ 5,
                           village == "lambayama" & visit == 4 ~ 6,
                           TRUE ~ visit)) %>%
  select(village, grid_number, trap_number, visit) %>%
  mutate(site = paste0(village, "_", grid_number)) %>%
  left_join(site_numbers, .) %>%
  group_by(site_number, visit) %>%
  summarise(tn = n() * 4)

season <- trap_locations %>%
  filter(visit != 7) %>%
  mutate(grid_number = case_when(grid_number %in% c("6", "7") ~ "6",
                                 TRUE ~ as.character(grid_number)),
         visit = case_when(village == "baiama" & visit == 1 ~ 3,
                           village == "baiama" & visit == 2 ~ 4,
                           village == "baiama" & visit == 3 ~ 5,
                           village == "baiama" & visit == 4 ~ 6,
                           village == "lambayama" & visit == 1 ~ 3,
                           village == "lambayama" & visit == 2 ~ 4,
                           village == "lambayama" & visit == 3 ~ 5,
                           village == "lambayama" & visit == 4 ~ 6,
                           TRUE ~ visit)) %>%
  distinct(village, grid_number, visit) %>%
  mutate(site = paste0(village, "_", grid_number)) %>%
  left_join(site_numbers, .) %>%
  mutate(season = case_when(visit %in% c("1", "2", "5", "6") ~ "dry",
                            TRUE ~ "wet")) %>%
  select(site_number, visit, season) %>%
  mutate(visit = paste0("x", visit)) %>%
  arrange(visit) %>%
  pivot_wider(names_from = visit, values_from = season)

tn <- long_sites %>%
  left_join(., trap_nights, by = c("site_number", "visit")) %>%
  mutate(visit = paste0("x", visit)) %>%
  arrange(visit) %>%
  pivot_wider(names_from = visit, values_from = tn)

obs_covs <- list(tn = as.matrix(tn %>%
                                 arrange(site_number) %>%
                                 select(-site_number)),
                season = as.matrix(season %>%
                                     arrange(site_number) %>%
                                     select(-site_number)))

rodents_ab <- read_rds(here("data", "unique_rodents.rds")) %>%
  filter(village != "bambawo") %>%
  drop_na(initial_species_id) %>%
  drop_na(grid_number) %>%
  mutate(grid_number = case_when(grid_number %in% c("6", "7") ~ "6",
                                 TRUE ~ as.character(grid_number))) %>%
  select(initial_species_id, village, grid_number, visit) %>%
  mutate(grid_number = as.numeric(grid_number),
         visit = as.numeric(as.character(visit))) %>%
  mutate(site = str_c(village, grid_number, sep = "_")) %>%
  group_by(site, visit, initial_species_id) %>%
  summarise(n = n()) %>%
  left_join(site_numbers)


# Crocidura abundance -----------------------------------------------------

croc_ab <- rodents_ab %>%
  ungroup() %>%
  filter(str_detect(initial_species_id, "crocidura")) %>%
  select(site_number, visit, n) %>%
  right_join(long_sites) %>%
  mutate(visit = paste0("x", visit),
         n = replace_na(n, 0)) %>%
  arrange(visit) %>%
  pivot_wider(names_from = visit, values_from = n) %>%
  arrange(site_number) %>%
  select(-site_number) %>%
  as.matrix()

crocidura_umf <- unmarkedFramePCount(y = croc_ab, siteCovs = site_covs, obsCovs = obs_covs)
summary(crocidura_umf)

crocidura_ab_p <- pcount(formula = ~tn + season ~ location + habitat, crocidura_umf, mixture = "P", K = 105)
crocidura_ab_nb <- pcount(formula = ~tn + season ~ location + habitat, crocidura_umf, mixture = "NB", K = 105)
crocidura_ab_zip <- pcount(formula = ~tn + season ~ location + habitat, crocidura_umf, mixture = "ZIP", K = 105)

crocidura_ab_nb_re <- ranef(crocidura_ab_nb)
plot(crocidura_ab_nb_re, subset=site %in% 1:22)

sum(bup(crocidura_ab_nb_re))
sum(croc_ab, na.rm = TRUE)

# Mus musculus abundance --------------------------------------------------

mus_mus_ab <- rodents_ab %>%
  ungroup() %>%
  filter(str_detect(initial_species_id, "musculus")) %>%
  select(site_number, visit, n) %>%
  right_join(long_sites) %>%
  mutate(visit = paste0("x", visit),
         n = replace_na(n, 0)) %>%
  arrange(visit) %>%
  pivot_wider(names_from = visit, values_from = n) %>%
  arrange(site_number) %>%
  select(-site_number) %>%
  as.matrix()

mus_mus_umf <- unmarkedFramePCount(y = mus_mus_ab, siteCovs = site_covs, obsCovs = obs_covs)
summary(mus_mus_umf)

mus_mus_ab_p <- pcount(formula = ~tn + season ~ location + habitat, mus_mus_umf, mixture = "P", K = 105)
mus_mus_ab_nb <- pcount(formula = ~tn + season ~ location + habitat, mus_mus_umf, mixture = "NB", K = 105)
mus_mus_ab_zip <- pcount(formula = ~tn + season ~ location + habitat, mus_mus_umf, mixture = "ZIP", K = 105)

mus_mus_ab_p_re <- ranef(mus_mus_ab_p)
plot(mus_mus_ab_p_re, subset=site %in% 1:22)

sum(bup(mus_mus_ab_p_re))
sum(mus_mus_ab, na.rm = TRUE)

# Mastomys natalensis abundance -------------------------------------------

mas_nat_ab <- rodents_ab %>%
  ungroup() %>%
  filter(str_detect(initial_species_id, "mastomys")) %>%
  select(site_number, visit, n) %>%
  right_join(long_sites) %>%
  mutate(visit = paste0("x", visit),
         n = replace_na(n, 0)) %>%
  arrange(visit) %>%
  pivot_wider(names_from = visit, values_from = n) %>%
  arrange(site_number) %>%
  select(-site_number) %>%
  as.matrix()

mas_nat_umf <- unmarkedFramePCount(y = mas_nat_ab, siteCovs = site_covs, obsCovs = obs_covs)
summary(mas_nat_umf)

mas_nat_ab_p <- pcount(formula = ~tn + season ~ location + habitat, mas_nat_umf, mixture = "P", K = 105)
mas_nat_ab_nb <- pcount(formula = ~tn + season ~ location + habitat, mas_nat_umf, mixture = "NB", K = 105)
mas_nat_ab_zip <- pcount(formula = ~tn + season ~ location + habitat, mas_nat_umf, mixture = "ZIP", K = 105)

mas_nat_ab_p_re <- ranef(mas_nat_ab_nb)
plot(mas_nat_ab_p_re, subset=site %in% 1:22)

sum(bup(mas_nat_ab_p_re))
sum(mas_nat_ab, na.rm = TRUE)
