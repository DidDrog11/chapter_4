source(here::here("R", "00_setup.R"))

rodent_models_summary <- read_rds(here("temp", "rodent_models_summary_2023-02-23.rds"))

# Odds ratios for null model ----------------------------------------------
# What is the probability of a tie forming between nodes
# i.e. what is the probability of contact between two individual individual rodents in a given landuse setting
null_model <- rodent_models_summary$null_model

null_model_plot <- bind_rows(as.data.frame(null_model[[1]]$coefficients),
          as.data.frame(null_model[[2]]$coefficients),
          as.data.frame(null_model[[3]]$coefficients)) %>%
  mutate(Landuse = c("Agriculture", "Forest", "Village")) %>%
  tibble() %>%
  select(Landuse, OR, Lower, Upper, `Pr(>|z|)`) %>%
  mutate(Landuse = fct_rev(factor(Landuse, levels = c("Forest", "Agriculture", "Village")))) %>%
  ggplot() + 
  geom_point(aes(x = OR, y = Landuse)) + 
  geom_segment(aes(x = Lower, xend = Upper, y = Landuse, yend = Landuse)) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  labs(title = "Odds of a tie being observed between two individuals")

main_effects <- rodent_models_summary$main_effects

main_effect_plot <- bind_rows(as.data.frame(main_effects[[1]]$coefficients),
          as.data.frame(main_effects[[2]]$coefficients),
          as.data.frame(main_effects[[3]]$coefficients)) %>%
  mutate(Landuse = rep(c("Agriculture", "Forest", "Village"), each = nrow(main_effects[[1]]$coefficients)),
         Coefficient = fct_rev(fct_inorder(rep(c("Edges", "Species - Crocidura", "Species - Mastomys", "Species - Praomys",
                             "Species - Lophuromys", "Species - Mus musculus", "Species - Rattus rattus", "Species - Mus minutoides"),
                           times = length(main_effects)))),
         Unstable_estimate = fct_rev(factor(case_when(is.finite(Estimate) ~ FALSE,
                                       TRUE ~ TRUE))),
         Significant = case_when(`Pr(>|z|)` <= 0.05 ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  tibble() %>%
  select(Coefficient, Landuse, OR, Lower, Upper, `Pr(>|z|)`, Unstable_estimate, Significant) %>%
  mutate(Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village"))) %>%
  ggplot() + 
  geom_point(aes(x = OR, y = Coefficient, alpha = Unstable_estimate, colour = Significant)) + 
  geom_segment(aes(x = Lower, xend = Upper, y = Coefficient, yend = Coefficient, alpha = Unstable_estimate, colour = Significant)) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  scale_alpha_discrete(range = c(0.1, 1)) +
  facet_wrap(~ Landuse, ncol = 1) +
  coord_cartesian(xlim = c(0, 8)) +
  theme_bw() +
  guides(alpha = "none") +
  labs(title = "Odds of a tie being observed between two individuals by species")

homophily_term <- rodent_models_summary$homophily

coeff_names <- rownames(as.data.frame(homophily_term[[1]]$coefficients)) %>%
  str_replace(., "nodefactor.Species.", "Species - ") %>%
  str_replace(., "nodematch.Species.", "Match - ") %>%
  str_replace(., "edges", "Edges")

homophily_model_plot <- bind_rows(as.data.frame(homophily_term[[1]]$coefficients),
          as.data.frame(homophily_term[[2]]$coefficients),
          as.data.frame(homophily_term[[3]]$coefficients)) %>%
  mutate(Landuse = rep(c("Agriculture", "Forest", "Village"), each = nrow(homophily_term[[1]]$coefficients)),
         Coefficient = fct_rev(fct_inorder(rep(c(coeff_names),
                                               times = length(homophily_term)))),
         Unstable_estimate = fct_rev(factor(case_when(is.finite(Estimate) & is.finite(Upper) & OR != 0 ~ FALSE,
                                                      TRUE ~ TRUE))),
         Significant = case_when(`Pr(>|z|)` <= 0.05 ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  tibble() %>%
  select(Coefficient, Landuse, OR, Lower, Upper, `Pr(>|z|)`, Unstable_estimate, Significant) %>%
  mutate(Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village"))) %>%
  filter(Unstable_estimate == FALSE) %>%
  mutate(Upper = case_when(Upper >= 500 ~ 500,
                           TRUE ~ Upper)) %>%
  ggplot() + 
  geom_point(aes(x = OR, y = Coefficient, colour = Significant)) + 
  geom_segment(aes(x = Lower, xend = Upper, y = Coefficient, yend = Coefficient, colour = Significant)) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  facet_wrap(~ Landuse, ncol = 1, scales = "free_x") +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  guides(alpha = "none") +
  labs(caption = "n.b. Log10 x axis, varies by facet",
       title = "Odds of a tie being observed between two individuals by species and by homophily",
       subtitle = "Unstable estimates removed")

save_plot(plot = main_effect_plot, filename = here("output", "main_effect.png"), base_height = 6, base_width = 8)
save_plot(plot = homophily_model_plot, filename = here("output", "homophily.png"), base_height = 8, base_width = 8)
