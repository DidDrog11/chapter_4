source(here::here("R", "00_setup.R"))

rodent_models_summary <- read_rds(here("temp", "rodent_models_summary_2023-03-02.rds"))
rodent_network <- read_rds(here("temp", "rodent_networks_2023-03-02.rds"))
rodent_data <- read_rds(here("data", "expanded_assemblages_2023-03-02.rds"))

`%s%` <- network::`%s%`

# Odds ratios for model ----------------------------------------------
# What is the probability of a tie forming between nodes
# i.e. what is the probability of contact between two individual rodents in a given landuse setting

homophily_term <- rodent_models_summary$homophily

coefficients <- list()

for(i in 1:length(homophily_term)) {
  
  if(is.list(homophily_term[[i]])) {
    
  coeff = homophily_term[[i]][["coefficients"]] %>%
    as.data.frame()
  
  coeff$coefficient = row.names(coeff) %>%
    str_replace(., "nodefactor.Species.", "Species - ") %>%
    str_replace(., "nodematch.Species.", "Match - ") %>%
    str_replace(., "edges", "Edges")
  
  coeff$Network = i
  coeff$Visit = unique(rodent_data$nodelist[[i]]$Visit)
  coeff$Landuse = unique(rodent_data$nodelist[[i]]$Landuse)
  coeff$`Observed M. natalensis` = rodent_data$nodelist[[i]] %>%
    filter(Species == "Mastomys natalensis" & Observed == TRUE) %>%
    nrow()
  coeff$`Observed all species` =  rodent_data$nodelist[[i]] %>%
    filter(Observed == TRUE) %>%
    nrow()
  coeff$`Unobserved M. natalensis` = rodent_data$nodelist[[i]] %>%
    filter(Species == "Mastomys natalensis" & Observed == FALSE) %>%
    nrow()
  coeff$N = nrow(rodent_data$nodelist[[i]])
  
  coefficients[[i]] <- coeff %>%
    tibble() %>%
    select(Network, Visit, Landuse, N, `Observed all species`, `Observed M. natalensis`, `Unobserved M. natalensis`, Coefficient = coefficient, Estimate, `Std. Error`)
  
  } else {
    
    coefficients[[i]] <- tibble(Network = i,
                                Visit = unique(rodent_data$nodelist[[i]]$Visit),
                                Landuse = unique(rodent_data$nodelist[[i]]$Landuse),
                                `Observed M. natalensis` = rodent_data$nodelist[[i]] %>%
                                  filter(Species == "Mastomys natalensis" & Observed == TRUE) %>%
                                  nrow(),
                                N = nrow(rodent_data$nodelist[[i]]),
                                `Observed all species` =  rodent_data$nodelist[[i]] %>%
                                  filter(Observed == TRUE) %>%
                                  nrow(),
                                `Unobserved M. natalensis` = rodent_data$nodelist[[i]] %>%
                                  filter(Species == "Mastomys natalensis" & Observed == FALSE) %>%
                                  nrow(),
                                Coefficient = NA,
                                Estimate = NA,
                                `Std. Error` = NA) %>%
      tibble() %>%
      select(Network, Visit, Landuse, N, `Observed all species`, `Observed M. natalensis`, `Unobserved M. natalensis`, Coefficient, Estimate, `Std. Error`)
    
    
  }
  
}

coefficients_df <- bind_rows(coefficients) %>%
  mutate(Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village")),
         Estimate = case_when(Estimate == -Inf ~ as.numeric(NA),
                              TRUE ~ Estimate),
         Variance = `Std. Error`^2) %>%
  group_by(Landuse) %>%
  mutate(Network = paste0(Landuse, " ", Visit)) %>%
  filter(`Observed M. natalensis` > 1) %>% 
  filter(!Network %in% Network[is.na(Estimate)])


# Meta-analysis of edges --------------------------------------------------

edges <- coefficients_df %>%
  filter(Coefficient == "Edges") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .,
      slab = Network)

  # for multi-level meta-analysis
  #   rma.mv(yi = Estimate, V = Variance, data = .,
  #   random = ~ 1 | Landuse/Visit,
  #   slab = Network)

edges_rma <- bind_rows(tibble(ES = edges$yi, SE = sqrt(edges$vi), Type = "Network", Network = factor(edges$data$Network), N = edges$data$`Observed M. natalensis`),
                       tibble(ES = edges$b[, 1], SE = edges$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

edges_combined <- ggplot(edges_rma) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed (all species)",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis")


# Meta-analysis of edges stratified by landuse ----------------------------

edges_ag <- coefficients_df %>%
  filter(Coefficient == "Edges" & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

edges_vil <- coefficients_df %>%
  filter(Coefficient == "Edges" & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

edges_rma_df <- bind_rows(tibble(ES = edges_ag$yi, SE = sqrt(edges_ag$vi), Type = "Network",
                                 Weight = weights(edges_ag)/10,
                                 Network = factor(edges_ag$data$Network),
                                 N = edges_ag$data$`Observed M. natalensis`),
                          tibble(ES = edges_ag$b[, 1], SE = edges_ag$se, Type = "Summary",
                                 Weight = 10,
                                 Network = factor("Summary - Agriculture"),
                                 N = as.numeric(NA)),
                          tibble(ES = edges_vil$yi, SE = sqrt(edges_vil$vi), Type = "Network",
                                 Weight = weights(edges_vil)/10,
                                 Network = factor(edges_vil$data$Network),
                                 N = edges_vil$data$`Observed M. natalensis`),
                          tibble(ES = edges_vil$b[, 1], SE = edges_vil$se, Type = "Summary",
                                 Weight = 10,
                                 Network = factor("Summary - Village"),
                                 N = as.numeric(NA)))

edges_stratified <- ggplot(edges_rma_df) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N, size = Weight)) +
  scale_size(range = c(0.1, 1)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed (all species)",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis")


# Meta-analysis of ties based on species type -----------------------------

species <- coefficients_df %>%
  filter(str_detect(Coefficient, "Species")) %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_rma <- bind_rows(tibble(ES = species$yi, SE = sqrt(species$vi), Type = "Network", Network = factor(species$data$Network), N = species$data$`Observed M. natalensis`),
                       tibble(ES = species$b[, 1], SE = species$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

species_combined <- ggplot(species_rma) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed for M. natalensis",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis")

species_ag <- coefficients_df %>%
  filter(str_detect(Coefficient, "Species") & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_vil <- coefficients_df %>%
  filter(str_detect(Coefficient, "Species") & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

species_rma_df <- bind_rows(tibble(ES = species_ag$yi, SE = sqrt(species_ag$vi), Type = "Network",
                                   Weight = weights(species_ag)/10,
                                   Network = factor(species_ag$data$Network),
                                   N = species_ag$data$`Observed M. natalensis`),
                            tibble(ES = species_ag$b[, 1], SE = species_ag$se, Type = "Summary", Weight = 10, Network = factor("Summary - Agriculture"), N = as.numeric(NA)),
                            tibble(ES = species_vil$yi, SE = sqrt(species_vil$vi), Type = "Network",
                                   Weight = weights(species_vil)/10,
                                   Network = factor(species_vil$data$Network), N = species_vil$data$`Observed M. natalensis`),
                            tibble(ES = species_vil$b[, 1], SE = species_vil$se, Type = "Summary", Weight = 10, Network = factor("Summary - Village"), N = as.numeric(NA)))

species_stratified <- ggplot(species_rma_df) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N, size = Weight)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  scale_size(range = c(0.1, 1)) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "N - M. natalensis",
       title = "Odds of a tie being observed for M. natalensis",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis")


# Meta-analysis of intra-specific ties ------------------------------------

match <- coefficients_df %>%
  filter(str_detect(Coefficient, "Match")) %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_rma <- bind_rows(tibble(ES = match$yi, SE = sqrt(match$vi), Type = "Network", 
                              Network = factor(match$data %>%
                                                 drop_na(Estimate) %>%
                                                 pull(Network)), 
                              N = match$data %>%
                                drop_na(Estimate) %>%
                                pull(`Observed M. natalensis`)),
                       tibble(ES = match$b[, 1], SE = match$se, Type = "Summary", Network = factor("Summary"), N = as.numeric(NA)))

match_combined <- ggplot(match_rma) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = "Network number",
       colour = "N M. natalensis",
       title = "Odds of a tie between two M. natalensis being observed",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis") +
  coord_cartesian(xlim = c(0, 40))

match_ag <- coefficients_df %>%
  filter(str_detect(Coefficient, "Match") & Landuse == "Agriculture") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_vil <- coefficients_df %>%
  filter(str_detect(Coefficient, "Match") & Landuse == "Village") %>%
  rma(yi = Estimate, sei = `Std. Error`, data = .)

match_rma_df <- bind_rows(tibble(ES = match_ag$yi, SE = sqrt(match_ag$vi), Type = "Network", 
                                 Weight = weights(match_ag)/10,
                                 Network = factor(match_ag$data %>%
                                                    drop_na(Estimate) %>%
                                                    pull(Network)),
                                 N = match_ag$data %>%
                                   drop_na(Estimate) %>%
                                   pull(`Observed M. natalensis`)),
                          tibble(ES = match_ag$b[, 1], SE = match_ag$se, Type = "Summary", Weight = 10, Network = factor("Summary - Agriculture"), N = as.numeric(NA)),
                          tibble(ES = match_vil$yi, SE = sqrt(match_vil$vi), Type = "Network",
                                 Weight = weights(match_vil)/10,
                                 Network = factor(match_vil$data %>%
                                                    drop_na(Estimate) %>%
                                                    pull(Network)),
                                 N = match_vil$data %>%
                                   drop_na(Estimate) %>%
                                   pull(`Observed M. natalensis`)),
                          tibble(ES = match_vil$b[, 1], SE = match_vil$se, Type = "Summary", Weight = 10, Network = factor("Summary - Village"), N = as.numeric(NA)))

match_stratified <- ggplot(match_rma_df) +
  geom_pointrange(aes(y = fct_rev(Network), x = exp(ES), xmin = exp(ES - (1.96 * SE)), xmax = exp(ES + (1.96 * SE)), shape = Type, colour = N, size = Weight)) +
  geom_vline(aes(xintercept = 1), lty = 2, linewidth = 1) +
  scale_size(range = c(0.1, 1)) +
  theme_bw() +
  labs(x = "Odds Ratio",
       y = element_blank(),
       colour = "N - M. natalensis",
       title = "Odds of a tie between two M. natalensis being observed",
       subtitle = "Networks are limited to those containing more than 1 M. natalensis") +
  coord_cartesian(xlim = c(0, 40))
