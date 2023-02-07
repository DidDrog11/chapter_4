source(here::here("R", "00_setup.R"))

rodent_network <- read_rds(here("temp", "rodent_networks_2023-01-04.rds"))
rodent_models <- read_rds(here("temp", "rodent_models_2023-01-04.rds"))
rodent_models_summary <- read_rds(here("temp", "rodent_models_summary_2023-01-04.rds"))

# Odds ratios for null model ----------------------------------------------
# What is the probability of a tie forming between nodes
# i.e. what is the probability of contact between two individual individual rodents in a given landuse setting



