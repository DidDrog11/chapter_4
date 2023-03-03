# Load packages and project wide values
source(here::here("R", "00_setup.R"))

# Load cleaned data
source(here("R", "01_load_data.R"))

# Further clean the data for this analysis
source(here("R", "02_clean_data.R"))

# Descriptives
source(here("R", "03_descriptive.R"))

# Rodent networks
source(here("R", "05_site_networks.R"))

# Keep all 21/22 networks perform meta-analysis on the network by landuse type
# Limitations around the density of forest settings
# Further sampling needed
# Other rodents set as <50