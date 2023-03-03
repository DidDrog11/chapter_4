if (!require("pacman")) install.packages("pacman")

pkgs =
  c("cowplot",
    "ergm",
    "flextable",
    "ggforce",
    "ggraph",
    "here",
    "igraph",
    "intergraph",
    "lubridate",
    "mapview",
    "metafor",
    "osmdata",
    "patchwork",
    "RColorBrewer",
    "rosm",
    "sf",
    "statnet",
    "tidygraph",
    "tidyverse",
    "unmarked"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"

village_palette <- c("#7a0177", "#fec44f", "#005120", "#253494")
names(village_palette) <-  c("Lalehun", "Seilama", "Lambayama", "Baiama")

landuse_palette <- c("#00913a", "#FEC44F", "#a13b9e")
names(landuse_palette) <- c("Forest", "Agriculture", "Village")

# Allocate seasons to months
season <- tibble(month = 1:12, season = c(rep("Dry", 4), rep("Rainy", 6), rep("Dry", 2)))

visit_season <- tibble(visit = c(1:8),
                       season = c(rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2), rep("Rainy", 2)))
