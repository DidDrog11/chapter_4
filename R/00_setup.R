if (!require("pacman")) install.packages("pacman")

pkgs =
  c("cowplot",
    "ergm",
    "flextable",
    "ggdist",
    "ggforce",
    "ggraph",
    "here",
    "HomeRange",
    "igraph",
    "intergraph",
    "lubridate",
    "mapview",
    "metafor",
    "network",
    "osmdata",
    "patchwork",
    "RColorBrewer",
    "rosm",
    "sf",
    "statnet",
    "tidyfast",
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

species_palette <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")
names(species_palette) <- c("Mastomys natalensis", "Rattus rattus", "Mus musculus", "Praomys spp", "Crocidura spp", "Other spp")

# Allocate seasons to months
season <- tibble(month = 1:12, season = c(rep("Dry", 4), rep("Rainy", 6), rep("Dry", 2)))

visit_season <- tibble(visit = c(1:9),
                       season = c(rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2), rep("Rainy", 2), rep("Dry", 1)))

species_order_all <- c("Mastomys natalensis", "Mus musculus", "Rattus rattus", "Lemniscomys striatus",
                       "Lophuromys sikapusi", "Malacomys edwardsi", "Praomys rostratus", "Mus setulosus",
                       "Hylomyscus simus", "Hybomys planifrons", "Gerbilliscus guineae", "Dasymys spp",
                       "Crocidura spp")
