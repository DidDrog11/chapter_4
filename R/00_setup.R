if (!require("pacman")) install.packages("pacman")

pkgs =
  c("broom",
    "brms",
    "cowplot",
    "ergm",
    "flextable",
    "ggdist",
    "ggforce",
    "ggraph",
    "ggmap",
    "ggspatial",
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
    "rmapshaper",
    "rnaturalearth",
    "rosm",
    "sf",
    "statnet",
    "terra",
    "tidyfast",
    "tidygraph",
    "tidyterra",
    "tidyverse",
    "unmarked"
  )

pacman::p_load(pkgs, character.only = T)

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"

village_order <- c("baiama","lalehun", "lambayama", "seilama")

village_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#888888")
names(village_palette) <-  c("Lalehun", "Seilama", "Lambayama", "Baiama")

landuse_palette <- c("#00913A", "#FEC44F", "#A13B9E")
names(landuse_palette) <- c("Forest", "Agriculture", "Village")

all_species_order <- c("Mastomys natalensis", "Rattus rattus", "Mus musculus", "Crocidura olivieri", "Praomys rostratus", "Lophuromys sikapusi",
                       "Mus setulosus", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                       "Hylomyscus simus", "Hybomys planifrons", "Mastomys erythroleucus", "Crocidura theresae", "Gerbilliscus guineae", "Dasymys rufulus")

species_order_plots <- c("Mastomys natalensis", "Rattus rattus", "Mus musculus", "Crocidura olivieri", "Praomys rostratus", "Lophuromys sikapusi",
                         "Mus setulosus", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                         "Hylomyscus simus", "Hybomys planifrons", "Mastomys erythroleucus", "Crocidura theresae", "Gerbilliscus guineae", "Dasymys rufulus")

species_order_n <- c("Mastomys natalensis", "Crocidura olivieri", "Praomys rostratus", "Mus musculus", "Rattus rattus", "Lophuromys sikapusi",
                     "Mus setulosus", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                     "Hylomyscus simus", "Hybomys planifrons", "Crocidura theresae", "Mastomys erythroleucus", "Gerbilliscus guineae", "Dasymys rufulus")

group_landuse_palette <- c("#00913a", "#FEC44F", "#F7A820", "#A13B9E", "#5407A6")
names(group_landuse_palette) <- c("Forest - Rural", "Agriculture - Rural", "Agriculture - Peri-urban", "Village - Rural", "Village - Peri-urban")

species_palette <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")
names(species_palette) <- c("Mastomys natalensis", "Rattus rattus", "Mus musculus", "Praomys rostratus", "Crocidura olivieri", "Other spp")

# Allocate seasons to months
season <- tibble(month = 1:12, season = c(rep("Dry", 4), rep("Rainy", 6), rep("Dry", 2)))

visit_season <- tibble(visit = 1:10, season = c(rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2)))
