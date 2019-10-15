library(rnaturalearth)
library(sf)
library(maps)
library(rmapshaper)
library(tidyverse)
library(ggplot2)
library(sfdct)

Oz_cities <- world.cities %>%
  dplyr::filter(country.etc == "Australia",
                pop > 1000000 | name == "Darwin" | name == "Canberra")

Oz_coast <- rnaturalearth::ne_states("Australia", returnclass = "sf") %>%
  dplyr::filter(code_hasc %in% c("AU.NT",
                                 "AU.WA",
                                 "AU.CT",
                                 "AU.NS",
                                 "AU.SA",
                                 "AU.VI",
                                 "AU.QL")) %>%
  rmapshaper::ms_filter_islands(1e+12) %>% ## get only mainland Australia
  rmapshaper::ms_simplify(0.01) %>% ## simplify coast a little for easier computation later
  sf::st_union() %>% 
  sf::st_transform(3112) ## Geoscience Australia Lambert

ggplot(Oz_coast) +
  geom_sf(fill = "grey20") +
  theme_minimal()

Oz_tri <- sfdct::ct_triangulate(Oz_coast, a = 100000000000)

ggplot(Oz_tri) +
  geom_sf(fill = "grey80") +
  theme_minimal()

Oz_tri_lines <- Oz_tri %>%
  sf::st_cast() %>%
  sf::st_sf() %>%
  sf::st_cast("LINESTRING", do_split = TRUE) %>%
  sf::st_coordinates()
