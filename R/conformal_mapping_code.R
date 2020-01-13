library(rnaturalearth)
library(sf)
library(maps)
library(rmapshaper)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggforce)

Oz_cities <- world.cities %>%
  dplyr::filter(country.etc == "Australia",
                pop > 1000000 | name == "Darwin" | name == "Canberra") %>%
  sf::st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  sf::st_transform(3112) ## Geoscience Australia Lambert

Oz_coast <- rnaturalearth::ne_states("Australia", returnclass = "sf") %>%
  dplyr::filter(code_hasc %in% c("AU.NT",
                                 "AU.WA",
                                 "AU.CT",
                                 "AU.NS",
                                 "AU.SA",
                                 "AU.VI",
                                 "AU.QL")) %>%
  rmapshaper::ms_filter_islands(1e+12) %>% ## get only mainland Australia
  rmapshaper::ms_simplify() %>% ## simplify coast for easier computation later
  sf::st_union() %>% 
  sf::st_transform(3112) ## Geoscience Australia Lambert

oz <- ggplot(Oz_coast ) +
  geom_sf(fill = "grey20") +
  geom_sf(data = Oz_cities, colour = "red", size = 3) +
  theme_minimal()

oz +
  ggrepel::geom_label_repel(
    data = Oz_cities,
    aes(label = name, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0,
    segment.color = "red",
    force = 3,
    direction = "both",
    max.iter = 50000,
    nudge_x = 10000,
    nudge_y = -10000
  ) 
  #geom_sf_label(data = Oz_cities, aes(label = name)) +

## We start our circle packing by generating a hexagonal grid on our polygon

hexes <- Oz_coast %>%
  sf::st_make_grid(cellsize = 100000, square = FALSE)

oz +
  geom_sf(data = hexes, colour = "white", alpha = 0.6, fill = NA)

## Now we generate circles at each hexagon vertice and centre
hex_vert <- Oz_coast %>%
  sf::st_make_grid(cellsize = 100000, square = FALSE, what = "polygons") %>%
  sf::st_cast("POINT")
hex_cents <- Oz_coast %>%
  sf::st_make_grid(cellsize = 100000, square = FALSE, what = "centers")

circ_cents <- c(hex_vert, hex_cents)
circ_coords <- circ_cents %>%
  sf::st_coordinates() %>%
  dplyr::as_tibble() %>%
  dplyr::distinct()

rad <- (hexes[[1]] %>% sf::st_cast("LINESTRING") %>% sf::st_length()) / 12

oz +
  ggforce::geom_circle(data = circ_coords, aes(x0 = X, y0 = Y, r = rad), colour = "white")

## Now remove circles with centres outside coast
circle_centres <- circ_cents %>%
  sf::st_as_sf() %>%
  dplyr::distinct() %>%
  sf::st_intersection(Oz_coast)

ggplot(circle_centres %>% sf::st_coordinates() %>% dplyr::as_tibble(),
       aes(x0 = X, y0 = Y, r = rad)) +
  geom_circle() +
  coord_equal() +
  theme_minimal()