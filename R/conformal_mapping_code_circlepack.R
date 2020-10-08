library(rnaturalearth)
library(sf)
library(maps)
library(rmapshaper)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(distances)
library(spdep)
library(packcircles)

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
  rmapshaper::ms_filter_islands(1e+14) %>% ## get only mainland Australia
  rmapshaper::ms_simplify(0.05) %>% ## aggressive simplification to reduce number if 'insets'
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

cellsize <- 100000

hexes <- Oz_coast %>%
  sf::st_make_grid(cellsize = cellsize, square = FALSE)

oz +
  geom_sf(data = hexes, colour = "white", alpha = 0.6, fill = NA)

## Now we generate circles at each hexagon vertice and centre
hex_vert <- Oz_coast %>%
  sf::st_make_grid(cellsize = cellsize, square = FALSE, what = "polygons") %>%
  sf::st_cast("POINT") %>%
  lwgeom::st_snap_to_grid(25)
hex_cents <- Oz_coast %>%
  sf::st_make_grid(cellsize = cellsize, square = FALSE, what = "centers") %>%
  sf::st_centroid() %>%
  lwgeom::st_snap_to_grid(25)

circ_cents <- c(hex_vert, hex_cents)
circ_coords <- circ_cents %>%
  sf::st_coordinates() %>%
  dplyr::as_tibble() %>%
  dplyr::distinct()

rad <- (hexes[[1]] %>% sf::st_cast("LINESTRING") %>% sf::st_length()) / 12

oz +
  ggforce::geom_circle(data = circ_coords, 
                       aes(x0 = X, y0 = Y, r = rad), colour = "white")

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

## Now make a network of tangent circles, discard any not tangent to at least 3 others
circle_dists <- circle_centres %>%
  sf::st_coordinates() %>%
  distances::distances()

closest <- circle_dists %>%
  distances::nearest_neighbor_search(8)

# index <- closest[ , 1]
find_tangent_indices <- function(index, rad) {
  dists <- distances::distance_columns(circle_dists, index, closest[ , index]) 
  rownames(dists)[dists < (2.1 * rad)] %>% as.integer()
}

tangents <- purrr::map(seq_len(ncol(closest)), ~find_tangent_indices(.x, rad))
class(tangents) <- "nb"
tangent_mat <- spdep::nb2mat(tangents, style = "B")

nnum <- sapply(tangents, length) - 1
sum(nnum > 6)
sum(nnum < 3)

## remove circles with less than 3 tangent circles
num_less <- sum(nnum < 3)
while(num_less != 0) {
  circle_centres <- circle_centres[!nnum < 3, ]
  circle_dists <- circle_centres %>%
    sf::st_coordinates() %>%
    distances::distances()
  
  closest <- circle_dists %>%
    distances::nearest_neighbor_search(8)
  
  tangents <- purrr::map(seq_len(ncol(closest)), ~find_tangent_indices(.x, rad))
  nnum <- sapply(tangents, length) - 1
  
  num_less <- sum(nnum < 3)
}

sum(nnum < 3)
sum(nnum > 6)

ggplot(circle_centres %>% sf::st_coordinates() %>% dplyr::as_tibble(),
       aes(x0 = X, y0 = Y, r = rad)) +
  geom_circle() +
  coord_equal() +
  theme_minimal()


nnet <- lapply(seq_along(tangents), function(x) tangents[[x]][tangents[[x]] != x])

pnts <- circle_centres[nnet[[1]], ] %>%
  sf::st_coordinates()
order_points <- function(pnts, clockwise = FALSE) {
  cent <- apply(pnts, 2, mean)
  pnts2 <- t(t(pnts) - cent)
  angles <- atan2(pnts2[ , 2], pnts2[ , 1])
  if(clockwise) {
    indices <- order(angles, decreasing = TRUE)
  } else {
    indices <- order(angles)
  }
}

coords <- circle_centres %>%
  sf::st_coordinates()
all_pnts <- lapply(nnet, function(x) coords[x, ])
nnet <- lapply(seq_along(all_pnts), function(x) nnet[[x]][order_points(all_pnts[[x]])])
nnet <- lapply(seq_along(nnet), function(x) c(x, nnet[[x]]))




internal <- nnet[sapply(nnet, length) == 7]
external <- dplyr::tibble(id = which(sapply(nnet, length) < 7),
                          radius = 10000)

test <- circleGraphLayout(internal, external)

ggplot(test[test$radius != 10000, ],
       aes(x0 = x, y0 = y, r = radius)) + geom_circle() +
  theme_minimal()

ggplot(test,
       aes(x0 = x, y0 = y, r = radius)) + geom_circle() +
  theme_minimal()

