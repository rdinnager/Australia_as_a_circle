library(rnaturalearth)
library(sf)
library(maps)
library(rmapshaper)
library(tidyverse)
library(ggplot2)
library(sfdct)
library(tidygraph)
library(ggraph)
library(particles)

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

Oz_tri_edges <- Oz_tri %>%
  sf::st_cast() %>%
  sf::st_sf() %>%
  sf::st_cast("MULTILINESTRING", do_split = FALSE) %>%
  sf::st_geometry() %>%
  purrr::map(~as.matrix(.)) %>%
  purrr::map(~lapply(list(1:2, 2:3, 3:4), function(x) .[x, ])) %>%
  purrr::map(~sf::st_multilinestring(.)) %>%
  sf::st_sfc() %>%
  sf::st_sf() %>%
  dplyr::mutate(triangle_id = 1:n()) %>%
  sf::st_cast("LINESTRING", do_split = TRUE) %>%
  dplyr::mutate(edge_id = 1:n())

Oz_tri_nodes <- Oz_tri_edges %>%
  sf::st_coordinates() %>%
  dplyr::as_tibble() %>%
  dplyr::rename(edge_id = L1) %>%
  dplyr::group_by(edge_id) %>%
  dplyr::slice(c(1, n())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
  dplyr::mutate(node_id = dplyr::group_indices(., X, Y))


source_nodes <- Oz_tri_nodes %>%
  dplyr::filter(start_end == 'start') 

target_nodes <- Oz_tri_nodes %>%
  dplyr::filter(start_end == 'end')

Oz_tri_edges = Oz_tri_edges %>%
  dplyr::left_join(source_nodes %>%
                     dplyr::select(edge_id,
                                   from = node_id)) %>%
  dplyr::left_join(target_nodes %>%
                     dplyr::select(edge_id,
                                   to = node_id))

plot(Oz_tri_edges)
  
Oz_tri_nodes <- Oz_tri_nodes %>%
  dplyr::distinct(node_id, .keep_all = TRUE) %>%
  dplyr::select(-edge_id, -start_end) %>%
  sf::st_as_sf(coords = c('X', 'Y')) %>%
  sf::st_set_crs(st_crs(Oz_tri_edges))

plot(Oz_tri_nodes)

Oz_tri_graph = tidygraph::tbl_graph(nodes = Oz_tri_nodes, edges = as_tibble(Oz_tri_edges), directed = FALSE)

Oz_tri_graph

Oz_tri_graph <- Oz_tri_graph %>%
  tidygraph::activate(edges) %>%
  dplyr::mutate(length = sf::st_length(geometry))

ggraph(Oz_tri_graph) +
  geom_edge_link()

ggraph(Oz_tri_graph, "fr") +
  geom_edge_link()

ggraph(Oz_tri_graph, "kk") +
  geom_edge_link()
