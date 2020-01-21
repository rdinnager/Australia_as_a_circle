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
library(nimble)
library(nCompiler)

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
  rmapshaper::ms_simplify(0.01) %>% ## aggressive simplification to reduce number if 'insets'
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
  sf::st_make_grid(cellsize = 50000, square = FALSE)

oz +
  geom_sf(data = hexes, colour = "white", alpha = 0.6, fill = NA)

## Now we generate circles at each hexagon vertice and centre
hex_vert <- Oz_coast %>%
  sf::st_make_grid(cellsize = 50000, square = FALSE, what = "polygons") %>%
  sf::st_cast("POINT") %>%
  lwgeom::st_snap_to_grid(25)
hex_cents <- Oz_coast %>%
  sf::st_make_grid(cellsize = 50000, square = FALSE, what = "centers") %>%
  lwgeom::st_snap_to_grid(25)

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

######## start circle packing algorithm ##########
## now we transform our circle's radii and centres so that they pack a unit disk, instead of Australia
## to make this algorithm faster we will use the nCompiler package to compile our code to C

## first let's order our tangent circles in clockwise order
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

make_triangles <- function(index, tangent_indices) {
  triangles <- matrix(tangent_indices[1:2], nrow = 1)
  for(j in 2:(length(tangent_indices) - 1)) {
    triangles <- rbind(triangles, tangent_indices[j:(j + 1)])
  }
  triangles <- rbind(triangles, tangent_indices[c(length(tangent_indices), 1)])
  triangles <- cbind(index, triangles)
  #triangles <- t(apply(triangles, 1, function(x) x[order_points(coords[x, ])]))
  triangles
}

r <- current_radii[1, ]
cosine_law <- function(r, inf_tol = 1000) {
  infs <- r >= inf_tol
  if(infs[1]) {
    angle <- 0
  } else {
    if(infs[2] & infs[3]) {
      z <- 1 - 2 * exp(-2 * r[1])
      if(z < -1) {
        z <- -1
      }
      if(z > 1) {
        z <- 1
      }
      angle <- acos(z)
    } else {
      if(infs[2] & !infs[3]) {
        z <- (cosh(r[1] + r[3]) - exp(r[3] - r[1])) /
          sinh(r[1] + r[3])
        if(z < -1) {
          z <- -1
        }
        if(z > 1) {
          z <- 1
        }
        angle <- acos(z)
      }
      if(!infs[2] & infs[3]) {
        z <- (cosh(r[1] + r[2]) - exp(r[2] - r[1])) /
          sinh(r[1] + r[2])
        if(z < -1) {
          z <- -1
        }
        if(z > 1) {
          z <- 1
        }
        angle <- acos(z)
      }
      if(!infs[2] & !infs[3]) {
        z <- (cosh(r[1] + r[2]) * cosh(r[1] + r[3]) - cosh(r[2] + r[3])) /
          (sinh(r[1] + r[2]) * sinh(r[1] + r[3]))
        if(z < -1) {
          z <- -1
        }
        if(z > 1) {
          z <- 1
        }
        angle <- acos(z)
      }
    }
  }
  angle
}

r <- current_radii[1, ]
cosine_law_eff <- function(r) {
  angle <- 2*asin(sqrt((r[2] / (r[1] + r[2])) * (r[3] / (r[1] + r[3]))))
  angle
}

cosine_law_n <- nCompiler::nFunction(
  fun = cosine_law_eff,
  argTypes = list(r = 'numericVector()'),
  returnType = 'numericScalar()'
)

cosine_law_c <- nCompiler::nCompile_nFunction(cosine_law_n)

cosine_law_n <- nimbleFunction(run = function(r = double(1)){ ## type declarations
  angle <- 2*asin(sqrt((r[2] / (r[1] + r[2])) * (r[3] / (r[1] + r[3]))))
  return(angle)
  returnType(double(0)) ## return type declaration
})

cosine_law_c <- compileNimble(cosine_law_n)

cosine_law_c(c(0.1, 0.1, 0.1))
cosine_law_eff(c(0.1, 0.1, 0.1))

boundary <- nnum < 6
triangles <- lapply(seq_along(nnet), function(x) make_triangles(x, nnet[[x]]))
circle_pack <- function(triangles, boundary, tol = 0.0001, smallest_radius = 1e-10) {
  radii <- runif(length(triangles), 0.4, 0.68)
  radii[boundary] <- 1000
  error <- rep(99, sum(!boundary))
  
  internal <- which(!boundary)
  while(any(error > tol)) {
    for(i in seq_along(internal)) {
      index <- internal[i]
      current_circle <- triangles[[index]]  
      current_radii <- matrix(radii[current_circle], ncol = 3)
      
      angle_sum <- sum(apply(current_radii, 1, cosine_law_c))
      
      if(!is.finite(angle_sum)) {
        print(index)
        stop("angle sum not finite")
      }
      
      e <- 2 * pi - angle_sum
      if(abs(e) > tol) {
        
        beta <- sin((angle_sum) / (2 * 6))
        delta = sin((2 * pi) / (2 * 6))
        
        v <- (beta / (1 - beta)) * radii[index]
        u <- ((1 - delta) / delta) * v
        
        # d <- min(abs(e), 1)
        # if(e > 0) {
        #   radii[index] <- (1 - (1 / 10) * d) * radii[index]
        # } else {
        #   radii[index] <- (1 + (1 / 10) * d) * radii[index]
        # }
        radii[index] <- u
      }
      error[i] <- abs(e)
      #radii[index]
    }
    print(sum(error))
  }
}

###### Try out Nimble solution ###########

triangle_array <- do.call(function(...) abind(..., along = 0), triangles[!boundary])

circle_centre <- nimbleFunction(
  setup = function(triangle_array, boundary, tol = 0.001) {
    triangle_array <- triangle_array
    boundary <- boundary
    nonboundary <- which(!boundary)
    radii <- runif(length(boundary), 0.4, 0.68)
    error <- rep(99, length(boundary))
  },
  run = function() {
    
    while(any(error > tol)) {
     for(index in 1:dim(triangle_array)[1]) {
       current_circle <- triangle_array[index, , ]  
       current_radii <- matrix(radii[c(current_circle)], ncol = 3)
       
       angles <- numeric(dim(current_radii)[1])
       for(j in 1:(dim(current_radii)[1])) {
         angles[j] <- 2*asin(sqrt((current_radii[j, 2] / (current_radii[j, 1] + current_radii[j, 2])) * (current_radii[j, 3] / (current_radii[j, 1] + current_radii[j, 3]))))
       }
       angle_sum <- sum(angles)
       
       e <- 2 * pi - angle_sum
       if(abs(e) > tol) {
         
         beta <- sin((angle_sum) / (2 * 6))
         delta = sin((2 * pi) / (2 * 6))
         
         v <- (beta / (1 - beta)) * radii[index]
         u <- ((1 - delta) / delta) * v
         
         # d <- min(abs(e), 1)
         # if(e > 0) {
         #   radii[index] <- (1 - (1 / 10) * d) * radii[index]
         # } else {
         #   radii[index] <- (1 + (1 / 10) * d) * radii[index]
         # }
         radii[nonboundary[index]] <<- u
       }
       error[index] <<- abs(e)

     }
    print(sum(error))
       
  }
   
  }
)

circle_centre_setup <- circle_centre(triangle_array, boundary)
circle_centre_c <-  compileNimble(circle_centre_setup, showCompilerOutput = TRUE)

circle_centre_c$error

circle_centre_c$run()