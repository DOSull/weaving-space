library(dplyr)
library(ggplot2)
library(stringr)
library(sf)
library(tmap)

source("biaxial-weave-units.R")

get_coord <- function(a, b, c, dstep, c0 = 0) {
  return(c0 + a * dstep[1] + b * dstep[2] + c * dstep[3])
}

# See: Nagy, B. N. 2003. Shortest Paths in Triangular Grids with Neighbourhood 
# Sequences. Journal of Computing and Information Technology 11 (2):111.
get_triangle_grid <- function(L = 10000, n = 5, diamond = TRUE) {
  # Unit vectors in each axial direction
  angles <- c(3, 7, 11) * pi / 6
  dx <- L * cos(angles)
  dy <- L * sin(angles)
  if (diamond) {
    abc <- expand.grid(a = 1:(2*n), b = 1:n, c = 1:n) 
    # abc <- expand.grid(a = -(2*n):(2*n), b = (1-n):n, c = (1-n):n) 
    score <- (4 * n) %/% 2
  } else {
    abc <- expand.grid(a = 1:n, b=1:n, c = 1:n)
    # abc <- expand.grid(a = -n:(n-1), b = (1-n):n, c = (1-n):n)
    score <- (3 * (n + 1)) %/% 2
  }
  print(score)
  abc <- abc %>% 
    filter((a + b + c) %in% (0:1 + score)) %>%
    mutate(parity = ((a + b + c) %% 2 + (score %% 2)) %% 2,
           # a = a - min(a) + 1, 
           # b = b - min(b) + 1, 
           # c = c - min(c) + 1,
           xc = get_coord(a, b, c, dx), 
           yc = get_coord(a, b, c, dy), 
           label = str_c(a, b, c, sep = ","),
           ID = row_number()) %>%
    select(ID, label, a, b, c, xc, yc, parity)
  return(abc)
}

plot_triangle_grid <- function(g) {
  p <- ggplot(g) + 
    geom_point(aes(x = xc, y = yc, 
                   shape = as.factor(parity), 
                   colour = as.factor(parity)), size = 14) +
    scale_shape_manual(values = c(6, 2)) + # up and down triangles
    geom_text(aes(x = xc, y = yc, label = label), size = 2, vjust = -1) +
    coord_equal() +
    theme_minimal()
  return(p)
}


get_triangle_from_centre <- Vectorize(
  function(xc, yc, up = TRUE, L = 10000) {
    angles <- (seq(1, 12, 4) + ifelse(up, 2, 0)) * pi / 6
    corners <- (L * c(cos(angles), sin(angles))) + c(rep(xc, 3), rep(yc, 3))
    corners <- c(corners[1:3], corners[1], corners[4:6], corners[4])
    pts <- matrix(corners, 4, 2)
    return(list(pts))
  }
)


tg <- get_triangle_grid(L = 1, n = 4, diamond = FALSE) 
plot_triangle_grid(tg)

triangles <- get_triangle_from_centre(tg$xc, tg$yc, tg$parity == 1) %>%
  lapply(st_multipoint) %>%
  lapply(st_cast, to = "POLYGON") %>%
  st_sfc()

sf_triangles <- st_sf(tg, geometry = triangles, crs = 3857)
sf_triangles %>% plot(lwd = 0.2)

range_a <- max(tg$a) - min(tg$a) + 1 
range_b <- max(tg$b) - min(tg$b) + 1 
range_c <- max(tg$c) - min(tg$c) + 1 
A = array(0, dim = c(range_a, range_b, range_c))
A[cbind(tg$a, tg$b, tg$c)] <- 1

# It's a slice through a 3D matrix
# library(plotly)
# plot_ly(data = tg, x = ~a, y = ~b, z = ~c, size = I(30), marker = list(colorscale = "RdBu"),
#         color = ~parity, type = "scatter3d", mode = "markers")
