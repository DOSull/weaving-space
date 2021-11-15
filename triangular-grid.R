library(dplyr)
library(ggplot2)
library(stringr)
library(sf)
library(tmap)

get_coord <- function(a, b, c, dstep, c0 = 0) {
  return(c0 + a * dstep[1] + b * dstep[2] + c * dstep[3])
}

# See: Nagy, B. N. 2003. Shortest Paths in Triangular Grids with Neighbourhood 
# Sequences. Journal of Computing and Information Technology 11 (2):111.
get_triangle_grid <- function(L = 10000, tilt = 1, n = 5) {
  # Unit vectors in each axial direction
  angles <- (tilt + c(0, 8, 4)) * pi / 6
  dx <- L * cos(angles)
  dy <- L * sin(angles)
  abc <- expand.grid(a = -(2*n):(2*n + 1), b = -n:n, c = -n:n) %>% 
    filter((a + b + c) %in% 0:1) %>%
    mutate(a = a + n + 1, b = b + n + 1, c = c + n + 1) %>%
    mutate(xc = get_coord(a, b, c, dx), 
           yc = get_coord(a, b, c, dy), 
           parity = 1 - ((n + a + b + c) %% 2),
           label = str_c(a, b, c, sep = ","),
           ID = row_number()) %>%
    select(ID, label, a, b, c, xc, yc, parity)
  return(abc)
}

plot_triangle_grid <- function(g) {
  p <- ggplot(g) + 
    geom_point(aes(x = xc, y = yc, 
                   shape = as.factor(parity), 
                   colour = as.factor(parity)), size = 5) +
    scale_shape_manual(values = c(2, 6)) +
    geom_text(aes(x = xc, y = yc, label = label), size = 2, vjust = -1) +
    coord_equal() +
    theme_minimal()
  return(p)
}


get_triangle_from_centre <- Vectorize(
  function(xc, yc, up = TRUE, L = 10000) {
    angles <- (seq(1, 12, 4) + ifelse(up, 0, 2)) * pi / 6
    corners <- (L * c(cos(angles), sin(angles))) + c(rep(xc, 3), rep(yc, 3))
    corners <- c(corners[1:3], corners[1], corners[4:6], corners[4])
    pts <- matrix(corners, 4, 2)
    return(list(pts))
  }
)


tg <- get_triangle_grid(n = 1) 
plot_triangle_grid(tg)
triangles <- get_triangle_from_centre(tg$xc, tg$yc, tg$parity == 1) %>%
  lapply(st_multipoint) %>%
  lapply(st_cast, to = "POLYGON") %>%
  st_sfc()
sf_triangles <- st_sf(tg, geometry = triangles, crs = 3857)
sf_triangles %>% plot(lwd = 0.2)

# A = array(0, dim = rep(11, 3))
# A[cbind(abc$a, abc$b, abc$c)] <- 1



# It's a slice through a 3D matrix
# library(plotly)
# plot_ly(data = abc, x = ~a, y = ~b, z = ~c, size = I(30), marker = list(colorscale = "RdBu"),
#         color = ~parity, type = "scatter3d", mode = "markers")



