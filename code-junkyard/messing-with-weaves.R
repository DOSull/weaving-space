library(pracma)
library(tidyr)
library(sf)
library(tmap)

make_random_pattern <- function(warp = seq(1, 24, 6), weft = seq(27, 50, 6), 
                                repeat_warp = 4, repeat_weft = 4) {
  size <- length(warp) * length(weft)
  mat <- matrix(sample(0:1, size, replace = TRUE), 4, 4)
  treadling <- matrix(sample(0:1, size, replace = TRUE), 
                      nrow = length(weft), ncol = length(warp)) %>%
    pracma::repmat(m = repeat_weft, n = 1)
  threading <- matrix(sample(0:1, size, replace = TRUE), 
                      nrow = length(weft), ncol = length(warp)) %>%
    pracma::repmat(m = 1, n = repeat_warp)
  ww <- threading %*% mat %*% treadling > 0
  fabric <- matrix(0, nrow = nrow(ww), ncol = ncol(ww))
  for (i in seq_along(weft)) {
    rows <- seq(i, nrow(ww), length(weft))
    idxs <- which(ww[rows, ])
    fabric[rows, ][idxs] <- weft[i]
  }
  for (i in seq_along(warp)) {
    cols <- seq(i, ncol(ww), length(warp))
    idxs <- which(!ww[, cols])
    fabric[, cols][idxs] <- warp[i]
  }
  return(fabric)
}

# convenience function to make a sfc POLYGON from a vector
# of points as c(x0,y0,x1,y1,...xn,yn) - note not closed
# this function will close it
get_polygon <- function(pts) {
  mpts <- matrix(c(pts, pts[1:2]), ncol = 2, byrow = TRUE)
  return(st_polygon(list(mpts)))
}

shift_poly <- function(pt, poly) {
  return(poly + pt)
}

make_polygons_from_matrix <- function(w, m, crs) {
  dxdy <- expand_grid(x = 0:(ncol(m) - 1), y = 0:(nrow(m) - 1)) %>% 
    mutate(x = w * x, y = w * y) %>%
    t() %>%
    as.data.frame() %>%
    as.list()
  p <- get_polygon(w / 2 * c(-1, -1, 1, -1, 1, 1, -1, 1))
  ids <- c(m)
  return(
    lapply(dxdy, shift_poly, poly = p) %>% 
      st_as_sfc() %>%
      st_sf() %>%
      st_set_crs(crs) %>%
      mutate(id = ids) %>%
      group_by(id) %>%
      dplyr::summarise()
  )
}

# as vector on a map
make_polygons_from_matrix(10000, 
                          make_random_pattern(repeat_warp = 6, repeat_weft = 6), 
                          3857) %>%
  tm_shape() + 
  tm_polygons(col = "id", style = "cat", palette = "Set1")

# as an array of bitmaps
par(mfrow = c(4, 6)) 
for (i in 1:24) {
  image(make_random_pattern(), asp = 1, col = RColorBrewer::brewer.pal(8, "Spectral"))
}
