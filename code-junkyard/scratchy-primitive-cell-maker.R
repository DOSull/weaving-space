source("../tiling-space-utils.R")

library(magrittr)
library(dplyr)

s60 <- sin(pi / 3)

unit_square <- get_polygon(c(0, 0, 1, 0, 1, 1, 0, 1))

bp <- unit_square * t(matrix(c(2, 0, 0.5, s60), 2, 2)) - c(0.5, s60)
plot(bp)

base_trs <- sapply(seq(0, 240, 120), wk_affine_rotate)
base_trs

bps <- lapply(base_trs, FUN = wk_transform, handleable = bp) %>%
  sapply(st_as_sfc) %>%
  st_sfc()

plot(bps)

unit_trs <- c(
  0, 0,
  1.5, -3 * s60,
  3, 0, 
  1.5, 3 * s60
)

polys <- list()
item <- 1
for (p in bps) {
  for (i in seq(1, length(unit_trs), 2)) {
    polys[[item]] <- p + unit_trs[c(i, i + 1)]
    item <- item + 1
  }
}

polys <- polys %>% 
  st_as_sfc() %T>% 
  plot()

tile <- get_polygon(unit_trs)

funit <- polys %>% 
  st_intersection(tile) %>%
  st_sf() %>%
  dplyr::filter(st_geometry_type(.) == "POLYGON", st_area(.) > 1e-10) %T>%
  plot()
