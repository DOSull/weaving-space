# library(h3jsr)
library(h3)
library(sf)
library(nngeo)
library(tmap)
library(dplyr)
# library(units)
library(s2)

nz <- st_read("~/Documents/geodata/nz-2193.gpkg") %>%
  st_remove_holes() %>%
  st_make_valid() %>%
  # st_buffer(100) %>%
  # st_buffer(-100) %>%
  st_transform(4326)

cells <- nz %>%
  polyfill(res = 7) %>%
  compact() %>%
  h3_to_geo_boundary_sf()

# 
# while (TRUE) {
#   new_cells <- untiled %>%
#     polyfill(res = res) %>%
#     h3::h3_to_geo_boundary_sf() %>%
#     st_make_valid() %>%
#     st_filter(nz, .predicate = st_within)
#   cells <- cells %>%
#     bind_rows(new_cells) %>%
#     st_make_valid()
#   tiled <- cells %>%
#     st_union() %>%
#     st_make_valid()
#   untiled <- nz %>%
#     st_difference(tiled) %>%
#     st_union() %>%
#     st_make_valid()
#   if (res >= 7) break
#   # if (st_area(untiled) < set_units(1e8, "m^2")) break
#   res <- res + 2
# }

tmap_mode("view")

tm_shape(nz) + 
  tm_fill(col = "lightgreen") +
  tm_shape(cells) +
  tm_borders() 



# s2cells <- 
nz %>%
  s2_as_text() %>%
  as_s2_cell_union()

# %>%
#   plot()


    