library(tidyverse)
library(sf)

# credit to Jens von Bergmann for this algo https://github.com/mountainMath/dotdensity/blob/master/R/dot-density.R
random_round <- function(x) {
  v <- as.integer(x)
  r <- x - v
  test = runif(length(r), 0.0, 1.0)
  add = rep(as.integer(0), length(r))
  add[r > test] <- as.integer(1)
  value =v + add
  ifelse(is.na(value) | value < 0,0, value)
  return(value)
}

num_dots <- as.data.frame(region) %>% 
  select(European:Asian) %>% 
  mutate(across(where(is.numeric), ~ round(. / 25)))

# generates data frame with coordinates for each point + what group it is assiciated with
sf_dots <- map_df(names(num_dots), 
                  ~ st_sample(region, size = num_dots[, .x], type = "random") %>% # generate the points in each polygon
                    st_cast("POINT") %>%                                          # cast the geom set as 'POINT' data
                    st_coordinates() %>%                                          # pull out coordinates into a matrix
                    as_tibble() %>%                                               # convert to tibble
                    setNames(c("x", "y")) %>%                                     # set column names
                    mutate(group = .x)                                            # add categorical party variable
) %>% 
  slice(sample(1:n())) # once map_df binds rows randomise order to avoid bias in plotting order

dots <- sf_dots %>% 
  st_as_sf(coords = c("x", "y"), crs = 2193)

tmap_mode("plot")
basemap <- read_osm(st_bbox(region, zoom = 7, type = "stamen-toner"), raster = TRUE)

# tm_shape(basemap) + 
#   tm_rgb() + 
tm_shape(region) +
  tm_borders() +
  tm_shape(dots) + 
  tm_dots(col = "group", size = 0.04, border.lwd = 0, 
          palette = c("#31a354", "#252525", "#de2d26", "#756bb1"))
