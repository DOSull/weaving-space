library(sf)
library(dplyr)
library(anthromes)

source("weaving-space-utils.R")
source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weave-map.R")


# nz <- st_read("~/Documents/geodata/nz-2193.shp") %>% 
#   st_transform(4326)
# 
# an_baseline <- read.csv("anthromes/an12_dgg_baseline.csv") %>%
#   select(1, 39, 49, 59)
# 
# anthromes_dgg <- st_read("anthromes/an12_dgg_inputs.shp") %>%
#   st_filter(nz) %>%
#   st_transform(2193) %>% 
#   left_join(an_baseline) %>% 
#   left_join(anthrome_key, by = c("X1800AD" = "anthrome")) %>%
#   rename(class1800 = class, level1800 = level, type1800 = type) %>%
#   left_join(anthrome_key, by = c("X1900AD" = "anthrome")) %>%
#   rename(class1900 = class, level1900 = level, type1900 = type) %>%
#   left_join(anthrome_key, by = c("X2000AD" = "anthrome")) %>%
#   rename(class2000 = class, level2000 = level, type2000 = type)
# 
# 
# st_write(anthromes_dgg, "nz-anthromes-dgg.gpkg")

anthromes_dgg <- st_read("nz-anthromes-dgg.gpkg")

library(tmap)
tmap_mode("plot")

tm_shape(anthromes_dgg) + 
  tm_fill(col = paste("class", 18:20 * 100, sep = ""), 
          palette = anthrome_colors(), title = "Anthrome") +
  tm_layout(panel.labels = c("1800", "1900", "2000"), 
            legend.outside = TRUE)


st_write(anthrome_weave, "anthromes/anthrome-weave-cube.gpkg")

twill_tile <- get_biaxial_weave_unit(spacing = 5000, margin = 100, 
                                     aspect = 0.8, type = "twill", n = 2,
                                     ids = "a|b", crs = 2193)
plain_tile <- get_biaxial_weave_unit(spacing = 5000, margin = 100, 
                                     aspect = 0.65, ids = "a|b", crs = 2193)
cube_tile2 <- get_triaxial_weave_unit(spacing = 8000, type = "cube",
                                      ids = "a|b|c", crs = 2193)

anthrome_weave <- weave_layer(cube_tile2, anthromes_dgg, angle = 45)
anthrome_layers <- split(anthrome_weave, anthrome_weave$id)

tmap_mode("view")
tmap_options(check.and.fix = TRUE)
tm_shape(anthrome_layers$a) + 
  tm_fill(col = "class1800", palette = anthrome_colors(), 
          legend.show = FALSE) + 
  tm_shape(anthrome_layers$b) + 
  tm_fill(col = "class1900", palette = anthrome_colors(), 
          legend.show = FALSE) + 
  tm_shape(anthrome_layers$c) + 
  tm_fill(col = "class2000", palette = anthrome_colors(), 
          title = "Anthrome") + 
  tm_layout(legend.outside = TRUE, 
            main.title = "Transitions from 1800 to 2000")


####################################x

# library(anthromes)
# library(sf)
# require(terra)
# require(stars)
# require(rgdal)
# require(gdalcubes)
# library(ggplot2)
# 
# setwd("~/Documents/code/weaving-space")
# 
# nz <- st_read("~/Documents/geodata/nz-2193.gpkg") %>%
#   st_transform(4326)
# bbox <- nz %>%
#   st_bbox()
# 
# data("hyde_med")
# data("inputs_med")
# 
# anthromes <- anthrome_classify(hyde_med, inputs_med)
# 
# ggplot() +
#   geom_stars(data = anthromes) +
#   facet_wrap(~time) +
#   coord_quickmap() +
#   scale_fill_manual(values = anthrome_colors(), drop = TRUE) +
#   theme_bw() +
#   labs(title = 'Anthromes-12k', x = 'Latitude', y = 'Longitude')
# 
# w_inputs <- hyde_read(dir = "/home/osullid3/Downloads/Anthromes-12k-DGG", 
#                       vars = c("inputs")) %>%
#   st_crop(bbox)
# 
# 
# folder <- "/home/osullid3/Downloads/Anthromes-12k-DGG/raw-data/raw-data/HYDE/HYDE"
# tifs <- dir(path = folder, pattern = "*.tif$")
# 
# lyrs <- list()
# for (tif in tifs) {
#   lyrname <- strsplit(tif, ".", fixed = TRUE)[[1]][1]
#   lyr <- rast(file.path(folder, tif))[["1800AD"]]
#   names(lyr) <- lyrname
#   lyrs[[lyrname]] <- rast(file.path(folder, tif))[["1800AD"]]
# }
# 
# r1800 <- rast(lyrs) %>%
#   terra::crop(terra::ext(nz))
# 
# 
# hyde_nz <- hyde_read(dir = '..') %>%
#   st_crop(bbox) %>%
#   # current bug in stars loses offset info when filtered
#   #filter(time %in% time_steps_millennia) %>% 
#   # so subset instead
#   .[,,,c(38, 48, 58)] %>%
#   st_as_stars() %>% 
#   # change the names to something more readable
#   setNames(c('crops', 'grazing', 'rice', 'pop', 'irrigation', 'urban'))
# 
