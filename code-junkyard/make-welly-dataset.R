library(sf)
library(tidyr)
library(dplyr)
library(tmap)

sa2 <- st_read("../data/sa2-generalised.gpkg") %>%
  filter(TA2018_V1_00_NAME == "Wellington City")

# tm_shape(sa2) + 
#   tm_polygons(col = "population")

A <- runif(nrow(sa2))
B <- runif(nrow(sa2))
Ac <- floor(A * 3)
Bc <- floor(B * 2)

sa2 <- sa2 %>% 
  mutate(A = A, B = B, Ac = Ac, Bc = Bc)

st_write(sa2, "../data/weave-sample-data.gpkg", delete_dsn = TRUE)
