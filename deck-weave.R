library(classInt)
library(RColorBrewer)
library(deckgl)
library(hash)

x <- st_read("data/ak-weave.json") # this made from the auckland IMD data

lookup = hash("a" = list(var = "Rank_Emplo", pal = "BrBG"), 
              "b" = list(var = "Rank_Incom", pal = "PiYG"), 
              "c" = list(var = "Rank_Crime", pal = "RdBu"), 
              "d" = list(var = "Rank_Housi", pal = "RdGy"), 
              "e" = list(var = "Rank_Healt", pal = "PRGn"), 
              "f" = list(var = "Rank_Educa", pal = "RdYlBu"), 
              "g" = list(var = "Rank_Acces", pal = "PuOr"))

x$colour <- 0
for (strand in letters[1:7]) {
  info <- lookup[[strand]]
  var <- info$var
  palname <- info$pal
  ci <- classIntervals(x[[var]], n = 9)
  cols <- findColours(ci, pal = brewer.pal(9, palname) %>% rev())
  x$colour[which(x$strand == strand)] <- (cols[1:nrow(x)])[which(x$strand == strand)]
}

deck <- deckgl(latitude = -36.85, longitude = 174.75, 
               zoom = 10, pickingradius = 5) %>% 
  add_geojson_layer(data = x, 
                    filled = TRUE, 
                    getLineWidth = 0,
                    getFillColor = ~colour) %>% 
  add_basemap()
