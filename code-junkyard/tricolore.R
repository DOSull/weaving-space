library(tricolore)
library(ggplot2)
library(ggtern)

region <- st_read("../data/akcity-tb-06.gpkg")

eth_mix <- Tricolore(
  region, p1 = "EUR_P_06", p2 = "MAO_P_06", p3 = "ASI_P_06", breaks = Inf
)
region$eth_mix_tri <- eth_mix$rgb

ggplot(region) + 
  geom_sf(aes(fill = eth_mix_tri)) + 
  scale_fill_identity() +
  annotation_custom(
    ggplotGrob(eth_mix$key + labs(L = 'Pakeha', T = 'Maori', R = 'Asian')),
    xmin=1.748e6, xmax=1.753e6, ymin=-Inf, ymax=Inf) +
  coord_sf(datum = NA)
