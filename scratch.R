sf_use_s2(FALSE)

centers <- get_centres(anthromes_dgg, cube_tile2$primitive, TRUE) # 6694

z <- lapply(centers, sf_shift, shapes = cube_tile2$primitive) # 3 times above
z2 <- z %>% bind_rows()
z3 <- z2 %>% mutate(tile_id = row_number())
z3
z4 <- z3 %>% st_set_crs(2193) %>% st_cast()
z4
z5 <- z4 %>% st_intersection(anthromes_dgg) 
# z5 %>% filter(id.1 == 4716062) %>% plot()
z6 <- z5 %>% group_by(id.1)
 
z7 <- z6 %>% summarise(.groups = "keep")

tm_shape(z7) + tm_polygons(col = "cell2000")

