require(sf)

hex_types <- c("hex", "cube")

# delegates creation of tile offsets filling the region to
# sf::st_make_grid
get_centres <- function(region, tile, hexes) {
  bb <- st_bbox(tile)
  w <- bb$xmax - bb$xmin
  h <- bb$ymax - bb$ymin
  region_b <- region %>% st_union() %>% st_buffer(max(w, h))
  if (hexes) {
    pts <- region_b %>% st_make_grid(cellsize = h, flat_topped = TRUE,
                                     what = "centers", square = FALSE)
  } else {
    pts <- region_b %>% st_make_grid(cellsize = c(w, h),
                                     what = "centers")
  }
  pts %>%
    st_as_sf() %>%
    st_filter(region_b) %>%
    st_coordinates() %>%
    ## next steps make a list of pairs so we can lapply for speed
    t() %>%
    as.data.frame() %>%
    as.list()
}

## could use sf_transform for this, but just adding
## the point is much more convenient and readily compatible
## with lapply approach to the tiling
sf_shift <- function(pt, shapes) {
  st_geometry(shapes) <- st_geometry(shapes) + pt
  shapes
}


weave_layer <- function(
    weave_unit, region, angle = 0, tile = weave_unit$tile,
    transform = wk::wk_affine_identity(),
    merge_map_polys_in_data_polys = TRUE, # helps with leaflet, but caution
    region_cols_to_summarize_by_strand = c()) {

  if (merge_map_polys_in_data_polys &&
      (length(region_cols_to_summarize_by_strand) > 0)) {
    warning(stringr::str_c(
      "Merging map polygons within data polygons not supported when\n",
      "   adding data columns that summarize values across the strand unit.\n",
      "Here, summarizing data while NOT merging map polygons as requested."))
    merge_map_polys_in_data_polys <- FALSE
  }

  cntrd <- region %>% sf_get_centroid()
  cx <- cntrd[1]
  cy <- cntrd[2]

  to_tile <- region  %>%
    sf_transform(wk::wk_affine_invert(transform)) %>%
    sf_rotate(-angle, cx, cy) %>%
    sf_transform(weave_unit$transform) %>%
    dplyr::mutate(to_tile_id = row_number())

  the_unit <- weave_unit$primitive %>%
    sf_transform(weave_unit$transform)
  the_tile <- tile %>% 
    sf_transform(weave_unit$transform)

  pts <- get_centres(to_tile, tile = the_tile,
                     hex = weave_unit$type %in% hex_types)

  tiling <- lapply(pts, sf_shift, shapes = the_unit) %>%
    bind_rows() %>%
    mutate(strand_id = row_number()) %>%
    st_set_crs(st_crs(region)) %>%
    st_intersection(to_tile)

  if (length(region_cols_to_summarize_by_strand) > 0) {

    to_count_up <- tiling %>%
      mutate(area = st_area(.)) %>%
      st_drop_geometry() %>%
      as_tibble() %>%
      group_by(strand_id)

    for (col_to_summarize in region_cols_to_summarize_by_strand) {
      tiling <- tiling %>%
        left_join(
          to_count_up %>%
            count(!!as.name(col_to_summarize), wt = area, sort=TRUE) %>%
            top_n(1, n) %>%
            select(-n) %>%
            rename(!!stringr::str_c(col_to_summarize, "_by_strand"):=!!col_to_summarize),
          by = "strand_id"
        )
    }
  }

  if (merge_map_polys_in_data_polys) {
    tiling <- tiling %>%
      qgis::qgis_dissolve(FIELD = c("strand", "to_tile_id")) %>%
      st_as_sf() %>%
      st_set_crs(st_crs(region))
  }

  tiling %>%
    sf_transform(wk::wk_affine_invert(weave_unit$transform)) %>%
    sf_rotate(angle, cx, cy) %>%
    sf_transform(transform)
}


write_weave_layers <- function(weave, region, fname, var = "strand") {
  st_write(region, fname, "basemap", delete_dsn = TRUE)
  lyrs <- split(weave, weave[[var]])
  for (label in names(lyrs)) {
    lyrs[[label]] %>% st_write(fname, label, append = TRUE)
  }
  invisible(weave)
}
