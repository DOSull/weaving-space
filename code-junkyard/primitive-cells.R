require(tidyr)
require(stringr)
require(wk)
require(sf)
require(pracma)

S3 <- sqrt(3)


parse_labels <- function(ids) {
  lbls <- strsplit(ids, "|", fixed = TRUE)
  labels_1 <- lbls[[1]][1] 
  if (length(lbls[[1]]) > 1) {
    labels_2 <- lbls[[1]][2]
  } else {
    labels_2 <- "-"
  }
  if (length(lbls[[1]]) > 2) {
    labels_3 <- lbls[[1]][3]
  } else {
    labels_3 <- "-"
  }
  return(c(labels_1, labels_2, labels_3))
}

# these functions in primitive-cells.R
get_primitive_cell <- function(width, spacing, margin, ids, type,
                               square, n_twill, crs) {

  pl <- parse_labels(ids)

  cell <- switch(
    type,
    "open" = getPlainPrimitiveCell(
      width = width, spacing = spacing, margin = margin, 
      labels_1 = pl[1], labels_2 = pl[2], 
      square = square, crs = crs),
    
    "diamond" = getDiamondPrimitiveCell(
      width = width, spacing = spacing, margin = margin, 
      labels_1 = pl[1], labels_2 = pl[2], labels_3 = pl[3],
      crs = crs),
    
    "hex" = getHexPrimitiveCell(
      width = width, spacing = spacing, margin = margin, 
      labels_1 = pl[1], labels_2 = pl[2], labels_3 = pl[3],
      crs = crs)
  )
  return(cell)
}

## for local rotations of single geoms, this is more convenient
## than the wk::wk_ affine transforms option
get_rot_matrix <- function(angle) {
  a <- angle * pi / 180
  return(t(matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2, byrow = FALSE)))
}


# Makes a primitive cell for an orthogonal weave made up of
# alternating horizontal and vertical rectangles, the basic unit being
#    
#   ||==
#   ==||
#
# width is the width of the 'ribbons' in the weave
# spacing is the spacing from centre line to centre line of ribbons
# margin is an 'inset' may enhance visual impression of a weave
# labels_1 are the distinct labels in one direction, labels2 in the other
# the length of the labels vectors will determine the repeat of the tiling
# e.g. if labels_1 is length 1 and labels_2 length 2, we get
#    
#   ||==||==
#   ==||==||
#
# square is experimental TRUE will make the cell square regardless of
# the repeats in the two directions, and changing widths, spacings, and lengths
# accordingly
# crs is any sf compatible definition of a CRS
getPlainPrimitiveCell <- function(width, spacing, aspect, margin, 
                                  labels_1, labels_2, square = FALSE, crs) {

  lbls_1 <- strsplit(labels_1, "", fixed = TRUE)[[1]]
  lbls_2 <- strsplit(labels_2, "", fixed = TRUE)[[1]]
  short <- min(length(lbls_1), length(lbls_2))
  spacing_h <- ifelse(square, spacing / length(lbls_1) * short, spacing)
  spacing_v <- ifelse(square, spacing / length(lbls_2) * short, spacing)
  width_h <- ifelse(square, width / length(lbls_1) * short, width)
  width_v <- ifelse(square, width / length(lbls_2) * short, width)
  
  # calculate the lengths from the gaps in the transverse direction
  length_h <- width_h / aspect #2 * spacing_v - width_v
  length_v <- width_v / aspect #2 * spacing_h - width_h
  # make base polygons for the two directions, centred on (0, 0)
  base_poly_h <- get_base_rect(length_h, width_h, margin)
  base_poly_v <- get_base_rect(length_v, width_v, margin, orientation = "vertical")
  # determine the number of ribbons needed in each
  # direction in the 'base tiling' from which we 
  # will cut the repeating unit (primitive cell)
  reps_1 <- ifelse(length(lbls_1) %% 2 == 0, 1, 2)
  reps_2 <- ifelse(length(lbls_2) %% 2 == 0, 1, 2)
  n_rows <- reps_1 * length(lbls_1) + 2
  n_cols <- reps_2 * length(lbls_2) + 2
  # empty list for the polygons and vector for the labels
  polys <- list()
  ids <- c()
  item <- 1
  for (row in 0:(n_rows - 1)) {
    dy <- row * spacing_h
    for (col in 0:(n_cols - 1)) {
      dx <- col * spacing_v 
      # alternate adding horizontal and vertical polygons
      # at the offsets from the base given by dx and dy
      if ((col + row) %% 2 == 0) {
        this_label <- lbls_1[row %% length(lbls_1) + 1]
        if (this_label == "-") {
          polys[[item]] <- base_poly_v + c(dx, dy)
          ids <- c(ids, lbls_2[col %% length(lbls_2) + 1])
        } else {
          polys[[item]] <- base_poly_h + c(dx, dy)
          ids <- c(ids, this_label)
        }
      } 
      else {
        this_label <- lbls_2[col %% length(lbls_2) + 1]
        if (this_label == "-") {
          polys[[item]] <- base_poly_h + c(dx, dy)
          ids <- c(ids, lbls_1[row %% length(lbls_1) + 1])
        } else {
          polys[[item]] <- base_poly_v + c(dx, dy)
          ids <- c(ids, this_label)
        }
      }
      item <- item + 1
    }
  }
  all_tiles <- st_sf(id = ids, geometry = st_as_sfc(polys)) %>%
    filter(id != "-") %>%
    group_by(id) %>%
    dplyr::summarise()
  bb <- st_bbox(c(xmin = 0, ymin = 0,
                  xmax = reps_2 * length(lbls_2) * spacing_v,
                  ymax = reps_1 * length(lbls_1) * spacing_h))
  return(list(
    cell = all_tiles %>% 
      st_crop(bb) %>% 
      st_set_crs(crs),
    transform = wk::wk_affine_identity(),
    ids = unique(all_tiles$id)
  ))
}


# Returns a diamond weave primitive cell that is effectively equivalent
# to getHexPrimitveCell as currently written
getDiamondPrimitiveCell <- function(width, spacing, margin, 
                                    labels_1, labels_2, labels_3, crs) {

  lbls <- c(labels_1, labels_2, labels_3)
  spacings <- sapply(lbls, stringr::str_length)
  s_mult <- max(spacings)

  spacing_x <- s_mult * max(2 * S3, ceiling(spacing / width)) * width
  length <- spacing_x - width * 2 / S3
  base_p <- get_polygon(c(0, 0, length, 0,
                          length - width / S3, width,
                          -width / S3, width))
  base_polys <- list(base_p, 
                     base_p * get_rot_matrix(-120),
                     base_p * get_rot_matrix(120))
  translations <- c(0, 0, spacing_x, -spacing_x * S3,
                    2 * spacing_x, 0, spacing_x, spacing_x * S3) / 2
  polys <- list()
  item <- 1
  for (p in base_polys) {
    for (i in seq(1, 8, 2)) {
      polys[[item]] <- p + translations[i:(i + 1)]
      item <- item + 1
    }
  }
  ids <- c(rep(labels_1, 4), rep(labels_2, 4), rep(labels_3, 4))
  prototile <- get_polygon(translations) ## the diamond
  return(list(
    cell = st_sf(id = ids, geometry = st_as_sfc(polys)) %>% 
      st_buffer(-margin / 2) %>%
      st_intersection(prototile) %>%
      st_as_sf() %>%
      filter(st_geometry_type(.) == "POLYGON" & st_area(.) > 1e-10) %>%
      st_set_crs(crs),
    # the transform required to make this rectangular tile-able
    transform = wk::wk_affine_invert(
      affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)),
    ids = unique(ids)
  ))
}


getHexPrimitiveCell <- function(width = 200, spacing = 300, margin = 0,
                                labels_1, labels_2, labels_3, crs = 3857) {
  ## type == "hex" (which makes same weave as "diamond"... need to think this through)
  S3 <- sqrt(3)
  spacing_x <- max(spacing / width, S3 * 2) * width 
  length <- spacing_x - width * 2 / S3
  
  base_p1 <- get_polygon(c(0, 0, -length, 0, 
                           -length + width / S3, -width,
                           width / S3, -width))
  base_p2 <- base_p1 + c(length + 2 * width / S3, 0)
  
  polys <- list(base_p1, base_p1 * get_rot_matrix(-120), base_p1 * get_rot_matrix(120),
                base_p2, base_p2 * get_rot_matrix(-120), base_p2 * get_rot_matrix(120))
  ids <- rep(c(labels_1, labels_2, labels_3), 2)
  
  r <- sqrt((length - width / S3) ^ 2 + width ^ 2)
  angles <- seq(1, 12, 2) / 6 * pi
  dxdy <- c(matrix(c(cos(angles), sin(angles)), ncol = 6, byrow = TRUE)) * r
  prototile <- get_polygon(dxdy)
  
  return(list(
    cell = st_sf(id = ids, geometry = st_as_sfc(polys)) %>%
      st_buffer(-margin / 2) %>%
      st_intersection(prototile) %>%
      st_as_sf() %>%
      filter(st_geometry_type(.) == "POLYGON" & st_area(.) > 1e-10) %>%
      st_set_crs(crs),
    transform = wk::wk_affine_identity(),
    ids = unique(ids)
  ))
}




get_weave_unit <- function(spacing = 300, aspect = 1, width = 200, margin = 0,
                       ids = "a|b|c", type = "open", 
                       n_twill = 2, crs = 3857) {

  width = ifelse(type == "open", 2 * spacing * aspect / (1 + aspect), width)

  pc <- get_primitive_cell(
    spacing = spacing, aspect = aspect, width = width, margin = margin,
    ids = ids, type = type, square = square, n_twill = n_twill, crs = crs)
  
  return(
    list(
      primitive = pc$cell,
      transform = pc$transform,
      ids = pc$ids,
      type = type
    )
  )
}



