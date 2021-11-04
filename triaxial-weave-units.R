require(tidyr)
require(wk)
require(sf)
require(stringr)


S3 <- sqrt(3)


# these functions in primitive-cells.R
get_primitive_cell <- function(spacing, width, aspect, 
                               margin, ids, type,
                               square, n_twill, crs) {

  labels <- ids %>% 
    parse_labels() %>%
    lapply(string_to_chars)

  cell <- switch(
    type,
    "diamond" = getDiamondPrimitiveCell(
      width = width, spacing = spacing, margin = margin, 
      labels_1 = labels[[1]], labels_2 = labels[[2]], labels_3 = labels[[3]],
      crs = crs),
    
    "hex" = getHexPrimitiveCell(
      width = width, spacing = spacing, margin = margin, 
      labels_1 = labels[[1]], labels_2 = labels[[2]], labels_3 = labels[[3]],
      crs = crs),
    
    "cube" = getCubePrimitiveCell(
      spacing = spacing, margin = margin, 
      labels_1 = labels[[1]], labels_2 = labels[[2]], labels_3 = labels[[3]],
      crs = crs)
  )
  return(cell)
}

## for local rotations of single geoms, this is more convenient
## than the wk_ affine transforms option
get_rot_matrix <- function(angle) {
  a <- angle * pi / 180
  return(t(matrix(c(cos(a), sin(a), -sin(a), cos(a)), 
                  2, 2, byrow = FALSE)))
}


# Returns a diamond weave primitive cell that is effectively equivalent
# to getHexPrimitveCell as currently written
getDiamondPrimitiveCell <- function(width, spacing, margin, 
                                    labels_1, labels_2, labels_3, crs) {

  lbls <- c(labels_1, labels_2, labels_3)
  spacings <- sapply(lbls, str_length)
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
    transform = wk_affine_invert(
      affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)),
    ids = unique(ids)
  ))
}

get_parallelogram <- function(origin, L, W) {
  top_edge <- c(-L, 0)
  left_edge <- c(W / S3, -W)
  bottom_edge <- -top_edge
  return(get_polygon(
    c(origin, origin + top_edge, origin + top_edge + left_edge,
      origin + top_edge + left_edge + bottom_edge)))
}

getHexPrimitiveCell <- function(width = 200, spacing = 300, margin = 0,
                                labels_1 = letters[1:2], labels_2 = letters[3:4],
                                labels_3 = letters[5:6], crs = 3857) {
  ## type == "hex" (which makes same weave as "diamond"... 
  ## need to think this through)
  spacing_x <- max(spacing / width, S3 * 2) * width 
  length <- spacing_x - width * 2 / S3
  
  ns <- c(length(labels_1), length(labels_2), length(labels_3))
  angles <- 0:2 * 120
  polys <- list()
  i <- 1
  for (dir in 1:3) {
    rot <- get_rot_matrix(angles[dir])
    n <- ns[dir]
    for (p in 1:n) {
      p1 <- get_parallelogram(c(width / S3, -width) * (p - 1) / n,
                              length, width / n)
      p2 <- p1 + c(length + 2 * width / S3, 0)
      polys[[i]] <- p1 * rot
      polys[[i + 1]] <- p2 * rot
      i <- i + 2
    }
  }
  ids <- rep(c(labels_1, labels_2, labels_3), 2) %>%
    matrix(sum(ns), 2) %>% t() %>% c()
  r <- sqrt((length - width / S3) ^ 2 + width ^ 2)
  angles <- seq(1, 12, 2) / 6 * pi
  dxdy <- c(matrix(c(cos(angles), sin(angles)), 
                   ncol = 6, byrow = TRUE)) * r
  prototile <- get_polygon(dxdy)

  return(list(
    cell = st_sf(id = ids, geometry = st_as_sfc(polys)) %>%
      # filter(id != "-") %>% # missing threads are kinda tricky
      st_buffer(-margin) %>%
      group_by(id) %>% 
      summarise() %>% 
      st_intersection(prototile) %>%
      st_as_sf() %>%
      st_set_crs(crs),
    transform = wk_affine_identity(),
    ids = unique(ids)
  ))
}


getCubePrimitiveCell <- function(spacing = 300, margin = 0,
                                labels_1 = letters[1:2], labels_2 = letters[3:4],
                                labels_3 = letters[5:6], crs = 3857) {
  L <- spacing
  W <- spacing * S3 / 2
  ns <- c(length(labels_1), length(labels_2), length(labels_3))
  angles <- c(30, 150, 270)
  polys <- list()
  i <- 1
  for (dir in 1:3) {
    rot <- get_rot_matrix(angles[dir])
    n <- ns[dir]
    for (p in 1:n) {
      p1 <- get_parallelogram(c(W / S3, -W) * (p - 1) / n,
                              L, W / n)
      polys[[i]] <- p1 * rot
      i <- i + 1
    }
  }
  ids <- c(labels_1, labels_2, labels_3)
  angles <- seq(1, 12, 2) / 6 * pi
  dxdy <- c(matrix(c(cos(angles), sin(angles)), 
                   ncol = 6, byrow = TRUE)) * L
  prototile <- get_polygon(dxdy)
  
  return(list(
    cell = st_sf(id = ids, geometry = st_as_sfc(polys)) %>%
      # filter(id != "-") %>% # missing threads are kinda tricky
      st_buffer(-margin) %>%
      group_by(id) %>% 
      summarise() %>% 
      st_intersection(prototile) %>%
      st_as_sf() %>%
      st_set_crs(crs),
    transform = wk_affine_identity(),
    ids = unique(ids)
  ))
}


get_triaxial_weave_unit <- function(spacing = 300, aspect = 1, 
                                    width = 200, margin = 0, 
                                    ids = "a|b|c", type = "hex", 
                                    crs = 3857) {

  pc <- get_primitive_cell(
    spacing = spacing, aspect = aspect, width = width, margin = margin,
    ids = ids, type = type, crs = crs)
  
  return(
    list(
      primitive = pc$cell,
      transform = pc$transform,
      ids = pc$ids,
      type = type
    )
  )
}



