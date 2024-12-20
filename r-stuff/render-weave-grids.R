library(dplyr)
library(sf)

sf_use_s2(FALSE)
gPRECISION <- 1e8

apply_precision <- function(x, p = gPRECISION) {
  round(x * gPRECISION) / gPRECISION
}


# Returns a matrix of coordinates and a list of the orderings
# of axes at those sites
matrices_as_loom <- function(...) {
  matrices <- list(...)
  if (length(matrices) == 1) {
    mat <- matrices[[1]]
    indices <- expand.grid(weft = seq_len(nrow(mat)),
                           warp = seq_len(ncol(mat))) %>%
      as.matrix()
    orderings <- mat[indices] %>%
      lapply(decode_biaxial_to_order, axis = 0)
    parity <- NULL
    orientations <- c(0, -90)
    dimensions <- c(nrow(mat), ncol(mat))
  } else {
    m_1 <- matrices[[1]]
    m_2 <- matrices[[2]]
    m_3 <- matrices[[3]]

    n_A <- max(ncol(m_1), nrow(m_3))
    n_B <- max(ncol(m_2), nrow(m_1))
    n_C <- max(ncol(m_3), nrow(m_2))
    dimensions <- c(n_A, n_B, n_C)

    m_AB <- m_1 %>%
      pracma::repmat(n = pracma::Lcm(nrow(m_1), n_A) %/% nrow(m_1),
                     m = pracma::Lcm(ncol(m_1), n_A) %/% ncol(m_1))
    m_BC <- m_2 %>%
      pracma::repmat(n = pracma::Lcm(nrow(m_2), n_B) %/% nrow(m_2),
                     m = pracma::Lcm(ncol(m_2), n_B) %/% ncol(m_2))
    m_CA <- m_3 %>%
      pracma::repmat(n = pracma::Lcm(nrow(m_3), n_C) %/% nrow(m_3),
                     m = pracma::Lcm(ncol(m_3), n_C) %/% ncol(m_3))

    parity <- (3 + n_A + n_B + n_C) %/% 2
    orientations <- c(0, 120, 240)
    indices <- expand.grid(a = 1:n_A, b = 1:n_B, c = 1:n_C) %>%
      filter((a + b + c) %in% parity:(parity + 1)) %>%
      as.matrix()
    orderings <-
      mapply(
        combine_orderings,
        m_AB[indices[, 2:1]] %>% lapply(decode_biaxial_to_order, axis = 1),
        m_BC[indices[, 3:2]] %>% lapply(decode_biaxial_to_order, axis = 2),
        m_CA[indices[, c(1, 3)]] %>% lapply(decode_biaxial_to_order, axis = 3),
        SIMPLIFY = FALSE
      )
  }
  list(
    indices = indices,
    orderings = orderings,
    parity = parity,
    orientations = orientations,
    dimensions = dimensions)
}


# Returns a function that will generate the x-y coordinates for
# a supplied set of integer coordinates in a particular 'space'
# n_axes = 2 is a Cartesian grid with spacing S and reverse x-y order
# n_axes = 3 is a Nagy triangular grid with spacing S, and axes in
# vertical, down 120 SW and down 120 SE directions
# See: Nagy, B. N. 2003. Shortest Paths in Triangular Grids with
#      Neighbourhood Sequences. Journal of Computing and Information
#      Technology 11 (2):111.
grid_generator <- function(n_axes = 2, S = 1) {
  if (n_axes == 2) {
    ### NOTE that this reverse x-y coordinates (to match matrix convention)
    angles <- 1:0 * pi / 2
    dx <- S * cos(angles)
    dy <- S * sin(angles)
  }
  if (n_axes == 3) {
    angles <- seq(3, 11, 4) * pi / 6
    dx <- S * cos(angles) * 2 / 3 # every 3rd site missing from coordinate space
    dy <- S * sin(angles) * 2 / 3
  }
  basis <- matrix(c(dx, dy), nrow = 2, ncol = n_axes, byrow = TRUE)
  function(coords) {
    t(basis %*% coords) %>% c()
  }
}

# utility function to rotate sf shape through angle in degrees
# about a centre point as supplied
rotate_shape <- function(shape, angle, centre = c(0, 0)) {
  a <- angle * pi / 180
  m <- t(matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2))
  ((shape %>% translate_shape(-centre)) * m) %>%
    translate_shape(centre)
}

# utility function to translate an sf shape by the
# supplied displacement vector
translate_shape <- function(shape, dxdy = c(0, 0)) {
  shape + dxdy
}

# Returns a grid cell polygon centred at (0, 0) with
# the number of sides requested. One side of the polygon
# will lie horizontal below the x-axis symmetric about the
# y-axis. This means that radii to the corners of the polygon
# are arranged as either
#
# n even:      or, n odd:
#   \  /            |
# ___\/___       ___|___
#    /\            / \
#   /  \          /   \
# The parameter L is the face to face distance of the polygon from
# the base edge vertically to the opposite face (or if n is odd to
# the opposite corner). This means L is either:
# n even: L = 2Rcos(pi/n), or
# n odd:  L = R + Rcos(pi/n)
# where R is the radius of the circumcircle
#
# The polygon is generated by finding the n points equally spaced
# on this circumcircle
get_grid_cell_polygon <- function(face_to_face_distance = 1,
                                  n_sides = 4, parity = 0) {
  # Determine value of R
  if (n_sides %% 2 == 0) {
    R <- face_to_face_distance / (2 * cos(pi / n_sides))
  } else {
    R <- face_to_face_distance / (1 + cos(pi / n_sides))
  }
  # determine angles
  # we start at 6 o'clock (3pi/2), then add (pi/n), then add n more 2pi/n steps
  angles <- (3 * pi / 2) + (pi / n_sides) + (0:n_sides / n_sides) * 2 * pi
  angles[length(angles)] <- angles[1]
  # the corners are the cos and sin of these angles
  corners <- R * c(cos(angles), sin(angles)) %>%
    matrix(nrow = n_sides + 1, ncol = 2)
  polygon <- corners %>% list() %>% st_polygon()
  if (n_sides == 4 || parity %% 2 == 1) {
    return(polygon %>% st_sfc()) #precision = gPRECISION))
  } else {
    return(rotate_shape(polygon, 180) %>% st_sfc()) #precision = gPRECISION))
  }
}

# Returns 'slices' across a grid cell (i.e. horizontally) centered vertically
# relative to the cell, ie
#
#            /\
#           /  \
#   +------------------+
#   |     /      \     |
#   +------------------+
#   |   /          \   |
#   +------------------+
#     /              \
#    /________________\
#
# Horizontal extent is given by L, total width of the strips is W, they are
# 'sliced' horizontally in n_equal slices. An offset should be provided to
# center the slices vertically on the vertical extent of the cell (not its)
# centroid. This is supplied from get_cell_strands()
get_grid_cell_slices <- function(L = 1, W = 1, n_slices = 1, offset = c(0, 0)) {
  sW <- W / n_slices # slice widths
  # odd numbers from 1 to 2n-1
  odd_numbers <- seq(1, (2 * n_slices - 1), 2)
  slice_offsets <-  sW * odd_numbers / 2 - W / 2
  slices <- list()
  for (i in seq_along(slice_offsets)) {
    # L by W rectangle centred at 0,0
    slices[[i]] <- (matrix(0.5 * c(-L, -sW, L, -sW, L, sW, -L, sW, -L, -sW),
                           5, 2, byrow = TRUE) %>%
                      list() %>%
                      st_polygon()) %>%
      translate_shape(c(0, slice_offsets[i])) %>%
      translate_shape(offset)
  }
  slices %>% st_sfc() #precision = gPRECISION)
}


# Gets the cross grid cell strands running across a cell in the x direction
# optionally rotated by orientation for a grid cell spacing S, total strand
# width width (as a fraction of S), sliced into n_slices along its length
get_cell_strands <- function(n = 4, S = 1, width = 1, parity = 0,
                             orientation = 0, n_slices = 1) {
  W <- width * S
  # make expanded cell that reaches to the strands in neighbours 
  big_s <- ifelse(n == 4, 
                  S + S * (1 - width), 
                  S * (5 - 3 * width) / 2)
  expanded_cell <- get_grid_cell_polygon(face_to_face_distance = big_s, 
                                         n = n, parity = parity)
  bb <- st_bbox(expanded_cell)
  strand_offset <- c(bb$xmin + bb$xmax, bb$ymin + bb$ymax) / 2
  cell <- get_grid_cell_polygon(face_to_face_distance = S,
                                n = n, parity = parity)
  # determine its x-y centre (which may not be where its centroid is)
  bb <- st_bbox(cell)
  cell_offset <- c(bb$xmin + bb$xmax, bb$ymin + bb$ymax) / 2
  big_l <- ifelse(n == 4, big_s, 
                  big_s * 2 / sqrt(3) * (3 - width) / 2)
  get_grid_cell_slices(L = big_l, W = W, n_slices = n_slices, 
                       offset = strand_offset) %>%
    lapply(translate_shape, dxdy = -strand_offset + cell_offset) %>%
    st_sfc() %>% #precision = gPRECISION) %>% 
    st_intersection(expanded_cell) %>%
    lapply(rotate_shape, angle = orientation) %>%
    st_sfc() #precision = gPRECISION)
}


add_shapes_to_list <- function(lst, shapes) {
  for (s in shapes) {
    lst <- append(lst, list(s))
  }
  lst
}


# Essentially a wrapper for get_cell_strands that returns the strands in all
# the requested cross directions
get_all_cell_strands <- function(n = 4, S = 1, width = 1, parity = 0,
                                 orientations = c(0, 90),
                                 n_slices = rep(1, length(orientations))) {
  polys <- list()
  for (i in seq_along(orientations)) {
    next_strands <- get_cell_strands(n = n, S = S, width = width,
                                     parity = parity,
                                     orientation = orientations[i],
                                     n_slices = n_slices[i])

    polys <- add_shapes_to_list(polys, next_strands)
  }
  polys %>% st_sfc() #precision = gPRECISION)
}

# Returns the visible parts of the strands in a grid, given the spacing S
# strand width width, parity (for the triangular case), a vector of strand
# orders and matching vectors of orientations and the desired number of slices
get_visible_cell_strands <- function(n = 4, S = 1, width = 1, parity = 0,
                                     strand_order = 1:(6 - n),
                                     orientations = (0:(n - 1)) * 360 / n,
                                     n_slices = rep(1, length(orientations))) {
  all_polys <- list()
  for (i in seq_along(strand_order)) {
    next_polys <- get_cell_strands(n = n, S = S, width = width, parity = parity,
                                   orientation = orientations[strand_order[i]],
                                   n_slices = n_slices[strand_order[i]])
    if (i == 1) {
      all_polys <- add_shapes_to_list(all_polys, next_polys)
      # mask poly progressively builds the union of all polygons
      # so far, to mask out invisible parts of those underneath
      mask_poly <- next_polys %>%
        st_sf() %>% #precision = gPRECISION) %>% 
        st_union()
    } else {
      all_polys <- add_shapes_to_list(all_polys, next_polys %>%
                                        st_snap(mask_poly, 1 / gPRECISION) %>%
                                        st_difference(mask_poly))
      mask_poly <- mask_poly %>%
        st_union(next_polys %>%
                   st_sf() %>% #precision = gPRECISION) %>%
                   st_union())
    }
    # if the width is 1 then no lower polygons are visible
    if (width == 1) break # for efficiency?
  }
  all_polys %>% st_sfc() #precision = gPRECISION)
}

# returns a rectangular polygon matching a provided bounding box
sfc_from_bbox <- function(bb, crs) {
  st_polygon(
    list(matrix(c(bb$xmin, bb$ymin, bb$xmax, bb$ymin, bb$xmax, bb$ymax,
                  bb$xmin, bb$ymax, bb$xmin, bb$ymin), 5, 2, byrow = TRUE))) %>%
    st_sfc(crs = crs) #, precision = gPRECISION)
}

# determines translation vector required to centre shape on centre
centre_offset <- function(shape, centre = c(0, 0)) {
  bb <- st_bbox(shape)
  centre - c((bb$xmax + bb$xmin) / 2, (bb$ymax + bb$ymin) / 2)
}


# builds the sf associate with a given weave supplied as 'loom' which is a list
# containing the coordinates in an appropriate grid (Cartesian or triangular)
# and the orderings of the strands at each coordinate location
make_sf_from_coded_weave_matrix <- function(loom, spacing = 1, width = 1,
                                            margin = 0,
                                            axis1_threads = letters[1],
                                            axis2_threads = letters[2],
                                            axis3_threads = letters[3],
                                            crs = 3857) {
  # we need number of axes to make a grid generator function
  n_axes <- length(loom$dimensions)
  n_sides <- if_else(n_axes == 2, 4, 3)
  gg <- grid_generator(n_axes = n_axes, S = spacing)
  # the labels for axis1 and 2 are required in both cases and we repeat
  # them if required by the size of the grid relative to the number of ids
  ids1 <- rep(axis1_threads, loom$dimensions[1] / length(axis1_threads))
  ids2 <- rep(axis2_threads, loom$dimensions[2] / length(axis2_threads))
  if (n_axes == 3) {
    ids3 <- rep(axis3_threads, loom$dimensions[3] / length(axis3_threads))
  }
  parity <- 1
  # setup empty lists and vectors for the outputs
  weave_polys <- list()
  bb_polys <- list()
  strands <- c()
  # step through the loom index coordinates
  for (i in seq_len(nrow(loom$indices))) {
    coords <- loom$indices[i, ]
    ids <- c(ids1[coords[1]], ids2[coords[2]])
    if (n_axes == 3) { # triaxial - extra threads and also have to set parity
      ids <- c(ids, ids3[coords[3]])
      parity <- (sum(coords) - loom$parity) %% 2
    }
    # get the offset vector
    xy <- gg(coords)
    strand_order <- loom$orderings[[i]]      # order of strands from the top
    bb_polys[[i]] <- 
      get_grid_cell_polygon(face_to_face_distance = spacing, 
                            n_sides = n_sides, parity = parity) %>%
      translate_shape(xy)
    if (is.null(strand_order)) next          # nothing here so move on
    if (anyNA(strand_order)) {               # NA areas to add
      weave_polys <- 
        add_shapes_to_list(weave_polys, 
                           get_grid_cell_polygon(
                             face_to_face_distance = spacing, 
                             n_sides = n_sides, parity = parity) %>% 
                             translate_shape(xy))
      strands <- c(strands, "NA")
      warning(paste("Impossible to determine strand order at:", 
                   paste(coords, collapse = ","), collapse = " "))
      next
    }
    n_slices <- stringr::str_length(ids)   # number of slices in each direction
    next_polys <- get_visible_cell_strands(n = n_sides, S = spacing,
                                           width = width,
                                           strand_order = strand_order,
                                           parity = parity,
                                           orientations = loom$orientations,
                                           n_slices = n_slices)
    # add polygons one at a time to simplify later conversion to sfc
    # make up labels as a string in the order they'll be needed
    # TODO: this isn't working right just yet - when n_slices > 1 it gets
    # the segments that show in the lower layers wrong quite often - presumably
    # get_visible_cell_strands() does not guarantee the returned order
    labels <- ids[strand_order] %>% paste0(collapse = "")
    for (p in seq_along(next_polys)) {
      weave_polys <- append(weave_polys, 
                            list(next_polys[[p]] %>% translate_shape(xy)))
      strands <- c(strands, stringr::str_sub(labels, p, p))
    }
  }
  tile <- bb_polys %>% 
    sapply("[") %>%                 # convert to sfc then to sp
    st_sfc() %>%                    # so that we can use rmapshaper::  
    as("Spatial") %>%               # ms_dissolve - because sf 
    rmapshaper::ms_dissolve() %>%   # group_by sucks 
    st_as_sf() 
  # %>%                             # back to sf
  #   st_set_precision(gPRECISION)
  tile_centre_offset <- centre_offset(tile)
  tile <- wk::wk_transform(tile, wk::wk_affine_translate(tile_centre_offset[1],
                                                         tile_centre_offset[2]))
  list(
    weave_unit = weave_polys %>%
      st_as_sfc() %>%                          # convert to sfc
      wk::wk_transform(wk::wk_affine_translate(tile_centre_offset[1],
                                               tile_centre_offset[2])) %>%
      st_sf() %>%
      mutate(strand = strands) %>%             # add the strands information
      filter(strand != "-") %>%                # remove any tagged missing
      as("Spatial") %>%                        # convert to sp for rmapshaper 
      rmapshaper::ms_dissolve(field = "strand") %>%
      st_as_sf() %>%                           # back to sf
      st_buffer(-margin * spacing) %>%         # include a negative margin
      st_intersection(tile) %>%                # cookie cut to tile
      # st_set_precision(gPRECISION) %>%
      st_set_crs(crs),                         # set CRS
    tile = tile %>%
      st_set_crs(crs)                          # set CRS
  )
}

