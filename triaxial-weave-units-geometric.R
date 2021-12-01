require(sf)

# THIS version stashed in case of a need for future use
# The geometric approach to weave generation has likely application to the
# more geometric challenges of tessellations

# Functions to generate sf datasets that represent tileable to produce the
# appearance of a woven surface, with three directions of weaving thread.
# The methods used are geometric, i.e., there is not underlying mathematical
# or computational representation of a weaving process, unlike in
# biaxial-weave-units.R.


# convenience storage of sqrt(3)
S3 <- sqrt(3)


get_triaxial_weave_unit <- function(spacing = 500, aspect = 1, margin = 0,
                                    strands = "a|b|c", type = "hex", crs = 3857) {

  pc <- get_primitive_cell(
    spacing = spacing, aspect = aspect, margin = margin,
    strands = strands, type = type, crs = crs)

  list(
    primitive = pc$cell,
    transform = pc$transform,
    strands = pc$strands,
    tile = pc$tile,
    type = type
  )
}

# Delegates creation of the tileable unit to an appropriate function
# based on the supplied type "hex", "diamond" or "cube".
get_primitive_cell <- function(spacing = 500, aspect = 1, margin = 0,
                               strands = "a|b|c", type = "hex", crs = 3857) {
  # parse e.g. "a|bc|d" to list(c("a"), c("b", "c"), c("d"))
  labels <- strands %>%
    parse_labels() %>%
    lapply(string_to_chars)

  cell <- switch(
    type,
    "diamond" = getDiamondPrimitiveCell(
      spacing = spacing, aspect = aspect, margin = margin,
      labels_1 = labels[[1]], labels_2 = labels[[2]], labels_3 = labels[[3]],
      crs = crs),

    "hex" = getHexPrimitiveCell(
      spacing = spacing, aspect = aspect, margin = margin,
      labels_1 = labels[[1]], labels_2 = labels[[2]], labels_3 = labels[[3]],
      crs = crs),

    "cube" = getCubePrimitiveCell(
      spacing = spacing, margin = margin,
      labels_1 = labels[[1]], labels_2 = labels[[2]], labels_3 = labels[[3]],
      crs = crs)
  )
  cell
}

## for local rotations of single geoms (around 0, 0), this is more convenient
## than the wk affine transforms option
get_rot_matrix <- function(angle) {
  a <- angle * pi / 180
  t(matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2, byrow = FALSE))
}


# Returns a hexagonally tileable unit that produces a triangular weave
# optionally with more than one thread in any of the three directions. The
# base tile is something like
#
#      /    /
#     /    / \
#  __/____/   \___
#         \    \
#  ________\    \_
#           \    \
#
# where the tileable shape is the 'points up' hexagon enclosing this
# spacing, is the repeat width of the hexagon, aspect is the (normal) width
# of the weave strands relative to the spacing, margin is an inset
# labels_1/2/3 are vectors of identifiers
#
# Note: although the "diamond" option produces the same weave unit by default
# this function is more flexible, and preferred...
getHexPrimitiveCell <- function(spacing = 500, aspect = 1, margin = 0,
                                labels_1 = letters[1:2], labels_2 = letters[3:4],
                                labels_3 = letters[5:6], crs = 3857) {
  # aspect must be <= 1 / 2.S3
  width <- min(aspect, 1 / 2 / S3) * spacing
  length <- spacing - width * 2 / S3

  # the lengths of the label vectors
  ns <- c(length(labels_1), length(labels_2), length(labels_3))
  angles <- 0:2 * 120
  polys <- list()
  i <- 1
  for (dir in 1:3) {
    rot <- get_rot_matrix(angles[dir]) # required rotation matrix
    n <- ns[dir] # the number of lengthwise slices of the strand in this direction
    for (p in 1:n) {
      # parallelogram in the bottom left quadrant
      p1 <- get_parallelogram(c(width / S3, -width) * (p - 1) / n,
                              length, width / n)
      # and slid horizontally along by the spacing
      p2 <- p1 + c(spacing, 0)
      # add them to the list rotated by the required amount
      polys[[i]] <- p1 * rot
      polys[[i + 1]] <- p2 * rot
      i <- i + 2
    }
  }
  # make up a list of the labels
  ids <- rep(c(labels_1, labels_2, labels_3), 2) %>%
    matrix(sum(ns), 2) %>% t() %>% c()
  # and a hexagonal tile
  tile <- get_hexagon(spacing)
  list(
    cell = st_sf(strand = ids, geometry = st_as_sfc(polys)) %>%
      st_buffer(-margin) %>%          # inset margin
      group_by(strand) %>%                # dissolve by id
      summarise() %>%
      st_intersection(tile) %>%       # hex
      st_as_sf() %>%
      st_set_crs(crs),
    transform = wk::wk_affine_identity(), # no transform needed
    tile = tile,
    strands = unique(ids)
  )
}

# Returns a 'cube tile like this
#    ____
#   /   /\
#  /___/  \
#  \   \  /
#   \___\/   but rotated 90 degrees!
#
# and which doesn't really read as a weave!
getCubePrimitiveCell <- function(spacing = 500, margin = 0,
                                 labels_1 = letters[1:2], labels_2 = letters[3:4],
                                 labels_3 = letters[5:6], crs = 3857) {
  L <- spacing / S3
  W <- spacing / 2

  ns <- c(length(labels_1), length(labels_2), length(labels_3))
  angles <- c(30, 150, 270)
  polys <- list()
  i <- 1
  for (dir in 1:3) {
    rot <- get_rot_matrix(angles[dir])
    n <- ns[dir]
    for (p in 1:n) {
      p1 <- get_parallelogram(c(-W / S3, W) * (p - 1) / n, L, W / n, "UR")
      polys[[i]] <- p1 * rot
      i <- i + 1
    }
  }
  ids <- c(labels_1, labels_2, labels_3)
  tile <- get_hexagon(spacing)
  return(list(
    cell = st_sf(strand = ids, geometry = st_as_sfc(polys)) %>%
      st_intersection(tile) %>%
      st_buffer(-margin) %>%
      st_as_sf() %>%
      st_set_crs(crs),
    transform = wk::wk_affine_identity(),
    tile = tile,
    strands = unique(ids)
  ))
}


get_parallelogram <- function(origin, L, W, quadrant = "BL") {
  top_edge <- c(-L, 0)
  left_edge <- c(W / S3, -W)
  bottom_edge <- -top_edge
  right_edge <- -left_edge
  if (quadrant == "BL") {
    return(get_polygon(
      c(origin, origin + top_edge,
        origin + top_edge + left_edge,
        origin + top_edge + left_edge + bottom_edge)
    ))
  } else {
    return(get_polygon(
      c(origin, origin + bottom_edge,
        origin + bottom_edge + right_edge,
        origin + bottom_edge + right_edge + top_edge)
    ))
  }
}


# returns hexagon of the specified width (face to face normal distance)
# if point_up then point up orientation, else points left/right
get_hexagon <- function(w, point_up = TRUE) {
  r <- w / S3
  angles <- seq(1, 12, 2) / 6 * pi
  if (!point_up) {
    angles <- angles - pi / 6
  }
  get_polygon(c(matrix(c(cos(angles), sin(angles)), ncol = 6, byrow = TRUE)) * r)
}



## NOTE because the diamond tile makes a triangular weave identical to the
## "hex" type, it's not deprecated exactly, but it's not really needed
## It's why the 'transform' was added so historically important, but not
## no longer required (a useful mutation...)

# Returns a diamond weave primitive cell like
#
#     /\
#    /_/\
#   /__\ \
#   \/  \/
#    \__/
#     \/
#
# spacing is the total width of this diamond unit tile. If spacing is sqrt(3)
# times the width it will be increased to this. Between ...
# width is the normal width of the parallelograms (i.e. their height)
# margin is an inset requested to create gaps between the elements.
#
# labels_1, 2, 3 are lists of (character) labels distinguishing the
# three directions -- general "a", "b", "c"
getDiamondPrimitiveCell <- function(spacing = 500, aspect = 1, margin = 0,
                                    labels_1 = letters[1], labels_2 = letters[2],
                                    labels_3 = letters[3], crs = 3857) {

  lbls <- c(labels_1, labels_2, labels_3)

  # width requested by aspect * spacing must be <= spacing / 2.S3
  W <- min(aspect, 1 / 2 / S3 ) * spacing
  L <- spacing - W * 2 / S3

  # make a parallelogram in upper-right quadrant of length L, 'height' W
  base_p <- get_parallelogram(c(0, 0), L, W, "UR")
  # replicate 3 times, rotated around 0, 0
  base_polys <- list(base_p,                        # the base polygon
                     base_p * get_rot_matrix(120),  # rotated 120
                     base_p * get_rot_matrix(240))  # rotated 240
  # 4 translations to the corners of the tile unit
  translations <- c(0, 0,                           # in place
                    spacing / 2, -spacing * S3 / 2, # down to the right
                    spacing, 0,                     # across
                    spacing / 2, spacing * S3 / 2)  # up to the right
  tile <- get_polygon(translations)                 # the diamond tile

  # make up polygons, from the base polygons translated in all 4 directions
  polys <- list(); item <- 1
  for (p in base_polys) {
    for (i in seq(1, 8, 2)) {
      polys[[item]] <- p + translations[i:(i + 1)]
      item <- item + 1
    }
  }
  # and ids in the same order
  ids <- c(rep(labels_1, 4), rep(labels_2, 4), rep(labels_3, 4))
  # now make an sf and further process as needed
  list(
    cell = st_sf(strand = ids, geometry = st_as_sfc(polys)) %>%
      st_buffer(-margin / 2) %>%  # inset margin
      st_intersection(tile) %>%   # intersect to the tile
      st_as_sf() %>%
      # because polygon edges lie on edges of the tile, we get
      # 0 area slivers and non POLYGON geoms, so remove them
      # this is why "hex" type may be preferable
      filter(st_geometry_type(.) == "POLYGON" & st_area(.) > 1e-10) %>%
      st_set_crs(crs),
    # the transform required to make this rectangular tile-able
    transform = wk::wk_affine_invert(
      affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)),
    tile = tile,
    strands = unique(ids)
  )
}
