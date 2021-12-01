require(sf)

# Functions to generate sf datasets that are tileable to produce the
# appearance of a woven surface, with three directions of weaving thread.


# convenience storage of sqrt(3)
S3 <- sqrt(3)


get_triaxial_weave_pattern_matrices <- function(type = "hex",
                                                strands_1 = c("a", "-", "-"),
                                                strands_2 = c("b", "-", "-"),
                                                strands_3 = c("c", "-", "-")) {

  if (type == "hex") {
    loom <- matrices_as_loom(
      get_weave_pattern_matrix(type = "this", tie_up = ones(6),
                               warp = strands_1, weft = strands_2),
      get_weave_pattern_matrix(type = "this", tie_up = ones(6),
                               warp = strands_2, weft = strands_3),
      get_weave_pattern_matrix(type = "this", tie_up = ones(6),
                               warp = strands_3, weft = strands_1))
  }
  if (type == "cube") {
    loom <- matrices_as_loom(
      get_weave_pattern_matrix(type = "twill", n = rep(1:2, 2),
                               warp = strands_1, weft = strands_2),
      get_weave_pattern_matrix(type = "twill", n = rep(1:2, 2),
                               warp = strands_2, weft = strands_3),
      get_weave_pattern_matrix(type = "twill", n = rep(1:2, 2),
                               warp = strands_3, weft = strands_1))
  }
  loom
}



get_triaxial_weave_unit <- function(spacing = 500, aspect = 1, margin = 0,
                                    strands = "a--|b--|c--", type = "hex", crs = 3857) {

  parsed_labels <- strands %>%  # e.g. "a(bc)|ef-"
    parse_labels() %>%          # c("a(bc)", "ef-", "-")
    lapply(parse_strand_label)  # list(c("a", "bc"), c("e", "f", "-"), c("-"))

  strands_1 <- parsed_labels[[1]]
  strands_2 <- parsed_labels[[2]]
  strands_3 <- parsed_labels[[3]]

  cell <- get_triaxial_weave_pattern_matrices(type = type,
                                              strands_1, strands_2, strands_3) %>%
    make_sf_from_coded_weave_matrix(spacing = spacing,
                                    width = aspect, margin = margin,
                                    axis1_threads = strands_1,
                                    axis2_threads = strands_2,
                                    axis3_threads = strands_3, crs = crs)

  list(
    primitive = cell$weave_unit,
    transform = wk::wk_affine_identity(),
    strands = unique(cell$weave_unit$strand),
    tile = cell$tile,
    type = type
  )
}


# function to add values to a supplied list tri, from a supplied matrix M, where
# values in the list will be indexed by triangular coordinates (converted to a
# string). This will be called 3 times, once for each axis. The calling context
# should determine the parity, based on the sizes of the three matrices
add_biaxial_to_triaxial <- function(tri, M, axis = 1, parity = 0) {
  for (col in 1:ncol(M)) {
    for (row in 1:nrow(M)) {
      # get the target grid coordinates (there will be two)
      abc <- transform_ab_to_abc(col, row, parity = parity, axis = axis)
      for (par in 1:2) {
        # make a key (string) from each pair
        key <- abc[par, ] %>% stringr::str_c(collapse = ",")
        # if it is not in the target list, add a new thread order
        if (is.null(tri[[key]])) {
          tri[[key]] <- list(decode_biaxial_to_order(M[row, col], axis = axis)) # c(M[row, col])
        } else { # if there is something there, then append the value
          tri[[key]] <- append(tri[[key]], list(decode_biaxial_to_order(M[row, col], axis = axis)))
          # tri[[key]] <- c(tri[[key]], M[row, col] )
        }
      }
    }
  }
  tri
}

# combines a set of orderings on the values
# the orderings are a list of vectors (which may be empty)
# for example list(c(1, 2), c(2, 3), c(1, 3))
combine_orderings <- function(..., values = 1:3, verbose = FALSE) {
  orderings <- list(...)
  # assemble a matrix of the positions of entries
  # in each ordering among the values
  ranks <- matrix(0, length(orderings), length(values))
  for (i in seq_along(orderings)) {
    ranks[i, ] <- match(values, orderings[[i]])
  }
  # replace any missing matches with a high score
  ranks[which(is.na(ranks))] <- 100
  # sum the ranks of each value
  scores <- colSums(ranks)
  max_score <- length(orderings) * 100
  number_present <- sum(scores < max_score)
  if (number_present == 0) {
    result <- NULL
  } else {
    result <- as.vector(values[order(scores)[1:number_present]])
  }
  if (verbose) {
    return(list(result = result, ranks = ranks, scores = scores))
  } else {
    return(result)
  }
}

# function to supply missing third coordinate in a triangular grid given two
# other coordinates, the parity, and the axis. The axis tells us which is missing:
# axis = 1 --> 3rd, axis = 2 --> 1st, axis = 3 --> 2nd
transform_ab_to_abc <- function(z1, z2, parity = 0, axis = 1) {
  # comments on first case, others are cyclic shifts - most easily read down the code
  # when axis = 1 a_coord and b_coord are preserved, c_coord is the paired residual values
  a_coord <- switch(axis, z1, c(parity + 1 - z1 - z2, parity - z1 - z2), z2)
  b_coord <- switch(axis, z2, z1, c(parity + 1 - z1 - z2, parity - z1 - z2))
  c_coord <- switch(axis, c(parity + 1 - z1 - z2, parity - z1 - z2), z2, z1)
  result <- switch(
    axis, # return two triples that differ in the missing coordinate
    matrix(c(a_coord, b_coord, c_coord[1],
             a_coord, b_coord, c_coord[2]), nrow = 2, ncol = 3, byrow = TRUE),
    matrix(c(a_coord[1], b_coord, c_coord,
             a_coord[2], b_coord, c_coord), nrow = 2, ncol = 3, byrow = TRUE),
    matrix(c(a_coord, b_coord[1], c_coord,
             a_coord, b_coord[2], c_coord), nrow = 2, ncol = 3, byrow = TRUE))
  result
}



# returns hexagon of the specified width (face to face normal distance)
# if point_up then point up orientation, else points left/right
get_hexagon <- function(w, point_up = TRUE) {
  r <- w / S3
  angles <- seq(1, 12, 2) / 6 * pi
  if (!point_up) {
    angles <- angles - pi / 6
  }
  get_polygon(c(matrix(c(cos(angles), sin(angles)), 
                       ncol = 6, byrow = TRUE)) * r)
}
