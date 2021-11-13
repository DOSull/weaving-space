require(dplyr)
require(sf)
require(wk)

S3 <- sqrt(3)

# convenience function to make a sfc POLYGON from a vector
# of points as c(x0,y0,x1,y1,...xn,yn) - note not closed
# this function will close it
get_polygon <- function(pts) {
  mpts <- matrix(c(pts, pts[1:2]), ncol = 2, byrow = TRUE)
  return(st_polygon(list(mpts)))
}

# base rectangle with length, width, and specified orientation
# centred on the origin
get_base_rect <- function(L, W, orientation = "horizontal") {
  # coords of a unit square centred on origin
  base <- 0.5 * matrix(c(-1, -1, 1, -1, 1, 1, -1, 1), nrow = 2)
  if (orientation == "horizontal") {
    return(get_polygon(c(matrix(c(L, 0, 0, W), nrow = 2) %*% base)))
  } else {
    return(get_polygon(c(matrix(c(W, 0, 0, L), nrow = 2) %*% base)))
  }
}

# translates and returns the supplied sf shapes by the offsets dx and dy
sf_translate <- function(shapes, dx = 0, dy = 0) {
  return(shapes %>% sf_transform(wk_affine_translate(dx, dy)))
}

# rotates the supplied sf shapes by the specified angle in degrees, around
# the specified centre coordinates
sf_rotate <- function(shapes, angle, cx = 0, cy = 0) {
  return(shapes %>% sf_transform(affine_rotn_around_xy(angle, cx, cy)))
}

# affine transforms the supplied sf shapes such that the diamond
#
#      0.5,S3/2
#      /\
#  0,0/  \1,0
#     \  /
#      \/
#      0.5,S3/2
#
# is transformed to the unit square 0,0 1,0 1,1 0,1
sf_diamond_to_square <- function(shapes) {
  return(shapes %>% sf_transform(
    wk_affine_invert(
      affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)
    )
  ))
}

# affine transforms the supplied sf shapes such that the unit square
# becomes the diamond shape (see previous)
sf_square_to_diamond <- function(shapes) {
  return(shapes %>% sf_transform(
    affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)))
}

# returns coordinates of the centroid of the bounding box of the supplied sf shapes
sf_get_centroid <- function(shapes) {
  return(shapes %>%
           st_bbox() %>%
           st_as_sfc() %>%
           st_centroid() %>%
           st_coordinates())
}

## wrapper functions for wk transforms applied to sf geometries
geoms_transform <- function(geoms, transform) {
  return(
    geoms %>%
      lapply(wk_collection) %>%
      lapply(wk_transform, trans = transform) %>%
      sapply(st_as_sfc, simplify = TRUE) %>%
      st_as_sfc() %>%
      st_cast() %>%      ### I don't know why, but this is VITAL
      st_set_crs(st_crs(geoms))
  )
}

# wrapper function to transform the suppied sf by the provide wk transform
sf_transform <- function(shapes, transform) {
  st_geometry(shapes) <- geoms_transform(st_geometry(shapes), transform)
  return(shapes)
}

# returns the wk affine transform for the rotation by angle (in degrees)
# around centre of rotation cx, cy
affine_rotn_around_xy <- function(angle, cx, cy) {
  return(wk_affine_compose(
    wk_affine_translate(-cx, -cy),
    wk_affine_rotate(angle),
    wk_affine_translate(cx, cy)
  ))
}

# returns the wk affine transform that will transform the unit square
# basis vectors 0,1 and 1,0 to a,b and c,d, while preserving area
affine_abcd <- function(a, b, c, d) {
  areaScale = abs(a * d - b * c)
  return(wk_affine_compose(
    wk_trans_affine(matrix(c(a, b, 1, c, d, 1, 0, 0, 1), 3, 3)),
    wk_affine_scale(1 / sqrt(areaScale), 1 / sqrt(areaScale))
  ))
}

# parses an ids string "ab|c|de" (for example) by splitting at '|'s
# i.e. to c("ab", "c", "de")
parse_labels <- function(ids) {
  lbls <- strsplit(ids, split = "|", fixed = TRUE)[[1]]
  labels_1 <- lbls[[1]][1] 
  if (length(lbls) > 1) {
    labels_2 <- lbls[2]
  } else {
    labels_2 <- "-"
  }
  if (length(lbls) > 2) {
    labels_3 <- lbls[3]
  } else {
    labels_3 <- "-"
  }
  return(c(labels_1, labels_2, labels_3))
}

# converts a string to a vector of characters
string_to_chars <- function(s) {
  return(substring(s, 1:nchar(s), 1:nchar(s)))
}

parse_strand_label <- function(s) {
  clean_s <- gsub("[(]+", "(", s)
  clean_s <- gsub("[)]+", ")", clean_s)
  result <- c()
  combo <- FALSE
  current <- ""
  for (i in 1:nchar(clean_s)) {
    nextChar <- substr(clean_s, i, i)
    if (combo) {
      if (nextChar == ")") {
        result <- c(result, current)
        current <- ""
        combo <- FALSE
      } else {
        current <- paste(current, nextChar, sep = "")
      }
    } else {
      if (nextChar == "(") {
        combo <- TRUE
      } else {
        result <- c(result, nextChar)
      }
    }
  }
  return(result)
}
