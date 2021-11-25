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

# base rectangle with length, width, margin, and specified orientation
# centred on the origin
get_base_rect <- function(L, W, M, orientation = "horizontal") {
  # coords of a unit square centred on origin
  base <- 0.5 * matrix(c(-1, -1, 1, -1, 1, 1, -1, 1), nrow = 2)
  if (orientation == "horizontal") {
    return(get_polygon(c(matrix(c(L - M, 0, 0, W - M), nrow = 2) %*% base)))
  } else {
    return(get_polygon(c(matrix(c(W - M, 0, 0, L - M), nrow = 2) %*% base)))
  }
}

sf_translate <- function(shapes, dx = 0, dy = 0) {
  return(shapes %>% sf_transform(wk_affine_translate(dx, dy)))
}

sf_rotate <- function(shapes, angle, cx = 0, cy = 0) {
  return(shapes %>% sf_transform(affine_rotn_around_xy(angle, cx, cy)))
}

sf_diamond_to_square <- function(shapes) {
  return(shapes %>% sf_transform(
    wk_affine_invert(
      affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)
    )
  ))
}

sf_square_to_diamond <- function(shapes) {
  return(shapes %>% sf_transform(
    affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)))
}

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

sf_transform <- function(shapes, transform) {
  st_geometry(shapes) <- geoms_transform(st_geometry(shapes), transform)
  return(shapes)
}

affine_rotn_around_xy <- function(angle, cx, cy) {
  return(wk_affine_compose(
    wk_affine_translate(-cx, -cy),
    wk_affine_rotate(angle),
    wk_affine_translate(cx, cy)
  ))
}

affine_abcd <- function(a, b, c, d) {
  areaScale = abs(a * d - b * c)
  return(wk_affine_compose(
    wk_trans_affine(matrix(c(a, b, 1, c, d, 1, 0, 0, 1), 3, 3)),
    wk_affine_scale(1 / sqrt(areaScale), 1 / sqrt(areaScale))
  ))
}

