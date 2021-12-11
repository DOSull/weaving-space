require(sf)

S3 <- sqrt(3)

get_weave_unit <- function(type = "plain", spacing = 100, aspect = 1, 
                           margin = 0, n = c(2, 2), strands = "a|b|c",
                           tie_up = this_tu, tr = diag(nrow(tie_up)),
                           th = diag(ncol(tie_up)), crs = 3857) {
  aspect_messages(aspect)
  margin_messages(margin, max_margin = (1 - aspect) / 2)

  if (type %in% c("hex", "cube")) {
    unit <- get_triaxial_weave_unit(spacing = spacing, aspect = aspect, 
                                    margin = margin, type = type,
                                    strands = strands, crs = crs)
  } else {
    unit <- get_biaxial_weave_unit(spacing = spacing, aspect = aspect,
                                   margin = margin, type = type, n = n, 
                                   strands = strands, crs = crs,
                                   tie_up = tie_up, tr = tr, th = th)
  }
  # return(unit)
  list(
    weave_unit = unit$weave_unit,
    transform = wk::wk_affine_identity(),
    strands = unique(unit$weave_unit$strand),
    tile = unit$tile,
    type = type
  )
}


# convenience function to make a sfc POLYGON from a vector
# of points as c(x0,y0,x1,y1,...xn,yn) - note not closed
# this function will close it
get_polygon <- function(pts) {
  mpts <- matrix(c(pts, pts[1:2]), ncol = 2, byrow = TRUE)
  st_polygon(list(mpts))
}

# base rectangles with length, width, and specified orientation
# centred on the origin, and sliced lengthwise into n_split elements
get_base_rect <- function(L, W, orientation = "horizontal", n_split = 1) {
  # coords of a unit square centred on origin
  base <- 0.5 * matrix(c(-1, -1, 1, -1, 1, 1, -1, 1), nrow = 2)
  rects <- list()
  offsets <- seq(1, 2 * n_split - 1, 2) / n_split / 2 * W - W / 2
  for (i in 1:n_split) {
    if (orientation == "horizontal") {
      rects[[i]] <- get_polygon(c(matrix(c(L, 0, 0, W / n_split),
                                         nrow = 2) %*% base)) + c(0, offsets[i])
    } else {
      rects[[i]] <- get_polygon(c(matrix(c(W / n_split, 0, 0, L),
                                         nrow = 2) %*% base)) + c(offsets[i], 0)
    }
  }
  rects %>% st_as_sfc()
}


# translates and returns the supplied sf shapes by the offsets dx and dy
sf_translate <- function(shapes, dx = 0, dy = 0) {
  shapes %>% sf_transform(wk::wk_affine_translate(dx, dy))
}

# rotates the supplied sf shapes by the specified angle in degrees, around
# the specified centre coordinates
sf_rotate <- function(shapes, angle, cx = 0, cy = 0) {
  shapes %>% sf_transform(affine_rotn_around_xy(angle, cx, cy))
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
  shapes %>% sf_transform(
    wk::wk_affine_invert(affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2)))
}

# affine transforms the supplied sf shapes such that the unit square
# becomes the diamond shape (see previous)
sf_square_to_diamond <- function(shapes) {
  shapes %>% sf_transform(affine_abcd(0.5, -S3 / 2, 0.5, S3 / 2))
}

# returns coordinates of the centroid of the bounding box of supplied sf shapes
sf_get_centroid <- function(shapes) {
  shapes %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_centroid() %>%
    st_coordinates()
}

## wrapper functions for wk transforms applied to sf geometries
geoms_transform <- function(geoms, transform) {
  geoms %>%
    lapply(wk::wk_collection) %>%
    lapply(wk::wk_transform, trans = transform) %>%
    sapply(st_as_sfc, simplify = TRUE) %>%
    st_as_sfc() %>%
    st_cast() %>%      ### I don't know why, but this is VITAL
    st_set_crs(st_crs(geoms))
}

# wrapper function to transform the suppied sf by the provide wk transform
sf_transform <- function(shapes, transform) {
  st_geometry(shapes) <- geoms_transform(st_geometry(shapes), transform)
  shapes
}

# returns the wk affine transform for the rotation by angle (in degrees)
# around centre of rotation cx, cy
affine_rotn_around_xy <- function(angle, cx, cy) {
  wk::wk_affine_compose(
    wk::wk_affine_translate(-cx, -cy),
    wk::wk_affine_rotate(angle),
    wk::wk_affine_translate(cx, cy)
  )
}

# returns the wk affine transform that will transform the unit square
# basis vectors 0,1 and 1,0 to a,b and c,d, while preserving area
affine_abcd <- function(a, b, c, d) {
  area_scale <- abs(a * d - b * c)
  wk::wk_affine_compose(
    wk::wk_trans_affine(matrix(c(a, b, 1, c, d, 1, 0, 0, 1), 3, 3)),
    wk::wk_affine_scale(1 / sqrt(area_scale), 1 / sqrt(area_scale))
  )
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
  c(labels_1, labels_2, labels_3)
}

# converts a string to a vector of characters
string_to_chars <- function(s) {
  stringr::str_sub(s, 1:stringr::str_length(s), 1:stringr::str_length(s))
}

parse_strand_label <- function(s) {
  clean_s <- s %>%
    stringr::str_replace("[(]+", "(") %>%
    stringr::str_replace("[)]+", ")")
  result <- c()
  combo <- FALSE
  current <- ""
  for (i in 1:stringr::str_length(clean_s)) {
    next_char <- stringr::str_sub(clean_s, i, i)
    if (combo) {
      if (next_char == ")") {
        result <- c(result, current)
        current <- ""
        combo <- FALSE
      } else {
        current <- stringr::str_c(current, next_char, sep = "")
      }
    } else {
      if (next_char == "(") {
        combo <- TRUE
      } else {
        result <- c(result, next_char)
      }
    }
  }
  result
}

get_strand_ids <- function(strands_spec) {
  strands_spec %>%              # e.g. "a(bc)|ef-"
    parse_labels() %>%          # c("a(bc)", "ef-", "-")
    lapply(parse_strand_label)  # list(c("a", "bc"), c("e", "f", "-"), c("-"))
}


plot_unit <- function(unit, bg = "white") {
  if (any(unit$weave_unit$strand == "NA")) {
    main <- "Red areas are unresolved"
  } else {
    main <- "Different colours can be split on the strand variable"
  }
  unit$tile %>%
    plot(border = "red", col = rgb(0, 0, 0, 0,5), lwd = 3, lty = 2, 
         reset = FALSE, main = main)
  unit$weave_unit %>% 
    filter(strand != "NA") %>% 
    plot(add = TRUE, lwd = 0.5, bg = bg)
  unit$weave_unit %>% 
    filter(strand == "NA") %>%
    plot(add = TRUE, border = "black", col = rgb(1, 0, 0, 0.75))
  invisible(unit)
}

margin_messages <- function(margin, max_margin) {
  if (margin > max_margin) {
    warning(strwrap(
      paste("The largest margin that won't show gaps at tile boundaries with",
        " this spacing and aspect is around ", signif(max_margin, 3),
        ". Instead, consider setting margin to that value or less, or set it",
        " to 0, and apply a negative buffer to the final woven",
        " map elements.", sep = ""), prefix = " ", initial = ""))
  } 
  if (margin < 0) {
    warning(strwrap(
      paste("Negative margin settings are unlikely to produce patterns that",
            " read as weaves. They also probably contain overlapping shapes.",
            " And they may even cause topology exceptions, which is a fail.",
            " The tiles look cool plotted with tmap and alpha values < 1!",
            " Knock yourself out, it's all good.",
            sep = ""), prefix = " ", initial = ""))
  } 
  invisible(margin)
}

aspect_messages <- function(aspect) {
  if (aspect == 0) {
    stop(strwrap(
      paste("We have to divide things by the aspect setting sometimes, so",
            " setting aspect to 0 is not allowed. Feel free to try values",
            " outside the range 0 to 1. They won't look like weaves, but they",
            " might be fun anyway.", sep = ""), prefix = " ", initial = ""))
  } 
  if (aspect < 0 || aspect > 1) {
    warning(strwrap(
      paste("Here be dragons! Setting aspect to negative values or to values",
            " greater than 1 is unlikely to produce tiles that yield woven",
            " patterns. At extreme values you'll probably get topology errors,",
            " but between -1 and 2 you get interesting effects, especially if",
            " combined with negative margin settings. Not a bad idea to save",
            " your work as unexpected crashes are also a possibility :-)",
            sep = ""), prefix = " ", initial = ""))
  }
  invisible(aspect)
}



############ TRY THESE!
# the below options make the pattern give or take an offset in the colours)
# from Figure 22 of Glassner 2002
this_tu <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1,
                    0, 0, 0, 1, 0, 1, 1, 1,
                    0, 0, 1, 0, 1, 0, 1, 1,
                    0, 1, 0, 1, 0, 1, 0, 1,
                    1, 0, 1, 0, 1, 0, 1, 0,
                    1, 1, 0, 1, 0, 1, 0, 0,
                    1, 1, 1, 0, 1, 0, 0, 0,
                    1, 1, 1, 1, 0, 0, 0, 0), 8, 8, byrow = TRUE)
ids <- "aaaaaaaabbbbbbbb|aaaaaaaacccccccc"
