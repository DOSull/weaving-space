require(pracma)
require(abind)
require(tidyr)
require(dplyr)
require(sf)
require(wk)

# Functions that can be used to generate sf data 'weave units' i.e. a 
# tileable repeating element that when tiled gives the appearance of a 
# biaxial woven surface composed of criss-crossing rectangular
# elements. Implementation is based on ideas discussed in variously
#
# Glassner, A. 2002. Digital weaving. 1. IEEE Computer Graphics and 
#    Applications 22 (6):108–118.
# ———. 2003a. Digital weaving. 3. IEEE Computer Graphics and 
#    Applications 23 (2):80–83.
# ———. 2003b. Digital weaving. 2. IEEE Computer Graphics and 
#    Applications 23 (1):77–90.
#
# and (unpublished)
# 
# Griswold, R. 2006. Mathematical and Computational Topics in Weaving. 
# www2.cs.arizona.edu/patterns/weaving/webdocs/mo/Griswold-MO.pdf 
# (last accessed 29 October 2021).
#
# where weaving is shown to be essentially a matrix multiplication of 
# tie-up, threading and treadling matrices. An accessible introduction 
# can be found at
#
# https://www.youtube.com/watch?v=oMOSiag3dxg
#


# returns a matrix giving where 1 indicates warp on top, 2
# indicates weft on top
# types: 
# "plain"    : over 1 under 1 in both directions
# "twill"    : over n under n each weft thread, shifting 
#              one along between rows
# n          : the over-under for twill patterns
# warp, weft : a vector of distinct values (ints or chars) where each 
#              indicates a different thread colour; repeats are allowed, 
#              and "-" indicates that a thread should be skipped
get_weave_pattern_matrix <- function (type = "plain", n = 2, 
                                      warp = letters[1:2], 
                                      weft = letters[3:4], 
                                      tie_up = this_tu,
                                      th = this_th, tr = this_tr) {
  # make a sequence of ints from the warp/weft, 
  # substituting -1 for any "-"
  idxs <- 1:(length(c(warp, weft)))
  missing <- c(which(warp == "-"), length(warp) + which(weft == "-"))
  idxs[missing] <- -1

  width <- length(warp)
  height <- length(weft)
  p <- switch(
    type,
    "random" = make_random_pattern(warp = width, weft = height),

    "plain" = make_plain_pattern(warp_n = width, weft_n = height),

    "twill" = make_twill_pattern(n = n, 
                                 warp_n = width, weft_n = height),
    
    "basket" = make_basket_pattern(n = n, 
                                   warp_n = width, weft_n = height),
    
    "this" = make_this_pattern(tie_up = tie_up, th = th, tr = tr, 
                               warp_n = width, weft_n = height)
  )
  nc <- ncol(p)
  nr <- nrow(p)
  # make matrices of columns/rows for the warp and weft threads
  warp_threads <- matrix(head(idxs, width), nr, nc, byrow = TRUE)
  weft_threads <- matrix(tail(idxs, height), nr, nc)
  # handle missing threads (which will change the warp/weft assignment)
  return(p %>% 
           modify_pattern_for_missing_threads(
             warp_threads, weft_threads) %>% augment())
}

# add a row and a column to a matrix by copying the first to the last
augment <- function(mat, x = 1) {
  rc <- dim(mat)
  result <- repmat(mat, 2)
  return(result[1:(rc[1] + x), 1:(rc[2] + x)])
}

# returns number of repetitions of vector 
# length m needed to match some number of 
# repetitions of vector of length n
reps_needed <- function(n, m) {
  return(Lcm(n, m) / m)
}

# given tie up treadling and threading matrix returns 
# a 1/2 pattern matrix
get_pattern <- function(tie_up, treadling, threading, warp_n, weft_n) {
  rep_weft <- reps_needed(weft_n, nrow(tie_up)) # repeat to match weft
  rep_warp <- reps_needed(warp_n, ncol(tie_up)) # repeat to match warp
  return(
    (((repmat(treadling, n = rep_weft, m = 1)) %*%
        tie_up %*%     
        (repmat(threading, n = 1, m = rep_warp))) > 0) + 1)
}



# given warp weft pattern and column and row matrices, handles 
# missing threads
# pattern : matrix of 1 = warp on top 2 = weft on top 
# warp    : column matrix of ints where -1 is a missing thread
# weft    : row matrix of ints where -1 is a missing thread
modify_pattern_for_missing_threads <- function(pattern, warp, weft) {
  # threads missing in a direction mean the other direction 
  # should be on top
  pattern[which(warp < 0)] <- 2
  pattern[which(weft < 0)] <- 1
  return(pattern)
}

# NO LONGER NEEDED - but useful for debugging - pass output from 
# get_weave_pattern_matrix to get a matrix that can be imaged, i.e.
#
# get_weave_pattern_matrix %>% assign_indexes %>% image()
#
# returns an int matrix picking from the warp and weft based on the 
# values in the pattern matrix
# pattern : matrix of 1/2 encoding warp/weft on top
# warp    : matrix of int columns encoding distinct colours
# weft    : matrix of int rows encoding distinct colours
assign_indexes <- function(pattern, warp, weft) {
  # stack the warp and weft matrices into a 2 layer array
  loom <- abind(list(warp, weft), along = 3)
  nr <- nrow(pattern)
  nc <- ncol(pattern)
  ij <- expand.grid(1:nr, 1:nc)
  # use values in the pattern to select from either layer 1 or layer 2
  return(loom[cbind(ij[, 1], ij[, 2], c(pattern))] %>%
           matrix(nr, nc))
}


# simple over-under weave
make_plain_pattern <- function(warp_n = 1, weft_n = 1) {
  return(make_twill_pattern(n = 1, warp_n = warp_n, weft_n = weft_n))
}

# twill weave with n the number of over-unders
# note this is used with n = 1 to make plain weaves
make_twill_pattern <- function(n = 2, warp_n = 2, weft_n = 2) {
  ou <- n
  if (length(ou) == 1) {
    ou <- rep(n, 2)
  }
  dimension <- sum(ou)
  tie_up <- make_twill_matrix(ou, dimension)
  threading <- diag(nrow(tie_up)) 
  treadling <- diag(ncol(tie_up)) 
  return(get_pattern(tie_up, treadling, threading, warp_n, weft_n))
}


# returns a vector of runs of 1s and 0s 
# per the supplied vector. If n is a single
# value it is converted to c(n, n)
make_over_under_row <- function(n) {
  ou <- n
  if (length(n) == 1) {
    ou <- rep(n, 2)
  }
  x <- 1
  row <- c()
  for (y in ou) {
    row <- c(row, rep(x, y))
    x <- 1 - x
  }
  return(row)
}

# wraps a vector
# by : the number of positions to shift the row
# r  : the row
wrap_row <- function(by, r) {
  return(c(tail(r, by), head(r, length(r) - by)))
}

# makes a matrix like
# 1 1 0 0
# 0 1 1 0
# 0 0 1 1
# 1 0 0 1
# where the repeat runs in each row are length n
make_twill_matrix <- function(over_under, d) {
  row <- make_over_under_row(over_under)
  out <- row
  for (by in 2:d) {
    row <- wrap_row(1, row)
    out <- c(out, row)
  }
  return(matrix(out, d, d, byrow = TRUE))
}


make_basket_pattern <- function(n = 2, warp_n = 2, weft_n = 2) {
  tie_up <- make_basket_matrix(n)
  threading <- diag(nrow(tie_up)) 
  treadling <- diag(ncol(tie_up)) 
  return(get_pattern(tie_up, treadling, threading, warp_n, weft_n))
}

# makes a matrix like
# 1 1 0 0
# 1 1 0 0
# 0 0 1 1
# 0 0 1 1
# where the repeat runs in each row are length n
make_basket_matrix <- function(n) {
  return(
    matrix(
      c(rep(make_over_under_row(n), n), 
        rep(rev(make_over_under_row(n)), n)), n * 2, n * 2))
}


# stuff it let's see what happens!
make_this_pattern <- function(tie_up = this_tu, 
                              th = this_th, tr = this_tr,
                              warp_n = 2, weft_n = 2) {
  rep_warp <- reps_needed(weft_n, nrow(tie_up))
  threading <- th %>% 
    repmat(n = rep_warp, m = 1) 
  
  rep_weft <- reps_needed(warp_n, ncol(tie_up))
  treadling <- tr %>%
    repmat(n = 1, m = rep_weft)
  
  return(get_pattern(tie_up, treadling, threading, warp_n, weft_n))
}

# This function makes a random pattern (as a matrix of values) 
# with the number of different warp and weft threads specfied by
# warp and weft. Default values will make a 2 x 2 repeating unit. 
make_random_pattern <- function(n = 4, warp_n = 1, weft_n = 1) {
  width <- Lcm(n, warp_n)
  height <- Lcm(n, weft_n)
  tie_up <- matrix(sample(0:1, width * height, replace = TRUE), 
                   height, width)
  treadling <- make_matrix_from_seq(sample(1:height, width))
  threading <- make_matrix_from_seq(sample(1:width, height))
  return(
    get_pattern(tie_up, treadling, threading, warp_n, weft_n)
  )
}

zeros_with_a_one <- function(idx, n) {
  z <- rep(0, n)
  z[idx] <- 1
  return(z)
}

make_matrix_from_seq <- function(row_picks) {
  nrows <- max(row_picks)
  ncols <- length(row_picks)
  values <- sapply(row_picks, zeros_with_a_one, n = nrows)
  return(matrix(values, nrows, ncols))
}

# translate a polygon by the supplied (x, y). 
# The ordering of parameters facilitates use in an `lapply`.
translate_poly <- function(pt, poly) {
  return(poly + pt)
}

# makes polygons
make_polys <- function(L, W, wow, dx, dy) {
  orientations <- c("vertical", "horizontal")
  polys <- list()
  polys[[1]] <- get_base_rect(L, W, orientations[wow]) + c(dx, dy)
  if (L == W) {
    return(polys)
  }
  gap <- (L - W) / 2
  ddx <- ifelse(wow == 2, 0, (W + gap) / 2)
  ddy <- ifelse(wow == 2, (W + gap) / 2, 0)
  orientation <- orientations[3 - wow]
  polys[[2]] <- get_base_rect(gap, W, orientation) + 
    c(dx - ddx, dy - ddy)
  polys[[3]] <- get_base_rect(gap, W, orientation) + 
    c(dx + ddx, dy + ddy)
  return(polys)
}

make_polygons_from_matrix <- function(ww, spacing, aspect, margin,
                                        warp, weft, crs) {
  h <- nrow(ww)
  w <- ncol(ww)
  L <- 2 * spacing / (1 + aspect)
  W <- L * aspect
  bb <- c(xmin = 0, ymin = 0, xmax = spacing * (w - 1), 
          ymax = spacing * (h - 1))
  weft_ids <- rep(weft, ceiling(h / length(weft)))
  warp_ids <- rep(warp, ceiling(w / length(warp)))
  rc <- expand_grid(row = 1:h, col = 1:w)
  polys <- list()
  out_ids <- c()
  n <- 1
  for (k in 1:(dim(rc)[1])) { #c(mats$ids))) {
    row <- rc$row[k]
    col <- rc$col[k]
    next_polys <- make_polys(L, W, ww[row, col], 
                              spacing * (col - 1), spacing * (row - 1))
    polys[[n]] <- next_polys[[1]]
    out_ids <- append(out_ids, ifelse(ww[row, col] == 1, 
                                      warp_ids[col], weft_ids[row]))
    n <- n + 1
    if (length(next_polys) > 1) {
      polys[[n]] <- next_polys[[2]]
      polys[[n + 1]] <- next_polys[[3]]
      out_ids <- append(out_ids, 
                        rep(ifelse(ww[row, col] == 1, 
                                   weft_ids[row], warp_ids[col]), 2))
      n <- n + 2
    }
  }
  polys <- polys %>%
    st_as_sfc() %>% st_sf() %>%
    # add in the id attribute and dissolve
    mutate(id = out_ids) %>%    # the indices into the thread names
    filter(id != "-") %>%       # throw away the missing ones coded -1
    group_by(id) %>%            # dissolve
    summarise() %>%
    st_buffer(-margin) %>%
    st_crop(bb) %>%
    st_set_crs(crs)
  return(polys)
}



# hard-coded defaults for the arbitrary weave
this_tr <- matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0, 
                    0, 0, 1, 0,
                    0, 0, 0, 1,
                    0, 0, 1, 0,
                    0, 1, 0, 0), 6, 4, byrow = TRUE)
this_tu <- matrix(c(0, 0, 1, 
                    0, 1, 0,
                    1, 1, 0, 
                    1, 0, 0), 4, 3, byrow = TRUE)
this_th <- matrix(c(0, 0, 1,
                    0, 1, 0, 
                    1, 0, 0), 3, 3, byrow = TRUE)

get_biaxial_weave_unit <- function(spacing = 10000, aspect = 1, 
                                   margin = 0, type = "plain", 
                                   n = c(2, 2), # used by twill
                                   ids = "ab|cd", crs = 3857,
                                   tie_up = this_tu, 
                                   tr = this_tr, th = this_th) {
  
  parsed_labels = ids %>% parse_labels() %>% lapply(string_to_chars)
  warp_threads = parsed_labels[[1]]
  weft_threads = parsed_labels[[2]]
  
  if (type == "basket") {
    n = n[1]
  }
  cell <- get_weave_pattern_matrix(type = type, n = n, 
                                   warp_threads, weft_threads, 
                                   tie_up = tie_up, 
                                   tr = tr, th = th) %>%
    make_polygons_from_matrix(
      spacing = spacing, margin = margin, aspect = aspect,
      warp_threads, weft_threads, crs = crs)
  return(
    list(
      primitive = cell,
      transform = wk_affine_identity(),
      ids = unique(cell$id),
      type = type
    )
  )
}
