require(dplyr)
require(sf)

## ---- UTILITY FUNCTIONS ----

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
get_weave_pattern_matrix <- function(type = "plain", n = 2,
                                     warp = letters[1:2],
                                     weft = letters[3:4],
                                     tie_up = this_tu,
                                     th = diag(ncol(tie_up)),
                                     tr = diag(nrow(tie_up))) {
  # make a sequence of ints from the warp/weft,
  # substituting -1 for any "-"
  warps <- seq_along(warp)
  wefts <- seq_along(weft)
  warps[which(warp == "-")] <- -1
  wefts[which(weft == "-")] <- -1

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
  warp_threads <- matrix(warps, nr, nc, byrow = TRUE)
  weft_threads <- matrix(wefts, nr, nc)
  # encode to reflect missing threads
  p %>% encode_biaxial_weave(warp_threads, weft_threads)
}

# add a row and a column to a matrix by copying the first to the last
augment <- function(mat, x = 1) {
  rc <- dim(mat)
  result <- pracma::repmat(mat, 2)
  result[1:(rc[1] + x), 1:(rc[2] + x)]
}

augment_with_values <- function(mat, n = 1, values = 0) {
  rc <- dim(mat)
  result <- pracma::repmat(mat, 2)
  result[, (rc[2] + 1):(rc[2] + n)] <- values
  result[(rc[1] + 1):(rc[1] + n), ] <- values
  result[1:(rc[1] + n), 1:(rc[2] + n)]
}

# returns the lowest common multiple of v1 and v2 divided by each
# used to determine how many repetitions of v1 and v2 required to
# make matched length vectors
# Returns a vector with item 1 being the required repeats of v1
# and item 2 the required repeats of v2
reps_needed <- function(v1, v2) {
  n <- pracma::Lcm(v1, v2)
  c(n / v1, n / v2)
}

# Returns a 1/2 encoded weave matrix given tie_up, treadling and
# threading matrices. The following conditions must be satisfied to
# avoid non-conformable matrix error:
# The "twill", "random", "basket" and "plain" options should guarantee
# this, but the "this" option requires the user to make this happen
# If the warp_n and weft_n values are not factors of nrow(treadling) and
# ncol(threading) respectively, the output matrix will be repeated as
# needed to make this match
get_pattern <- function(tie_up, treadling, threading,
                        warp_n, weft_n, rep = 1) {
  pat <- ((treadling %*% tie_up %*% threading) > 0) + 0
  # now determine repetitions needed to match warp_n and weft_n
  rep_warp <- reps_needed(warp_n, ncol(pat))
  rep_weft <- reps_needed(weft_n, nrow(pat))
  pracma::repmat(pat, n = rep_weft[2] * rep, m = rep_warp[2] * rep)
}


# Note that as currently written this function requires the warp and weft
# matrices to be the same size, which get_weave_pattern_matrix will ensure, b
# but which may not be the case if called from elsewhere
encode_biaxial_weave <- function(pattern, warp, weft) {
  pattern[which(pattern == 1)] <- 5        # warp present and on top
  pattern[which(pattern == 0)] <- 4        # weft present and on top
  pattern[which(warp < 0)] <- 1            # warp absent
  pattern[which(weft < 0)] <- 2            # weft absent
  pattern[which(weft < 0 & warp < 0)] <- 3 # both absent
  pattern
}

decode_biaxial_to_order <- function(code, axis = 0) {
  if (axis == 0) { # biaxial
    return(switch(
      code, 1, 2, NULL, 1:2, 2:1
    ))
  }
  if (axis == 1) { # A&B triaxial
    return(switch(
      code, 2, 1, NULL, 2:1, 1:2
    ))
  }
  if (axis == 2) { # B&C triaxial
    return(switch(
      code, 3, 2, NULL, 3:2, 2:3
    ))
  }
  if (axis == 3) { # C&A triaxial
    return(switch(
      code, 1, 3, NULL, c(1, 3), c(3, 1)
    ))
  }
}

## ----- THE WEAVES -----

# simple over-under weave
make_plain_pattern <- function(warp_n = 1, weft_n = 1) {
  make_twill_pattern(n = 1, warp_n = warp_n, weft_n = weft_n)
}

# twill weave with n the number of over-unders
# note this is used with n = 1 to make plain weaves
make_twill_pattern <- function(n = 2, warp_n = 2, weft_n = 2) {
  ou <- n
  if (length(ou) == 1) {
    ou <- rep(n, 2)
  }
  tie_up <- make_twill_matrix(ou)
  threading <- diag(nrow(tie_up))
  treadling <- diag(ncol(tie_up))
  get_pattern(tie_up, treadling, threading, warp_n, weft_n)
}

# returns a vector of runs of 1s and 0s per the supplied vector.
# If the length of n is odd then it is doubled to produce an even
# length over-under sequence that repeats. If we don't do this then,
# e.g, 1:3 becomes 100111 which repeated is 1001111001111, i.e. a 2-4
# over-under pattern. Doubling it makes 100111011000 which has the
# requested pattern
make_over_under_row <- function(n) {
  ou <- n
  if (length(n) %% 2 != 0) {
    ou <- rep(n, 2)
  }
  x <- 1
  row <- c()
  for (y in ou) {
    row <- c(row, rep(x, y))
    x <- 1 - x
  }
  row
}

# wraps a vector
# by : the number of positions to shift the row
# r  : the row
wrap_row <- function(by, r) {
  c(tail(r, by), head(r, length(r) - by))
}

# makes a matrix like
# 1 1 0 0
# 0 1 1 0
# 0 0 1 1
# 1 0 0 1
# where the repeat runs in each row are length n
make_twill_matrix <- function(over_under) {
  row <- make_over_under_row(over_under)
  d <- length(row)
  out <- row
  for (i in 2:d) {
    row <- wrap_row(1, row)
    out <- c(out, row)
  }
  matrix(out, d, d, byrow = TRUE)
}

make_basket_pattern <- function(n = 2, warp_n = 2, weft_n = 2) {
  tie_up <- make_basket_matrix(n)
  threading <- diag(nrow(tie_up))
  treadling <- diag(ncol(tie_up))
  get_pattern(tie_up, treadling, threading, warp_n, weft_n)
}

# makes a matrix like
# 1 1 0 0
# 1 1 0 0
# 0 0 1 1
# 0 0 1 1
# where the repeat runs in each row are length n
make_basket_matrix <- function(n) {
  matrix(c(rep(make_over_under_row(n), n),
           rep(rev(make_over_under_row(n)), n)),
         n * 2, n * 2)
}

# This is just a pass through function. Could try to enforce
#   ncol(treadling) == nrow(tie_up) and ncol(tie_up) == nrow(threading)
# but unsure what would be an appropriate way to do this...
make_this_pattern <- function(tie_up = this_tu,
                              threading = diag(ncol(tie_up)),
                              treadling = diag(nrow(tie_up)),
                              warp_n = 2, weft_n = 2) {
  get_pattern(tie_up, treadling, threading, warp_n, weft_n)
}

# This function makes a random pattern (as a matrix of values)
# with the number of different warp and weft threads specfied by
# warp and weft. Default values will make a 2 x 2 repeating unit.
make_random_pattern <- function(n = 4, warp_n = 1, weft_n = 1) {
  width <- pracma::Lcm(n, warp_n)
  height <- pracma::Lcm(n, weft_n)
  tie_up <- matrix(sample(0:1, width * height, replace = TRUE),
                   height, width)
  treadling <- make_matrix_from_seq(sample(1:height, width))
  threading <- make_matrix_from_seq(sample(1:width, height))
  
  get_pattern(tie_up, treadling, threading, warp_n, weft_n)
}

zeros_with_a_one <- function(idx, n) {
  z <- rep(0, n)
  z[idx] <- 1
  z
}

make_matrix_from_seq <- function(row_picks) {
  nrows <- max(row_picks)
  ncols <- length(row_picks)
  values <- sapply(row_picks, zeros_with_a_one, n = nrows)
  matrix(values, nrows, ncols)
}



## ---- EXTERNAL API ----

get_biaxial_weave_unit <- function(spacing = 10000, aspect = 1, margin = 0,
                                   type = "plain", n = c(2, 2), # used by twill
                                   strands = "ab|cd", crs = 3857,
                                   tie_up = this_tu, tr = diag(nrow(tie_up)),
                                   th = diag(ncol(tie_up))) {
  margin_messages(margin, (1 - aspect) / 2)
  strand_ids <- get_strand_ids(strands)  
  warp_threads <- strand_ids[[1]]
  weft_threads <- strand_ids[[2]]

  if (type == "basket") {
    n <- n[1]
  }
  get_weave_pattern_matrix(type = type, n = n,
                           warp_threads, weft_threads,
                           tie_up = tie_up, tr = tr, th = th) %>%
    matrices_as_loom() %>%
    make_sf_from_coded_weave_matrix(spacing = spacing, width = aspect,
                                    margin = margin,
                                    axis1_threads = weft_threads,
                                    axis2_threads = warp_threads, crs = crs)
  # list(
  #   primitive = cell$weave_unit,
  #   transform = wk::wk_affine_identity(),
  #   strands = unique(cell$weave_unit$strand),
  #   tile = cell$tile,
  #   type = type
  # )
}
