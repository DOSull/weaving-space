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
                                      th = diag(ncol(tie_up)), 
                                      tr = diag(nrow(tie_up))) {
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

augment_with_values <- function(mat, n = 1, values = 0) {
  rc <- dim(mat) 
  result <- repmat(mat, 2)
  result[, (rc[2] + 1):(rc[2] + n)] <- values
  result[(rc[1] + 1):(rc[1] + n), ] <- values
  return(result[1:(rc[1] + n), 1:(rc[2] + n)])
}

# returns the lowest common multiple of v1 and v2 divided by each 
# used to determine how many repetitions of v1 and v2 required to 
# make matched length vectors
# Returns a vector with item 1 being the required repeats of v1 
# and item 2 the required repeats of v2
reps_needed <- function(v1, v2) {
  n <- Lcm(v1, v2)
  return(c(n / v1, n / v2))
}

# Returns a 1/2 encoded weave matrix given tie_up, treadling and 
# threading matrices. The following conditions must be satisfied to
# avoid non-conformable matrix error:
# ncol(treadling) == nrow(tie_up)
# ncol(tie_ip) == nrow(threading)
# The "twill", "random", "basket" and "plain" options should guarantee
# this, but the "this" option requires the user to make this happen
# If the warp_n and weft_n values are not factors of nrow(treadling) and 
# ncol(threading) respectively, the output matrix will be repeated as
# needed to make this match
get_pattern <- function(tie_up, treadling, threading, 
                        warp_n, weft_n, rep = 1) {
  pat <- ((treadling %*% tie_up %*% threading) > 0) + 1
  # now determine repetitions needed to match warp_n and weft_n
  rep_warp <- reps_needed(warp_n, ncol(pat))
  rep_weft <- reps_needed(weft_n, nrow(pat))
  return(repmat(pat, n = rep_weft[2] * rep, m = rep_warp[2] * rep))
}



# given warp weft pattern and column and row matrices, handles 
# missing threads
# pattern : matrix of 1 = warp on top 2 = weft on top 
# warp    : column matrix of ints where -1 is a missing thread
# weft    : row matrix of ints where -1 hhis a missing thread
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
  tie_up <- make_twill_matrix(ou)
  threading <- diag(nrow(tie_up)) 
  treadling <- diag(ncol(tie_up)) 
  return(get_pattern(tie_up, treadling, threading, warp_n, weft_n))
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
make_twill_matrix <- function(over_under) {
  row <- make_over_under_row(over_under)
  d <- length(row)
  out <- row
  for (i in 2:d) {
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


# this is just a pass through function
# could try to enforce ncol(treadling) == nrow(tie_up) and
# ncol(tie_up) == nrow(threading) 
# but unsure what would be an appropriate way to do this...
make_this_pattern <- function(tie_up = this_tu, 
                              threading = diag(ncol(tie_up)), 
                              treadling = diag(nrow(tie_up)),
                              warp_n = 2, weft_n = 2) {
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

# makes polygons at a particular spot in a weave unit (dx, dy)
# with length L, width W, orientation determined from
# wow == 1 --> horizontal, wow == 2 --> vertical, and with
# rectangles sliced lengthwise if required based on the supplied 
# warp or weft id (e.g. "ab" produces two narrower rectangles)
make_polys <- function(L, W, wow, dx, dy, warp_id, weft_id) {
  orientations <- c("vertical", "horizontal")
  n_slices <- c(str_length(warp_id), str_length(weft_id))
  over_polys <- get_base_rect(L, W, orientations[wow], 
                              n_slices[wow]) + c(dx, dy)
  if (L == W) { # no gaps for the cross (under) strand to show
    return(over_polys)
  }
  under_polys <- (get_base_rect(L, W, orientations[3 - wow], 
                                n_slices[3 - wow]) + c(dx, dy)) %>% 
    st_difference(st_union(over_polys))
  return(c(over_polys, under_polys))
}


make_polygons_from_matrix <- function(ww = matrix(c(1, 2, 1, 2, 1, 2, 1, 2, 1), 3, 3), 
                                      spacing = 10000, aspect = 0.6, margin = 0, 
                                      warp = letters[1:2], weft = letters[3:4],
                                      crs = 3857) {
  # set height, width, length and width of thread elements, and bounding box
  h <- nrow(ww)
  w <- ncol(ww)
  L <- spacing + spacing * (1 - aspect) #2 * spacing / (1 + aspect)
  W <- spacing * aspect
  bb <- c(xmin = 0, ymin = 0, xmax = spacing * (w - 1), ymax = spacing * (h - 1))
  # extend the lists of thread IDs in case we need to run over
  weft_ids <- rep(weft, ceiling(h / length(weft)))
  warp_ids <- rep(warp, ceiling(w / length(warp)))
  # empty list for the polygons and vector for the strand ids
  polys <- list()
  strand_ids <- c()
  for(row in 1:h) {
    for(col in 1:w) {
      # get the next set of polygons
      next_polys <- make_polys(spacing, W, ww[row, col],
                               spacing * (col - 1), spacing * (row - 1),
                               warp_ids[col], weft_ids[row])
      # get number of ids in strand on top
      thread_ids <- c(warp_ids[col], weft_ids[row])
      n_on_top <- str_length(thread_ids)[ww[row, col]]
      for (i in seq_along(next_polys)) {
        # add to the list of polygons
        polys <- append(polys, list(next_polys[[i]]))
        # strand id is from the spec on top, or not
        id <- ifelse(i <= n_on_top, thread_ids[ww[row, col]] %>% substr(i, i),
                     thread_ids[3 - ww[row, col]] %>% substr(i - n_on_top, i - n_on_top))
        strand_ids <- c(strand_ids, id)
      }
    }
  }
  polys <- polys %>%
    st_as_sfc() %>% st_sf() %>%     # make into an sf
    mutate(strand = strand_ids) %>% # the indices into the thread names
    filter(strand != "-") %>%       # throw away the missing ones coded -1
    st_crop(bb) %>%                 # crop to bounding box and remove any slivers
    filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) %>%
    group_by(strand) %>%            # dissolve on the strand
    dplyr::summarise() %>%
    st_buffer(-margin) %>%          # do the margin inset
    st_set_crs(crs)                 # set the CRS    
  
  return(list(weave_unit = polys,
              tile = bb))
}


get_biaxial_weave_unit <- function(spacing = 10000, aspect = 1, margin = 0, 
                                   type = "plain", n = c(2, 2), # used by twill
                                   strands = "ab|cd", crs = 3857, tie_up = this_tu, 
                                   tr = diag(nrow(tie_up)), th = diag(ncol(tie_up))) {
  
  parsed_labels <- strands %>%  # e.g. "a(bc)|ef-"
    parse_labels() %>%          # c("a(bc)", "ef-", "-")
    lapply(parse_strand_label)  # list(c("a", "bc"), c("e", "f", "-"), c("-"))
  warp_threads <- parsed_labels[[1]]
  weft_threads <- parsed_labels[[2]]
  
  if (type == "basket") {
    n = n[1]
  }
  cell <- get_weave_pattern_matrix(type = type, n = n, warp_threads, weft_threads, 
                                   tie_up = tie_up, tr = tr, th = th) %>%
    make_polygons_from_matrix(spacing = spacing, margin = margin, aspect = aspect,
                              warp_threads, weft_threads, crs = crs)
  return(
    list(
      primitive = cell$weave_unit,
      transform = wk_affine_identity(),
      strands = unique(cell$strand),
      tile = cell$tile,
      type = type
    )
  )
}


############ TRY THESE!
# the below options make the pattern give or take an offset in the colours)
# from Figure 22 of Glassner 2002
# type <- "this"
this_tu <- matrix(c(0,0,0,0,1,1,1,1,
                    0,0,0,1,0,1,1,1,
                    0,0,1,0,1,0,1,1,
                    0,1,0,1,0,1,0,1,
                    1,0,1,0,1,0,1,0,
                    1,1,0,1,0,1,0,0,
                    1,1,1,0,1,0,0,0,
                    1,1,1,1,0,0,0,0), 8, 8, byrow = TRUE)
ids <- "aaaaaaaabbbbbbbb|aaaaaaaacccccccc"

