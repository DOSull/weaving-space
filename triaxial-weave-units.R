require(sf)

sf_use_s2(FALSE)

# Functions to generate sf datasets that are tileable to produce the
# appearance of a woven surface, with three directions of weaving thread.


# convenience storage of sqrt(3)
S3 <- sqrt(3)


get_triaxial_weave_matrices <- function(type = "cube", 
                                        strands_1 = c("a"),
                                        strands_2 = c("b"),
                                        strands_3 = c("c")) {

  if (type == "hex") {
    loom <- matrices_as_loom(
      get_weave_pattern_matrix(type = "this", tie_up = pracma::ones(6),
                               warp = strands_1, weft = strands_2),
      get_weave_pattern_matrix(type = "this", tie_up = pracma::ones(6),
                               warp = strands_2, weft = strands_3),
      get_weave_pattern_matrix(type = "this", tie_up = pracma::ones(6),
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
                                    strands = "a|b|c", type = "cube",
                                    crs = 3857) {
  margin_messages(margin, (1 - aspect) / 2 * sqrt(3) / 2)
  strand_ids <- get_strand_ids(strands)  
  strands_1 <- strand_ids[[1]]
  strands_2 <- strand_ids[[2]]
  strands_3 <- strand_ids[[3]]

  get_triaxial_weave_matrices(type = type,
                              strands_1, strands_2, strands_3) %>%
    make_sf_from_coded_weave_matrix(spacing = spacing,
                                    width = aspect, margin = margin,
                                    axis1_threads = strands_1,
                                    axis2_threads = strands_2,
                                    axis3_threads = strands_3, crs = crs)

  # list(
  #   primitive = cell$weave_unit,
  #   transform = wk::wk_affine_identity(),
  #   strands = unique(cell$weave_unit$strand),
  #   tile = cell$tile,
  #   type = type
  # )
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
  } else if (length(unique(scores)) != length(scores) && number_present != 1) {
    warning(paste0("Unable to determine a unique ordering on (",
                 paste0(orderings, collapse = ") ("), ")", collapse = ""))
    result <- NA
  } else {
    result <- as.vector(values[order(scores)[1:number_present]])
  }
  if (verbose) {
    return(list(result = result, ranks = ranks, scores = scores))
  } else {
    return(result)
  }
}
