#' Select the most abundant profiles from an expression matrix.
#'
#' @param expression The expression matrix.
#' @param n The number of profiles to select.
#' @param cut_mean If the mean expression in a profile is lower than this value,
#'                 ignore the entire profile.
#' @return A selection of rows from the input matrix.
#'
#' @export
#'
most_variable_profiles <- function(expression, n = 100L, cut_mean = 100) {
  present <- rowMeans(expression) >= cut_mean
  if (sum(present) == 0) {
    abort("expression levels are too low for this filter to work",
          class = "omicsdata_levels_too_low")
  }

  sel <- expression[present, ] %>%
    apply(1, function(profile) stats::sd(profile) / mean(profile)) %>%
    sort(decreasing = TRUE) %>%
    utils::head(n = n) %>%
    names()

  expression[sel, , drop = FALSE]
}


#' Compute the per-sample sequencing depth.
#'
#' @param counts The count matrix.
#' @return An integer vector, holding the sequencing depth for each sample.
#'
#' @export
#'
sequencing_depth <- function(counts) {
  colSums(counts)
}
