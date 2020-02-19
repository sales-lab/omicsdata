#' Select the most abundant profiles from an expression matrix.
#'
#' @param expression The expression matrix.
#' @param n The number of profiles to select (100, by default).
#' @return A selection of rows from the input matrix.
#'
#' @export
#'
most_variable_profiles <- function(expression, n = 100L) {
  sel <- expression %>%
    apply(1, function(profile) stats::sd(profile) / mean(profile)) %>%
    sort(decreasing = TRUE) %>%
    utils::head(n = n) %>%
    names()

  expression[sel, , drop = FALSE]
}
