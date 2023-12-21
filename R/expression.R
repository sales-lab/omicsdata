#' Select the most abundant profiles from an expression matrix.
#'
#' @param expression The expression matrix.
#' @param n The number of profiles to select.
#' @param cut_mean If the mean expression in a profile is lower than this value,
#'                 ignore the entire profile.
#' @return A selection of rows from the input matrix.
#'
#' @examples
#' library(airway)
#' data(airway)
#' most_variable_profiles(assays(airway)$counts)
#'
#' @export
#'
most_variable_profiles <- function(expression, n = 100L, cut_mean = 100) {
  present <- rowMeans(expression) >= cut_mean
  if (sum(present) == 0) {
    abort("expression levels are too low for this filter to work",
          class = "omicsdata_levels_too_low")
  }

  sel <- expression[present, , drop = FALSE] %>%
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
#' @examples
#' library(airway)
#' data(airway)
#' sequencing_depth(assays(airway)$counts)
#'
#' @export
#'
sequencing_depth <- function(counts) {
  tibble(sample = colnames(counts),
         depth = colSums(counts))
}


#' Normalize expression profiles by sequencing depth.
#'
#' @param counts The raw count matrix.
#' @return Trasformed counts.
#'
#' @examples
#' library(airway)
#' data(airway)
#' normalize_by_depth(assays(airway)$counts)
#'
#' @export
#'
normalize_by_depth <- function(counts) {
  t(t(counts) / colSums(counts)) * 1e6
}


#' Compute median expression levels for all genes.
#'
#' @param expr The expression matrix.
#' @return A named vector with gene medians.
#'
#' @examples
#' library(airway)
#' data(airway)
#' gene_medians(assays(airway)$counts)
#'
#' @export
#'
gene_medians <- function(expr) {
  setNames(matrixStats::rowMedians(expr), rownames(expr))
}


#' Compute deviations of gene expression for reference levels.
#'
#' @param expr The expression matrix.
#' @param reference_levels A vector with a reference level for each gene.
#' @return A table of deviations, with columns \code{sample} and
#'         \code{deviation}.
#'
#' @examples
#' library(airway)
#' data(airway)
#' counts <- assays(airway)$counts
#' m <- gene_medians(counts)
#' gene_deviations(counts, m)
#'
#' @export
#'
gene_deviations <- function(expr, reference_levels) {
  assert_that(are_equal(rownames(expr), names(reference_levels)))

  for (i in seq_len(nrow(expr))) {
    expr[i, ] <- expr[i, ] - reference_levels[[i]]
  }

  as.data.frame(expr) %>%
    utils::stack() %>%
    select(sample = "ind", deviation = "values")
}
