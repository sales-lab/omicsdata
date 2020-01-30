#' @importFrom curatedTCGAData curatedTCGAData
#' @importFrom dplyr as_tibble filter select
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom MultiAssayExperiment colData mergeReplicates sampleMap
#' @importFrom purrr safely
#' @importFrom rlang .data abort
#' @importFrom SummarizedExperiment assays
#' @importFrom TCGAutils splitAssays
NULL

#' Package development tools for R.
#'
#' More details on the course at \url{https://sales.bio.unipd.it/teaching.html}
#'
#' @docType package
#' @name omicsdata
"_PACKAGE"
