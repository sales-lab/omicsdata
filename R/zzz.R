#' @importFrom assertthat assert_that are_equal
#' @importFrom dplyr as_tibble select filter
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom rlang .data abort
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay assays assayNames assays<- colData
#'             colData<- rowData rowData<-
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare
#' @importFrom tibble tibble
NULL

#' Support package for the "Omics Data" course.
#'
#' @docType package
#' @name omicsdata
"_PACKAGE"
