#' Retrieve gene expression datasets for a disease.
#'
#' @param disease A disease code. See \code{\link{list_diseases}}.
#'
#' @export
#'
fetch_tcga_dataset <- function(disease) {
  listing <- safely(curatedTCGAData)(
    diseaseCode = disease,
    assays = "RNASeq2GeneNorm",
    dry.run = TRUE)
  if (!is.null(listing$error)) {
    abort(glue("no RNA-seq data for disease: {disease}"))
  }

  mdata <- curatedTCGAData(
    diseaseCode = disease,
    assays = "RNASeq2GeneNorm",
    dry.run = FALSE) %>%
    splitAssays("01") %>%
    mergeReplicates()

  expr <- assays(mdata[[1]])
  stopifnot(length(expr) == 1)
  expr <- expr[[1]]

  smap <- sampleMap(mdata)
  perm <- match(colnames(expr), smap$colname, nomatch = NA_integer_)
  stopifnot(all(!is.na(perm)))
  colnames(expr) <- smap$primary[perm]

  cdata <- colData(mdata) %>%
    as_tibble() %>%
    filter(.data$patientID %in% colnames(expr))
  stopifnot(nrow(cdata) == ncol(expr))

  list(
    metadata = as_tibble(cdata),
    expression = expr[, cdata$patientID]
  )
}
