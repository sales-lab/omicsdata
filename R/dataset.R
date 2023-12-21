#' Retrieve gene expression datasets for a disease.
#'
#' @param disease A disease abbreviation. See \code{\link{list_diseases}}.
#' @return A \code{\link{RangedSummarizedExperiment}} instance.
#'
#' @export
#'
fetch_tcga_dataset <- function(disease) {
  if (!(disease %in% list_diseases()$abbreviation)) {
    abort(glue("invalid disease abbreviation: {disease}"))
  }

  dest_dir <- file.path(rappdirs::user_cache_dir("omicsdata"), disease)
  prepared_path <- file.path(dest_dir, "prepared.rda")
  if (file.exists(prepared_path))
    return(readRDS(prepared_path))

  dset <- gdc_download(disease, dest_dir)
  se <- select_assays(
    dset,
    c("unstranded", "tpm_unstrand"),
    c("raw_count", "normalized_count")
  )
  saveRDS(se, file = prepared_path)

  return(se)
}

gdc_download <- function(disease, dest_dir) {
  query <- GDCquery(
    project = glue("TCGA-", disease),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  if (is.null(query))
    abort(glue("no gene expression data for disease: {disease}"))

  GDCdownload(query, method = "api", directory = dest_dir)
  GDCprepare(remove_duplicates(query), directory = dest_dir)
}

remove_duplicates <- function(query) {
  results <- query$results[[1]]
  results <- results[which(!duplicated(results$cases)), ]

  query$results[[1]] <- results
  return(query)
}

select_assays <- function(dataset, selection, names_) {
  assays_ <- lapply(selection, \(n) assay(dataset, n))
  names(assays_) <- names_
  assays(dataset) <- assays_
  return(dataset)
}


#' Delete cached data about TCGA datasets.
#'
#' @return None
#'
#' @examples
#' purge_dataset_cache()
#'
#' @export
#'
purge_dataset_cache <- function() {
  base <- rappdirs::user_cache_dir("omicsdata")
  entries <- list.files(base)
  for (entry in entries) {
    unlink(file.path(base, entry), recursive = TRUE)
  }
}
