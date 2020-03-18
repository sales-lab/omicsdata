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
  if (file.exists(prepared_path)) {
    se <- readRDS(prepared_path)
  } else {
    se <- harmonize_datasets(
      list(gdc_download(disease, dest_dir, "results"),
           gdc_download(disease, dest_dir, "normalized_results")),
      c("raw_count", "normalized_count")
    )
    saveRDS(se, file = prepared_path)
  }

  se
}

gdc_download <- function(disease, dest_dir, type) {
  # Workaround for non-uniform sample.type labels.
  res <- try(gdc_download2(disease, dest_dir, "Primary solid Tumor", type))
  if (inherits(res, "try-error")) {
    gdc_download2(disease, dest_dir, "Primary Tumor", type)
  } else {
    return(res)
  }
}

gdc_download2 <- function(disease, dest_dir, sample_type, file_type) {
  query <- GDCquery(
    project = glue("TCGA-", disease),
    data.category = "Gene expression",
    data.type = "Gene expression quantification",
    sample.type = sample_type,
    file.type = file_type,
    legacy = TRUE
  )

  GDCdownload(query, directory = dest_dir)
  GDCprepare(query, directory = dest_dir)
}

harmonize_datasets <- function(datasets, assays) {
  samples <- sort(colnames(datasets[[1]]))
  meta <- colnames(colData(datasets[[1]]))
  template <- NULL
  selected_assays <- list()

  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    dataset <- dataset[, samples]
    colData(dataset) <- colData(dataset)[samples, meta]
    dataset <- simplify_rowdata(dataset)

    if (i == 1) {
      template <- dataset
    }

    if (!identical(colnames(dataset), samples)) {
      abort("mismatch in sample names")
    } else if (!identical(rowData(dataset), rowData(template))) {
      abort("mismatch in rowData")
    }

    dataset_cdata <- colData(dataset)
    template_cdata <- colData(template)
    if (!identical(rownames(dataset_cdata), rownames(template_cdata))) {
      abort("mismatch in colData row names")
    } else if (!identical(colnames(dataset_cdata), colnames(template_cdata))) {
      abort("mismatch in colData column names")
    }

    found_assays <- assayNames(dataset)
    common_assays <- intersect(found_assays, assays)
    if (length(common_assays) == 0) {
      abort("dataset with no usable assay")
    }

    for (assay_name in common_assays) {
      selected_assays[[assay_name]] <- assay(dataset, common_assays)
    }
  }

  assays(template) <- selected_assays
  template
}

simplify_rowdata <- function(se) {
  dat <- rowData(se)
  dat <- dat[, c("gene_id", "entrezgene", "ensembl_gene_id")]
  rowData(se) <- dat
  rownames(se) <- dat[, "gene_id"]
  se
}


#' Delete cached data about TCGA datasets.
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
