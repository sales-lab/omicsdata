#' Retrieve information about diseases.
#'
#' @export
#'
list_diseases <- function() {
  TCGAutils::diseaseCodes %>%
    as_tibble() %>%
    select(
      abbreviation = .data$Study.Abbreviation,
      description = .data$Study.Name
    ) %>%
    filter(!(.data$abbreviation %in% c("CNTL", "FPPP", "MISC")))
}
