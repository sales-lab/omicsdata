#' Retrieve information about diseases.
#'
#' @export
#'
list_diseases <- function() {
  TCGAutils::diseaseCodes %>%
    select(
      abbreviation = .data$Study.Abbreviation,
      description = .data$Study.Name
    )
}
