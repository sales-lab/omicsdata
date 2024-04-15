#' Retrieve information about diseases.
#'
#' @return A data.frame containing disease abbreviations and descriptions.
#'
#' @examples
#' dplyr::filter(list_diseases(), abbreviation == "LUAD")
#'
#' @export
#'
list_diseases <- function() {
  url <- "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations"
  rvest::read_html(url) %>%
    rvest::html_element("#main table") %>%
    rvest::html_table() %>%
    select(
      abbreviation = "Study Abbreviation",
      description = "Study Name"
    ) %>%
    filter(!(.data$abbreviation %in% c("CNTL", "FPPP", "MISC")))
}
