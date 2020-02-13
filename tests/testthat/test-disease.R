test_that("list_diseases returns a two-column tibble", {
  x <- list_diseases()
  expect_s3_class(x, "tbl_df")
  expect_identical(colnames(x), c("abbreviation", "description"))
})

test_that("the list of diseases includes LUAD", {
  x <- list_diseases()
  expect_true("LUAD" %in% x$abbreviation)
})
