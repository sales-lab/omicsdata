context("Filter for most variable profiles")

test_that("the most variable profiles are selected", {
  profs <- rbind(
    rnorm(100, 6, 1),
    rnorm(100, -10, 2),
    rnorm(100, 7, 2),
    rep.int(10, 100)
  )
  rownames(profs) <- c("A", "B", "C", "D")

  sel <- most_variable_profiles(profs, n = 2L)
  expect_is(sel, "matrix")
  expect_equal(rownames(sel), c("C", "A"))
})

test_that("low profiles are ignored", {
  profs <- rbind(
    rep.int(3, 100),
    rnorm(100, 1, 0.01),
    rnorm(100, 2, 30)
  )

  expect_error(most_variable_profiles(profs),
               class = "omicsdata_levels_too_low")
})
