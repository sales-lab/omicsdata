test_that("the most variable profiles are selected", {
  profs <- rbind(
    rnorm(100, 0, 1),
    rnorm(100, 0, 2),
    rep.int(0, 100)
  )
  rownames(profs) <- c("A", "B", "C")

  sel <- most_variable_profiles(profs, n = 2L)
  expect_is(sel, "matrix")
  expect_equal(rownames(sel), c("B", "A"))
})
