context("Filter for most variable profiles")

test_that("the most variable profiles are selected", {
  profs <- rbind(
    rnorm(100, 105, 1),
    rnorm(100, -110, 2),
    rnorm(100, 107, 2),
    rep.int(101, 100)
  )
  rownames(profs) <- c("A", "B", "C", "D")

  sel <- most_variable_profiles(profs, n = 2L)
  expect_is(sel, "matrix")
  expect_equal(rownames(sel), c("C", "A"))
})

test_that("low profiles are ignored", {
  profs <- rbind(
    rep.int(95, 100),
    rnorm(100, 50, 0.01),
    rnorm(100, 97, 30)
  )

  expect_error(most_variable_profiles(profs),
               class = "omicsdata_levels_too_low")
})


context("Sequencing Depth")

test_that("depth is labeled with sample names", {
  counts <- cbind(
    rep.int(1, 10),
    rep.int(2, 10),
    rep.int(3, 10)
  )
  colnames(counts) <- c("A", "B", "C")

  depth <- sequencing_depth(counts)
  expect_vector(depth, size = 3)
  expect_named(depth, colnames(counts))
  expect_equivalent(depth, c(10, 20, 30))
})
