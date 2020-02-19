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
    rep.int(90, 100),
    rnorm(100, 50, 0.01),
    rnorm(100, 90, 30)
  )

  expect_error(most_variable_profiles(profs),
               class = "omicsdata_levels_too_low")
})


context("Sequencing depth")

test_that("depth is labeled with sample names", {
  counts <- cbind(
    rep.int(1, 10),
    rep.int(2, 10),
    rep.int(3, 10)
  )
  colnames(counts) <- c("A", "B", "C")

  sd <- sequencing_depth(counts)
  expect_s3_class(sd, "tbl_df")
  expect_identical(sd$sample, colnames(counts))
  expect_equivalent(sd$depth, c(10, 20, 30))
})


context("Depth normalization")

test_that("raw counts are normalized by sequencin depth", {
  counts <- cbind(
    rep.int(1, 10),
    rep.int(2, 10),
    rep.int(3, 10)
  )
  colnames(counts) <- c("A", "B", "C")

  nc <- normalize_by_depth(counts)
  expect_is(nc, "matrix")
  expect_identical(colnames(nc), colnames(counts))
  expect_equivalent(colSums(nc), rep.int(1e6, 3))
})
