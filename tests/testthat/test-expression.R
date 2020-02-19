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


context("Expression medians")

test_that("gene medians are computed correctly", {
  counts <- rbind(
    1:3,
    2:4,
    4:6
  )
  rownames(counts) <- c("A", "B", "C")

  ms <- gene_medians(counts)
  expect_identical(names(ms), rownames(counts))
  expect_equivalent(ms, c(2, 3, 5))
})


context("Gene deviations")

test_that("gene deviations are computed correctly", {
  counts <- rbind(
    1:3,
    2:4,
    4:6
  )
  colnames(counts) <- c("S", "T", "U")
  rownames(counts) <- c("A", "B", "C")

  refs <- setNames(c(2, 3, 0), rownames(counts))

  ds <- gene_deviations(counts, refs)
  expect_is(ds, "data.frame")
  expect_named(ds, c("sample", "deviation"))
  expect_identical(as.character(unique(ds$sample)), unique(colnames(counts)))
})

test_that("mismatched reference levels are rejected", {
  counts <- rbind(
    1:3,
    2:4,
    4:6
  )
  rownames(counts) <- c("A", "B", "C")

  refs <- setNames(1:3, c("S", "T", "U"))

  testthat::expect_error(gene_deviations(counts, refs))
})
