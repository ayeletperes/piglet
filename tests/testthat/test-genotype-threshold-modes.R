# Three ways to set the presence threshold. NULL used to load the bundled table and now
# means no table, so these pin the distinction rather than leaving it to the docs.

gd <- function() data.frame(
  v_call = c(rep("IGHV1-2*02", 60), rep("IGHV1-3*01", 40)),
  sequence_alignment = "ACGT", stringsAsFactors = FALSE)

test_that("NULL uses one constant threshold and only the alleles seen", {
  r <- suppressMessages(inferGenotypeAllele(gd(), find_unmutated = FALSE))
  expect_equal(nrow(r), 2L)                      # no pseudo-count rows for unseen alleles
  expect_equal(unique(r$threshold), 1e-4)
  expect_equal(unique(r$depth), 100)             # exactly the reads, no smoothing mass
})

test_that("default_allele_threshold sets the constant", {
  r <- suppressMessages(inferGenotypeAllele(gd(), find_unmutated = FALSE,
                                            default_allele_threshold = 5e-3))
  expect_equal(unique(r$threshold), 5e-3)
})

test_that("depth adaptation floors the constant at 1/N", {
  # the point of pairing NULL with adaptation: at depth N the null expects one read, so a
  # constant far below 1/N cannot make a single read look overwhelming
  r <- suppressMessages(inferGenotypeAllele(gd(), find_unmutated = FALSE,
                                            default_allele_threshold = 1e-9,
                                            depth_adjusted_threshold = TRUE))
  expect_equal(unique(r$threshold), 1 / 100)
  expect_true(all(is.finite(r$z_score)))
})

test_that("'default' asks for the bundled table explicitly", {
  r <- suppressMessages(inferGenotypeAllele(gd(), allele_threshold_table = "default",
                                            find_unmutated = FALSE))
  expect_gt(nrow(r), 100L)                       # the whole bundled reference
  expect_gt(length(unique(r$threshold)), 1L)     # per-allele, not constant
})

test_that("a data.frame still takes precedence", {
  att <- data.frame(allele = c("IGHV1-2*02", "IGHV1-3*01"), threshold = 5e-4,
                    stringsAsFactors = FALSE)
  r <- suppressMessages(inferGenotypeAllele(gd(), allele_threshold_table = att,
                                            find_unmutated = FALSE))
  expect_equal(unique(r$threshold), 5e-4)
  expect_equal(nrow(r), 2L)
})
