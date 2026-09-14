# The old names are called from roughly forty analysis scripts outside this package,
# most not under version control. These tests are what stop a future edit from
# silently changing their numbers.

test_that("renamed functions keep the old behaviour", {
  germs  <- c("ACGTACGT", "ACGTACGT", "ACGTACGT")
  inputs <- c("ACGTTCGT", "ACGTTCGA", "ACGT.CGT")

  # a gap on the input side is not a mismatch
  expect_equal(allele_diff_paired(germs, inputs, return_count = TRUE), c(1L, 2L, 0L))
  expect_equal(allele_diff_paired(germs, inputs)[[2]], c(5L, 8L))

  expect_equal(insert_gaps(c("caggtc..aact", "ac..gt"), c("caggtcaact", "acgt")),
               c("caggtc..aact", "ac..gt"))
})

test_that("unequal lengths are padded, not read past the end", {
  # The pre-rename allele_diff_indices_parallel looped to the germline length while
  # indexing the input, so a short input was read out of bounds and returned a count
  # drawn partly from adjacent memory (160 vs 8 characters gave 144).
  germ  <- paste(rep("ACGT", 40), collapse = "")
  short <- "ACGTACGT"
  expect_equal(allele_diff_paired(germ, short, return_count = TRUE), 0L)
  expect_equal(allele_diff_paired(short, germ, return_count = TRUE), 0L)
})

test_that("deprecated aliases still work and still warn", {
  germs  <- c("ACGTACGT", "ACGTACGT")
  inputs <- c("ACGTTCGT", "ACGTTCGA")

  expect_warning(r <- allele_diff_indices_parallel2(germs, inputs, return_count = TRUE),
                 "deprecated")
  expect_equal(r, allele_diff_paired(germs, inputs, return_count = TRUE))

  # the dead `parallel` argument is still accepted so existing calls do not error
  expect_warning(r2 <- allele_diff_indices_parallel2(germs, inputs, parallel = TRUE,
                                                     return_count = TRUE))
  expect_equal(r2, r)

  expect_warning(g <- insert_gaps2_vec(c("ac..gt"), c("acgt")), "deprecated")
  expect_equal(g, insert_gaps(c("ac..gt"), c("acgt")))

  # this one's counts deliberately changed; the warning has to say so
  expect_warning(allele_diff_indices_parallel(germs, inputs, return_count = TRUE),
                 "counts")
})
