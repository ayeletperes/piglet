# inferGenotypeAllele(): the threshold table may be a data.frame or a data.table.

class_data <- data.frame(
  v_call = c(rep("IGHV1-2*01", 12), rep("IGHV1-2*02", 8), rep("IGHV1-3*01", 20)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
class_att <- data.frame(
  allele    = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-3*01"),
  threshold = rep(1e-04, 3),
  stringsAsFactors = FALSE
)

test_that("a threshold table is accepted as a data.frame or a data.table", {
  from_df <- inferGenotypeAllele(class_data, allele_threshold_table = class_att,
                                 find_unmutated = FALSE)
  from_dt <- inferGenotypeAllele(class_data,
                                 allele_threshold_table = data.table::as.data.table(class_att),
                                 find_unmutated = FALSE)

  expect_gt(nrow(from_df), 0)
  expect_equal(as.data.frame(from_df), as.data.frame(from_dt))
})

test_that("the caller's threshold table is left as it was", {
  before_df <- class_att
  before_dt <- data.table::as.data.table(class_att)

  invisible(inferGenotypeAllele(class_data, allele_threshold_table = class_att,
                                find_unmutated = FALSE))
  invisible(inferGenotypeAllele(class_data, allele_threshold_table = before_dt,
                                find_unmutated = FALSE))

  expect_equal(class_att, before_df)
  expect_equal(names(before_dt), c("allele", "threshold"))
})
