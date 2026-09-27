# inferGenotypeAllele() joins the threshold table into the genotype counts with
# all.x = TRUE. A repeated join key therefore produced one output row per threshold
# row, each carrying the whole count, which inflated the per-locus depth. Depth is a
# shared denominator, so every row was wrong and z_score moved with it, not only the
# duplicated ones.

geno_data <- function(n1 = 22L, n2 = 20L) data.frame(
  v_call = c(rep("IGHVF1-G1*01", n1), rep("IGHVF1-G2*01", n2)),
  sequence_alignment = "ACGT", stringsAsFactors = FALSE)

share_table <- function(n_sharing) data.frame(
  allele     = c(paste0("IGHV1-2*0", seq_len(n_sharing)), "IGHV1-3*01"),
  asc_allele = c(rep("IGHVF1-G1*01", n_sharing), "IGHVF1-G2*01"),
  threshold  = 1e-4, stringsAsFactors = FALSE)

test_that("alleles sharing a cluster give one row, and depth is the read count", {
  # the cluster is the unit here, so none of this may vary with how many alleles
  # collapse into it
  for (n in 1:4) {
    r <- suppressMessages(inferGenotypeAllele(
      geno_data(), allele_threshold_table = share_table(n),
      asc_annotation = TRUE, find_unmutated = FALSE))

    expect_equal(nrow(r), 2L)
    expect_false(any(duplicated(r[, c("gene", "allele")])))
    expect_equal(unique(r$depth), 42)          # 22 + 20, not 22n + 20
    expect_equal(r[r$allele == "IGHVF1-G1*01"]$count, 22)
  }
})

test_that("a threshold table repeating an allele does not duplicate rows", {
  att <- data.frame(allele = c("IGHV1-2*02", "IGHV1-2*02", "IGHV1-3*01"),
                    asc_allele = c("IGHV1-2*02", "IGHV1-2*02", "IGHV1-3*01"),
                    threshold = 1e-4, stringsAsFactors = FALSE)
  d <- data.frame(v_call = c(rep("IGHV1-2*02", 20), rep("IGHV1-3*01", 20)),
                  sequence_alignment = "ACGT", stringsAsFactors = FALSE)

  r <- suppressMessages(inferGenotypeAllele(d, allele_threshold_table = att,
                                            find_unmutated = FALSE))
  expect_equal(nrow(r), 2L)
  expect_equal(unique(r$depth), 40)
})

test_that("an unlisted allele is added even when the table has no asc_allele column", {
  # the appended row used to name asc_allele literally, so a three-column table got a
  # four-column row and rbind failed from inside rbindlist, naming neither the caller
  # nor the column
  d <- data.frame(v_call = c(rep("IGHV1-2*02", 20), rep("IGHV1-3*01", 20),
                             rep("IGHV1-8*01", 20), rep("IGHV1-18*01", 20)),
                  sequence_alignment = "ACGT", stringsAsFactors = FALSE)
  att <- data.frame(allele = c("IGHV1-2*02", "IGHV1-3*01"), threshold = 1e-4,
                    stringsAsFactors = FALSE)

  r <- suppressWarnings(inferGenotypeAllele(d, allele_threshold_table = att,
                                            find_unmutated = FALSE))
  expect_equal(sum(unique(d$v_call) %in% r$allele), 4L)   # nothing dropped
  expect_equal(unique(r$depth), 80)
})
