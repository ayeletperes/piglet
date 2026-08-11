# New 1.4.0 flags on inferGenotypeAllele():
#   - depth_adjusted_threshold: raise each threshold to a depth-aware floor.
#   - z_score_threshold / in_genotype: flag (do not drop) alleles passing a z-score.
#   - observed: mark alleles actually seen; genotypeToTigger() excludes unseen ones.
#   - find_unmutated on D/J calls: warn and skip (biased for trimmed segments).

# --- depth_adjusted_threshold --------------------------------------------------

# IGH-only, 50 reads => per-locus depth N = 50, so the depth floor 1/N = 0.02
# dominates the tiny per-allele threshold 1e-04.
igh50 <- data.frame(
  v_call = c(rep("IGHV1-2*01", 40), rep("IGHV3-7*01", 10)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
att50 <- data.table::data.table(
  allele    = c("IGHV1-2*01", "IGHV3-7*01"),
  threshold = rep(1e-04, 2)
)

test_that("depth_adjusted_threshold raises thresholds to the depth floor", {
  g  <- inferGenotypeAllele(igh50, allele_threshold_table = att50,
                            find_unmutated = FALSE)
  ga <- inferGenotypeAllele(igh50, allele_threshold_table = att50,
                            find_unmutated = FALSE, depth_adjusted_threshold = TRUE)

  # without the flag, thresholds stay at the supplied 1e-04
  expect_equal(unique(g$threshold), 1e-04)
  # with the flag, every threshold is lifted to max(1e-04, 1/50) = 0.02
  expect_equal(unique(ga$threshold), 1 / 50)
  expect_true(all(ga$threshold >= g$threshold))
  expect_true(any(ga$threshold > g$threshold))
})

# --- z_score_threshold / in_genotype ------------------------------------------

# Mixed-locus table with a range of counts so a z-score cutoff can separate
# high- from low-count alleles. IGH depth = 12+2+1+40+9 = 64, IGK depth = 40.
flag_data <- data.frame(
  v_call = c(rep("IGHV1-2*01", 12), rep("IGHV1-2*02", 2), rep("IGHV1-2*06", 1),
             rep("IGHV1-3*01", 40), rep("IGHV1-3*04", 9),
             rep("IGKV1-5*01", 30), rep("IGKV3-20*01", 10)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
flag_att <- data.table::data.table(
  allele    = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-2*06", "IGHV1-3*01",
                "IGHV1-3*04", "IGKV1-5*01", "IGKV3-20*01"),
  threshold = rep(1e-04, 7)
)

test_that("z_score_threshold flags alleles without dropping rows", {
  g0 <- inferGenotypeAllele(flag_data, allele_threshold_table = flag_att,
                            find_unmutated = FALSE)
  g1 <- inferGenotypeAllele(flag_data, allele_threshold_table = flag_att,
                            find_unmutated = FALSE, z_score_threshold = 50)

  expect_true("in_genotype" %in% names(g0))
  # every allele clears the default z_score_threshold = 0 (tiny threshold, count >= 1)
  expect_true(all(g0$in_genotype))
  # the flag is exactly z_score >= cutoff, and it actually discriminates here
  expect_equal(g1$in_genotype, g1$z_score >= 50)
  expect_true(any(g1$in_genotype) && any(!g1$in_genotype))
  # a low-count allele drops out while a high-count one stays flagged
  expect_false(g1$in_genotype[g1$allele == "IGHV1-2*06"])
  expect_true(g1$in_genotype[g1$allele == "IGHV1-3*01"])
  # flag only: no rows are removed
  expect_equal(nrow(g1), nrow(g0))
})

# --- observed + genotypeToTigger exclusion ------------------------------------

# IGH-only data; the threshold table also lists IGHV1-2*07, which never appears
# in the data.
obs_data <- data.frame(
  v_call = c(rep("IGHV1-2*01", 12), rep("IGHV1-2*02", 2), rep("IGHV1-2*06", 1),
             rep("IGHV1-3*01", 40), rep("IGHV1-3*04", 9)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
obs_att <- data.table::data.table(
  allele    = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-2*06", "IGHV1-3*01",
                "IGHV1-3*04", "IGHV1-2*07"),
  threshold = rep(1e-04, 6)
)

test_that("observed marks unseen alleles and genotypeToTigger excludes them", {
  g <- inferGenotypeAllele(obs_data, allele_threshold_table = obs_att,
                           find_unmutated = FALSE)

  expect_true("observed" %in% names(g))
  # the allele present only in the threshold table is not observed
  expect_false(g$observed[g$allele == "IGHV1-2*07"])
  # everything actually in the data is observed
  expect_true(all(g$observed[g$allele != "IGHV1-2*07"]))

  tab <- genotypeToTigger(g, level = "gene", z_threshold = 0)
  ighv12 <- tab[tab$gene == "IGHV1-2", ]
  # the unobserved *07 is not listed in the gene's allele/count summary
  expect_false(grepl("(^|,)07(,|$)", ighv12$alleles))
  expect_equal(ighv12$alleles, "01,02,06")
  expect_equal(ighv12$counts, "12,2,1")
})

# --- find_unmutated guard on D/J calls ----------------------------------------

d_data <- data.frame(
  d_call = c(rep("IGHD3-10*01", 30), rep("IGHD2-2*01", 10)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
d_att <- data.table::data.table(
  allele    = c("IGHD3-10*01", "IGHD2-2*01"),
  threshold = rep(1e-04, 2)
)

test_that("find_unmutated warns and is skipped for D/J segments", {
  expect_warning(
    gd <- inferGenotypeAllele(d_data, allele_threshold_table = d_att,
                              call = "d_call", find_unmutated = TRUE),
    "only meaningful for V"
  )
  # genotyping still completes (the unmutated step is skipped, not fatal)
  expect_true(nrow(gd) >= 1)
  expect_true(all(c("allele", "z_score", "in_genotype") %in% names(gd)))
})
