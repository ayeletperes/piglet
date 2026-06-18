# Per-locus repertoire depth in genotype inference.
#
# When heavy and light chains are present in the same table (e.g. single-cell
# AIRR data), each locus must be normalised against its own repertoire depth
# rather than a single global depth shared across all loci.

# A mixed-locus, single-assignment table: 100 heavy (IGH) reads and 40 light
# (IGK) reads. Per locus the depths are therefore 100 and 40.
mixed_data <- data.frame(
  v_call = c(rep("IGHV1-2*01", 80), rep("IGHV3-7*01", 20),
             rep("IGKV1-5*01", 30), rep("IGKV3-20*01", 10)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)

test_that("inferGenotypeAllele computes depth per locus", {
  att <- data.table::data.table(
    allele    = c("IGHV1-2*01", "IGHV3-7*01", "IGKV1-5*01", "IGKV3-20*01"),
    threshold = rep(1e-04, 4)
  )

  geno <- inferGenotypeAllele(mixed_data, allele_threshold_table = att,
                              find_unmutated = FALSE)

  igh_depth <- unique(geno$depth[grepl("^IGH", geno$allele)])
  igk_depth <- unique(geno$depth[grepl("^IGK", geno$allele)])

  # each locus gets its own depth, not one global value of 140
  expect_equal(igh_depth, 100)
  expect_equal(igk_depth, 40)
  expect_false(isTRUE(all.equal(igh_depth, igk_depth)))
})

test_that("inferGenotypeAllele depth is unchanged for single-locus input", {
  # IGH-only data with an IGH-only threshold table reproduces the previous
  # global behaviour: one depth equal to the total read count.
  igh_data <- mixed_data[grepl("^IGH", mixed_data$v_call), ]
  att_igh <- data.table::data.table(
    allele    = c("IGHV1-2*01", "IGHV3-7*01"),
    threshold = rep(1e-04, 2)
  )

  geno <- inferGenotypeAllele(igh_data, allele_threshold_table = att_igh,
                              find_unmutated = FALSE)

  expect_equal(unique(geno$depth), 100)
})

test_that("inferGenotypeAllele_asc normalises fractions per locus", {
  act <- data.frame(
    new_allele  = c("IGHV1-2*01", "IGHV3-7*01", "IGKV1-5*01", "IGKV3-20*01"),
    func_group  = c("IGHV1-2", "IGHV3-7", "IGKV1-5", "IGKV3-20"),
    iuis_allele = c("IGHV1-2*01", "IGHV3-7*01", "IGKV1-5*01", "IGKV3-20*01"),
    thresh      = rep(1e-04, 4),
    stringsAsFactors = FALSE
  )

  asc <- inferGenotypeAllele_asc(mixed_data, alleleClusterTable = act,
                                 find_unmutated = FALSE)

  frac <- setNames(as.numeric(asc$absolute_fraction), asc$gene)

  # light-chain alleles are normalised against the light-chain depth (40),
  # not the global depth (140): 30/40 = 0.75, 10/40 = 0.25.
  expect_equal(unname(frac["IGKV1-5"]), 0.75)
  expect_equal(unname(frac["IGKV3-20"]), 0.25)
  # heavy-chain alleles are normalised against the heavy-chain depth (100).
  expect_equal(unname(frac["IGHV1-2"]), 0.8)
  expect_equal(unname(frac["IGHV3-7"]), 0.2)
})
