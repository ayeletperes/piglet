# genotypeToTigger(): convert a PIgLET genotype to the TIgGER/VDJbase table layout.

# A mixed-locus, single-assignment table with a couple of low-count alleles so
# that a z-score threshold can move alleles in and out of the genotype.
table_data <- data.frame(
  v_call = c(rep("IGHV1-2*01", 12), rep("IGHV1-2*02", 2), rep("IGHV1-2*06", 1),
             rep("IGHV1-3*01", 40), rep("IGHV1-3*04", 9),
             rep("IGKV1-5*01", 30), rep("IGKV3-20*01", 10)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
table_att <- data.table::data.table(
  allele    = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-2*06", "IGHV1-3*01",
                "IGHV1-3*04", "IGKV1-5*01", "IGKV3-20*01"),
  threshold = rep(1e-04, 7)
)

test_that("genotypeToTigger reshapes inferGenotypeAllele output", {
  g <- inferGenotypeAllele(table_data, allele_threshold_table = table_att,
                           find_unmutated = FALSE)
  tab <- genotypeToTigger(g, level = "gene", z_threshold = 0)

  expect_s3_class(tab, "data.frame")
  expect_equal(names(tab),
               c("gene", "alleles", "counts", "total", "depth", "threshold", "z_score", "genotyped_alleles"))
  # genotyped_alleles column is lower case
  expect_true("genotyped_alleles" %in% names(tab))
  expect_false("GENOTYPED_ALLELES" %in% names(tab))

  ighv12 <- tab[tab$gene == "IGHV1-2", ]
  # alleles are numeric suffixes, ordered by descending count
  expect_equal(ighv12$alleles, "01,02,06")
  expect_equal(ighv12$counts, "12,2,1")
  expect_equal(ighv12$total, 15)
  # depth is the per-locus depth (IGH = 12+2+1+40+9 = 64), one per-allele threshold each
  expect_equal(ighv12$depth, 64)
  expect_equal(ighv12$threshold, "1e-04,1e-04,1e-04")
  expect_equal(tab$depth[tab$gene == "IGKV1-5"], 40)  # IGK = 30+10
})

test_that("genotyped_alleles is driven by the z-score threshold", {
  g <- inferGenotypeAllele(table_data, allele_threshold_table = table_att,
                           find_unmutated = FALSE)

  # at z >= 0 all positive-z alleles are present
  lo <- genotypeToTigger(g, z_threshold = 0)
  expect_equal(lo$genotyped_alleles[lo$gene == "IGHV1-2"], "01,02,06")

  # a high threshold drops the low-confidence alleles (z 24.9 and 12.4)
  hi <- genotypeToTigger(g, z_threshold = 100)
  expect_equal(hi$genotyped_alleles[hi$gene == "IGHV1-2"], "01")
  # IGHV1-3 alleles both clear 100
  expect_equal(hi$genotyped_alleles[hi$gene == "IGHV1-3"], "01,04")
})

test_that("z_score column is rounded to 3 decimals", {
  g <- inferGenotypeAllele(table_data, allele_threshold_table = table_att,
                           find_unmutated = FALSE)
  tab <- genotypeToTigger(g)
  z1 <- strsplit(tab$z_score[tab$gene == "IGKV1-5"], ",")[[1]][1]
  expect_equal(z1, "474.302")
})

test_that("genotypeToTigger accepts inferGenotypeAllele_asc output", {
  act <- data.frame(
    new_allele  = table_att$allele,
    func_group  = gsub("[*].*", "", table_att$allele),
    iuis_allele = table_att$allele,
    thresh      = rep(1e-04, 7),
    stringsAsFactors = FALSE
  )
  asc <- inferGenotypeAllele_asc(table_data, alleleClusterTable = act,
                                 find_unmutated = FALSE)

  tab_asc  <- genotypeToTigger(asc, level = "asc")
  tab_gene <- genotypeToTigger(asc, level = "gene")

  expect_equal(names(tab_asc),
               c("gene", "alleles", "counts", "total", "depth", "threshold", "z_score", "genotyped_alleles"))
  expect_equal(tab_gene$total[tab_gene$gene == "IGHV1-2"], 15)
})

test_that("level = 'asc' errors when the genotype carries no ASC information", {
  g <- inferGenotypeAllele(table_data, allele_threshold_table = table_att,
                           find_unmutated = FALSE)
  expect_error(genotypeToTigger(g, level = "asc"), "ASC information")
})

test_that("file = writes a tab-separated table that round-trips", {
  g <- inferGenotypeAllele(table_data, allele_threshold_table = table_att,
                           find_unmutated = FALSE)
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)

  tab <- genotypeToTigger(g, file = tmp)
  expect_true(file.exists(tmp))
  back <- data.table::fread(tmp, colClasses = "character")
  expect_equal(nrow(back), nrow(tab))
  expect_equal(names(back), names(tab))
})

test_that("genotypeToTigger rejects unrecognised input", {
  expect_error(genotypeToTigger(data.frame(foo = 1, bar = 2)),
               "Unrecognised genotype input")
})
