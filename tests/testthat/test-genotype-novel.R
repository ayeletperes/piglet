# Novel-allele support in genotype inference (tigger::findNovelAlleles style).

# Build a novel allele from a real germline by flipping one nucleotide, and a
# data set whose reads are unmutated against that novel sequence.
make_novel_fixture <- function(base = "IGHV1-2*02") {
  utils::data("HVGERM", package = "piglet", envir = environment())
  gl <- HVGERM[[base]]
  sv <- strsplit(gl, "")[[1]]
  pos <- which(sv %in% c("A", "C", "G", "T"))[120]
  sv[pos] <- ifelse(sv[pos] == "A", "T", "A")
  novel_seq <- paste0(sv, collapse = "")
  list(
    base = base,
    germline_db = HVGERM,
    novel = data.frame(germline_call = base,
                       polymorphism_call = paste0(base, "_G120A"),
                       novel_imgt = novel_seq,
                       stringsAsFactors = FALSE),
    novel_seq = novel_seq
  )
}

test_that("inferGenotypeAllele incorporates a novel allele with the base threshold", {
  skip_on_cran()
  fx <- make_novel_fixture()
  data("allele_threshold_table", package = "piglet")
  data <- data.frame(
    v_call = c(rep(fx$base, 30), rep("IGHV1-3*01", 40)),
    sequence_alignment = c(rep(fx$novel_seq, 30), rep(fx$germline_db[["IGHV1-3*01"]], 40)),
    stringsAsFactors = FALSE
  )
  g <- inferGenotypeAllele(data, allele_threshold_table = allele_threshold_table,
                           germline_db = fx$germline_db, novel = fx$novel,
                           find_unmutated = TRUE)

  novel_name <- fx$novel$polymorphism_call
  expect_true(any(g$allele == novel_name))
  base_thr  <- allele_threshold_table[allele == fx$base][["threshold"]]
  expect_equal(g[allele == novel_name][["threshold"]], base_thr)
})

test_that("inferGenotypeAllele_asc incorporates a novel allele with the base threshold", {
  skip_on_cran()
  fx <- make_novel_fixture()
  act <- data.frame(
    new_allele  = c(fx$base, "IGHV1-3*01"),
    func_group  = c("IGHV1-2", "IGHV1-3"),
    iuis_allele = c(fx$base, "IGHV1-3*01"),
    thresh      = c(1e-03, 1e-04),
    stringsAsFactors = FALSE
  )
  data <- data.frame(
    v_call = c(rep(fx$base, 30), rep("IGHV1-3*01", 20)),
    sequence_alignment = c(rep(fx$novel_seq, 30), rep(fx$germline_db[["IGHV1-3*01"]], 20)),
    stringsAsFactors = FALSE
  )
  asc <- suppressWarnings(
    inferGenotypeAllele_asc(data, alleleClusterTable = act,
                            germline_db = fx$germline_db, novel = fx$novel,
                            find_unmutated = TRUE))
  expect_true(any(grepl("_G120A", asc$iuis_alleles)))
  # the novel inherited the base allele threshold (1e-03), not the default 1e-04
  row <- asc[grepl("_G120A", asc$iuis_alleles)]
  expect_true(any(grepl("0.001", row$absolute_threshold)))
})

test_that("genotypeToTigger preserves the novel allele suffix", {
  geno <- data.frame(
    gene = c("IGHV1-2", "IGHV1-2"),
    allele = c("IGHV1-2*02", "IGHV1-2*02_G120A"),
    count = c(2, 30), depth = c(32, 32), threshold = c(1e-3, 1e-3),
    z_score = c(1, 100), stringsAsFactors = FALSE
  )
  tab <- genotypeToTigger(geno)
  expect_true(grepl("02_G120A", tab$alleles[tab$gene == "IGHV1-2"]))
})

test_that("allelePalette renders novel alleles as a transparent base color", {
  pal <- allelePalette(c("02", "02_G120A"))
  expect_match(pal[["02"]], "^#[0-9A-Fa-f]{6}$")          # opaque base
  expect_match(pal[["02_G120A"]], "^#[0-9A-Fa-f]{8}$")    # base + alpha
  # same RGB as the base, just with an alpha channel appended
  expect_equal(toupper(substr(pal[["02_G120A"]], 1, 7)),
               toupper(pal[["02"]]))
})
