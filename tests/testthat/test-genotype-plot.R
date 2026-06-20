# plotGenotypeAllele(): genotype bar plot with a per-allele z-score panel.

plot_data <- data.frame(
  v_call = c(rep("IGHV1-2*01", 12), rep("IGHV1-2*02", 2), rep("IGHV1-2*06", 1),
             rep("IGHV1-3*01", 40), rep("IGHV1-3*04", 9),
             rep("IGKV1-5*01", 30), rep("IGKV3-20*01", 10)),
  sequence_alignment = "ACGT",
  stringsAsFactors = FALSE
)
plot_att <- data.table::data.table(
  allele    = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-2*06", "IGHV1-3*01",
                "IGHV1-3*04", "IGKV1-5*01", "IGKV3-20*01"),
  threshold = rep(1e-04, 7)
)

test_that("plotGenotypeAllele returns a grob from a raw genotype", {
  skip_on_cran()
  g <- inferGenotypeAllele(plot_data, allele_threshold_table = plot_att,
                           find_unmutated = FALSE)
  grob <- plotGenotypeAllele(g, silent = TRUE)
  expect_s3_class(grob, "gtable")
})

test_that("plotGenotypeAllele accepts a genotypeToTigger() table", {
  skip_on_cran()
  g <- inferGenotypeAllele(plot_data, allele_threshold_table = plot_att,
                           find_unmutated = FALSE)
  tab <- genotypeToTigger(g)
  grob <- plotGenotypeAllele(tab, silent = TRUE)
  expect_s3_class(grob, "gtable")
})

test_that("allelePalette matches VDJbase colors and generates many distinct colors", {
  # curated VDJbase colors for the common alleles
  pal <- allelePalette(c("01", "02", "03"))
  expect_equal(unname(pal[c("01", "02", "03")]),
               c("#f5bc6e", "#9d69f4", "#598200"))

  # scales beyond the curated 01-19 set without running out or returning NA
  many <- allelePalette(sprintf("%02d", 1:40))
  expect_length(unique(many), 40)
  expect_false(any(is.na(many)))

  # novel alleles share the base allele color, drawn semi-transparent
  nov <- allelePalette(c("01", "01_G123A"))
  expect_match(nov[["01_G123A"]], "^#[0-9A-Fa-f]{8}$")          # base + alpha
  expect_equal(toupper(substr(nov[["01_G123A"]], 1, 7)),
               toupper(nov[["01"]]))
})

test_that("plotGenotypeAllele errors when no allele passes the z threshold", {
  skip_on_cran()
  g <- inferGenotypeAllele(plot_data, allele_threshold_table = plot_att,
                           find_unmutated = FALSE)
  expect_error(plotGenotypeAllele(g, z_threshold = 1e6, silent = TRUE),
               "No genotyped alleles")
})
