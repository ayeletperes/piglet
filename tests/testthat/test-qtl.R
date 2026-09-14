# The two closed forms are the only places these scans can be silently wrong: everything
# else is bookkeeping that shows up as a missing column. So they are checked against the
# functions whose answer is not in doubt, lm() and stats::manova(), on planted effects.

make_data <- function(n = 120L, seed = 1L) {
  set.seed(seed)
  dosage <- rbind(
    v_effect = rep(c(0, 1, 2), length.out = n),
    v_null   = sample(c(0, 1, 2), n, replace = TRUE),
    v_miss   = c(NA, rep(c(0, 1, 2), length.out = n - 1L)))
  colnames(dosage) <- sprintf("s%03d", seq_len(n))
  pheno <- cbind(
    p_hit  = 0.8 * dosage["v_effect", ] + stats::rnorm(n),
    p_flat = stats::rnorm(n))
  rownames(pheno) <- colnames(dosage)
  list(dosage = dosage, pheno = pheno)
}

test_that("qtlScanUnivariate reproduces lm()", {
  d <- make_data()
  got <- qtlScanUnivariate(d$pheno, d$dosage, min_n = 60L)

  for (v in c("v_effect", "v_null")) {
    for (p in c("p_hit", "p_flat")) {
      fit <- summary(stats::lm(d$pheno[, p] ~ d$dosage[v, ]))$coefficients[2L, ]
      row <- got[got$variant == v & got$phenotype == p]
      expect_equal(row$beta, unname(fit[["Estimate"]]), tolerance = 1e-10)
      expect_equal(row$se, unname(fit[["Std. Error"]]), tolerance = 1e-10)
      expect_equal(row$t_stat, unname(fit[["t value"]]), tolerance = 1e-10)
      expect_equal(row$p_value, unname(fit[["Pr(>|t|)"]]), tolerance = 1e-10)
      expect_equal(row$n, ncol(d$dosage))
    }
  }
  # The planted effect is found. There is deliberately no assertion that the null
  # variant is non-significant: that is a draw from a uniform, and asserting it would
  # make the test fail on an unlucky seed rather than on a broken scan.
  expect_lt(got[got$variant == "v_effect" & got$phenotype == "p_hit"]$p_value, 1e-10)
})

test_that("a variant with a missing call is fitted on its own complete cases", {
  d <- make_data()
  got <- qtlScanUnivariate(d$pheno, d$dosage, min_n = 60L)
  row <- got[got$variant == "v_miss" & got$phenotype == "p_hit"]
  expect_equal(row$n, ncol(d$dosage) - 1L)

  keep <- !is.na(d$dosage["v_miss", ])
  fit <- summary(stats::lm(d$pheno[keep, "p_hit"] ~ d$dosage["v_miss", keep]))$coefficients[2L, ]
  expect_equal(row$beta, unname(fit[["Estimate"]]), tolerance = 1e-10)
  expect_equal(row$p_value, unname(fit[["Pr(>|t|)"]]), tolerance = 1e-10)
})

test_that("qtlScanMultivariate reproduces stats::manova() Pillai", {
  d <- make_data()
  pheno <- cbind(d$pheno, p_third = stats::rnorm(nrow(d$pheno)))
  rownames(pheno) <- rownames(d$pheno)
  got <- qtlScanMultivariate(pheno, d$dosage, min_n = 60L)

  for (v in c("v_effect", "v_null")) {
    g <- d$dosage[v, ]
    ref <- summary(stats::manova(pheno ~ g), test = "Pillai")$stats[1L, ]
    row <- got[got$variant == v]
    expect_equal(row$pillai, unname(ref[["Pillai"]]), tolerance = 1e-10)
    expect_equal(row$f_stat, unname(ref[["approx F"]]), tolerance = 1e-10)
    expect_equal(row$p_value, unname(ref[["Pr(>F)"]]), tolerance = 1e-10)
    expect_equal(row$n_response, ncol(pheno))
  }
})

test_that("qtlLDGroups collapses duplicated and mirrored variants but not distinct ones", {
  dosage <- rbind(a = c(0, 1, 2, 0, 1, 2), copy = c(0, 1, 2, 0, 1, 2),
                  mirror = c(2, 1, 0, 2, 1, 0), other = c(0, 0, 1, 2, 2, 1))
  colnames(dosage) <- sprintf("s%d", 1:6)
  g <- qtlLDGroups(dosage, rep("chr1", 4))
  expect_equal(g[["a"]], g[["copy"]])
  expect_equal(g[["a"]], g[["mirror"]])       # identical up to sign
  expect_false(g[["a"]] == g[["other"]])

  # A contig prefix stops the same pattern at two loci from merging: positions are
  # locus-relative and must never be pooled.
  expect_false(qtlLDGroups(dosage, c("chr1", "chr2", "chr1", "chr1"))[["a"]] ==
                 qtlLDGroups(dosage, c("chr1", "chr2", "chr1", "chr1"))[["copy"]])
})

test_that("qtlSmallestGenotypeClass ignores classes with no subject", {
  dosage <- rbind(all_three = c(0, 0, 0, 1, 1, 2),
                  two_only  = c(0, 0, 0, 0, 1, 1))
  colnames(dosage) <- sprintf("s%d", 1:6)
  got <- qtlSmallestGenotypeClass(dosage)
  expect_equal(unname(got[["all_three"]]), 1)   # one homozygote
  expect_equal(unname(got[["two_only"]]), 2)    # the empty class is not counted as 0
})

test_that("pairingScan refuses to run without a conditional", {
  rep_dt <- data.frame(
    subject = rep(sprintf("s%02d", 1:4), each = 4),
    d = rep(c("D1", "D1", "D2", "D2"), 4),
    j = rep(c("J1", "J2"), 8))
  pairs <- suppressMessages(pairingTable(rep_dt, anchor = "j", partner = "d"))
  dosage <- matrix(c(0, 1, 1, 2), nrow = 1, dimnames = list("v1", sprintf("s%02d", 1:4)))
  expect_error(pairingScan(pairs, dosage), "conditional")
  expect_error(pairingScan(pairs, dosage, conditional = ""), "conditional")
})

test_that("pairing enrichment is symmetric under swapping anchor and partner", {
  set.seed(2)
  rep_dt <- data.frame(
    subject = rep(sprintf("s%02d", 1:5), each = 40),
    d = sample(c("D1", "D2", "D3"), 200, replace = TRUE),
    j = sample(c("J1", "J2"), 200, replace = TRUE))
  a <- suppressMessages(pairingTable(rep_dt, anchor = "j", partner = "d"))
  b <- suppressMessages(pairingTable(rep_dt, anchor = "d", partner = "j"))
  m <- merge(a, b, by.x = c("subject", "anchor_gene", "partner_gene"),
             by.y = c("subject", "partner_gene", "anchor_gene"))
  expect_equal(nrow(m), nrow(a))
  expect_equal(m$enrichment.x, m$enrichment.y, tolerance = 1e-12)
  expect_equal(m$expected.x, m$expected.y, tolerance = 1e-12)
  # The conditionals are not symmetric, and are the ones that flip.
  expect_equal(m$p_partner_given_anchor.x, m$p_anchor_given_partner.y, tolerance = 1e-12)
})

test_that("the scalar and vectorised Pillai agree", {
  # The closed form is written twice in qtl_scan.R: vectorised over variants inside
  # qtlScanMultivariate, and for a single variant in .qtl_pillai, which the per-lead
  # follow-ups need. Keeping both is deliberate - delegating the scalar one to the
  # vectorised one measured 2x slower on the leave-one-out load - but a formula written
  # twice can be corrected once, and then pairingLeadCharacter's leave-one-out and
  # heterogeneity p-values would silently stop agreeing with the omnibus they follow up.
  set.seed(3)
  n <- 90L
  y <- matrix(stats::rnorm(n * 6L), n, 6L)
  rownames(y) <- sprintf("s%03d", seq_len(n))
  g <- rep(c(0, 1, 2), length.out = n)
  d <- matrix(g, nrow = 1L, dimnames = list("v1", rownames(y)))

  scalar <- piglet:::.qtl_pillai(y, g)
  vectorised <- qtlScanMultivariate(y, d, min_n = 60L)

  expect_equal(unname(scalar[["pillai"]]), vectorised$pillai, tolerance = 1e-12)
  expect_equal(unname(scalar[["f"]]), vectorised$f_stat, tolerance = 1e-12)
  expect_equal(unname(scalar[["p"]]), vectorised$p_value, tolerance = 1e-12)
})

test_that("a single-variant dosage matrix works", {
  # vapply simplifies to a vector when there is one row, and apply() / max.col() below
  # then have no dim to work on. One variant is the normal case when a caller inspects a
  # single lead, and every pairing function goes through both of these.
  d <- matrix(rep(c(0, 1, 2), length.out = 12L), nrow = 1L,
              dimnames = list("v1", sprintf("s%02d", 1:12)))
  expect_equal(unname(qtlSmallestGenotypeClass(d)), 4)

  y <- matrix(stats::rnorm(12L * 3L), 12L, 3L,
              dimnames = list(colnames(d), c("p1", "p2", "p3")))
  cm <- piglet:::.qtl_class_means(y, d)
  expect_equal(nrow(cm), 3L)             # one variant x three phenotypes
  expect_false(anyNA(cm$n_low))

  # and the two-variant result is unchanged
  d2 <- rbind(a = rep(c(0, 1, 2), length.out = 12L),
              b = rep(c(0, 0, 1, 2), length.out = 12L))
  colnames(d2) <- colnames(d)
  expect_equal(unname(qtlSmallestGenotypeClass(d2)), c(4, 3))
})
