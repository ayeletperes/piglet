# Experimental incremental / hysteresis ASC framework (R/asc_incremental.R).

# --- fixed transform is dataset-independent ------------------------------------

test_that(".asc_fixed_similarity depends only on pairwise distances", {
  # Equal-length sequences so hamming distance is a pure pairwise quantity.
  s3 <- c(A = "AAAACCCCGG", B = "AAAACCCCGT", C = "TTTTGGGGCC")
  s4 <- c(s3, D = "AAAACCCCTT")

  d3 <- as.matrix(igDistance(s3, method = "hamming"))
  d4 <- as.matrix(igDistance(s4, method = "hamming"))

  # existing pairwise distance unchanged when a 4th allele is added
  expect_equal(d4["A", "B"], d3["A", "B"])

  s3m <- .asc_fixed_similarity(d3)
  s4m <- .asc_fixed_similarity(d4)
  # the A-B edge weight is identical regardless of what else is in the set
  expect_equal(s4m["A", "B"], s3m["A", "B"])
  # linear transform stays in [0, 1] with a zero diagonal
  expect_true(all(s4m >= 0 & s4m <= 1))
  expect_equal(unname(diag(s4m)), rep(0, 4))
})

# --- retention -----------------------------------------------------------------

test_that("ascRetention computes max-overlap continuity per ASC", {
  old_assign <- c(a1 = "K1", a2 = "K1", a3 = "K1", a4 = "K1", b1 = "K2", b2 = "K2")
  # a4 moves out of K1's majority cluster; K2 stays intact (+ a new allele n1)
  new_assign <- c(a1 = "1", a2 = "1", a3 = "1", a4 = "2",
                  b1 = "3", b2 = "3", n1 = "1")

  ret <- ascRetention(old_assign, new_assign)
  k1 <- ret[ret$old_asc == "K1", ]
  k2 <- ret[ret$old_asc == "K2", ]

  expect_equal(k1$retention, 0.75)
  expect_equal(k1$best_new_cluster, "1")
  expect_equal(k1$n_moved, 1L)
  expect_equal(k1$n_new_spanned, 2L)     # clusters 1 and 2 both above min_share
  expect_equal(k2$retention, 1)
  expect_equal(k2$n_moved, 0L)
})

# --- change classification -----------------------------------------------------

test_that("ascPartitionChanges labels expansion / split / new_asc", {
  old_assign <- c(a1 = "k1", a2 = "k1", a3 = "k1", b1 = "k2", b2 = "k2")
  new_assign <- c(a1 = "10", a2 = "10", a3 = "11",   # k1 splits 10 | 11
                  b1 = "20", b2 = "20", n1 = "20",   # k2 expands with n1
                  n2 = "30", n3 = "30", n4 = "30")   # brand-new ASC 30
  new_alleles <- c("n1", "n2", "n3", "n4")

  ch <- ascPartitionChanges(old_assign, new_assign, new_alleles)

  k1 <- ch[ch$old_asc == "k1" & !is.na(ch$old_asc), ]
  expect_equal(k1$change_type, "candidate_split")
  expect_equal(k1$support, 0L)           # no new alleles drive the split

  k2 <- ch[ch$old_asc == "k2" & !is.na(ch$old_asc), ]
  expect_equal(k2$change_type, "simple_expansion")
  expect_equal(k2$support, 1L)           # n1

  na <- ch[ch$change_type == "new_asc", ]
  expect_equal(nrow(na), 1L)
  expect_equal(na$new_cluster, "30")
  expect_equal(na$support, 3L)           # n2, n3, n4
})

test_that("ascPartitionChanges detects a merge of two established ASCs", {
  old_assign <- c(a1 = "kA", a2 = "kA", b1 = "kB", b2 = "kB")
  new_assign <- c(a1 = "1", a2 = "1", b1 = "1", b2 = "1")   # kA and kB fused
  ch <- ascPartitionChanges(old_assign, new_assign, character(0))
  expect_true(all(ch$change_type == "candidate_merge"))
  expect_true(all(grepl("kA\\+kB", ch$members)))
})

# --- Leiden stability harness --------------------------------------------------

test_that("ascLeidenStability returns a valid co-assignment matrix", {
  data(HVGERM)
  d <- igDistance(HVGERM[1:12], method = "decipher")
  stab <- ascLeidenStability(d, resolution = 0.1, n_runs = 5)

  C <- stab$coassign
  expect_equal(dim(C), c(12, 12))
  expect_true(all(C >= 0 & C <= 1))
  expect_equal(unname(diag(C)), rep(1, 12))   # an allele always co-clusters with itself
  expect_true(isSymmetric(unname(C)))
  expect_length(stab$memberships, 5)
})

# --- within/between similarity -------------------------------------------------

test_that("ascWithinBetween separates tight groups from loose ones", {
  S <- matrix(0.2, 4, 4, dimnames = list(letters[1:4], letters[1:4]))
  S["a", "b"] <- S["b", "a"] <- 0.9
  S["c", "d"] <- S["d", "c"] <- 0.9
  diag(S) <- 0
  wb <- ascWithinBetween(S, list(c("a", "b"), c("c", "d")))
  expect_equal(wb$within, 0.9)
  expect_equal(wb$between, 0.2)
  expect_gt(wb$margin, 0)
})

# --- reconciliation / hysteresis ----------------------------------------------

test_that("ascReconcile preserves established ASCs unless confidence is High", {
  old_assign <- c(a1 = "k1", a2 = "k1", a3 = "k1")
  new_assign <- c(a1 = "10", a2 = "10", a3 = "11", n1 = "10")
  base_change <- data.frame(
    change_type = "candidate_split", old_asc = "k1", new_cluster = "10|11",
    n_old = 3L, n_new_alleles = 0L, retention = 2 / 3, support = 0L,
    members = "a1,a2,a3", n_pass = 0L, stringsAsFactors = FALSE)

  # Low confidence: the split is NOT applied; k1 stays intact.
  low <- transform(base_change, confidence = "Low")
  rec_low <- ascReconcile(old_assign, new_assign, low, "n1")
  expect_equal(rec_low$reconciled_asc[rec_low$allele == "a1"], "k1")
  expect_equal(rec_low$reconciled_asc[rec_low$allele == "a3"], "k1")
  # the new allele is still absorbed into its best original ASC
  expect_equal(rec_low$reconciled_asc[rec_low$allele == "n1"], "k1")

  # High confidence: the split IS applied; daughters get distinct labels.
  high <- transform(base_change, confidence = "High")
  rec_high <- ascReconcile(old_assign, new_assign, high, "n1")
  expect_equal(rec_high$reconciled_asc[rec_high$allele == "a1"], "k1.10")
  expect_equal(rec_high$reconciled_asc[rec_high$allele == "a3"], "k1.11")
})

test_that("ascClassifyConfidence is easy for new alleles, strict for established ASCs", {
  changes <- data.frame(
    change_type = c("simple_expansion", "new_asc", "candidate_split", "candidate_split"),
    old_asc = c("k1", NA, "k2", "k3"),
    new_cluster = c("1", "9", "2|3", "4|5"),
    n_old = c(2L, 0L, 4L, 4L), n_new_alleles = c(1L, 3L, 0L, 0L),
    retention = c(1, NA, 0.5, 0.5), support = c(1L, 3L, 5L, 5L),
    members = c("", "", "", ""), stringsAsFactors = FALSE)
  evidence <- data.frame(
    leiden_stability = c(NA, NA, 0.99, 0.30),
    resolution_stability = c(NA, NA, 1.0, 0.20),
    margin = c(NA, NA, 0.4, 0.01),
    silhouette_delta = c(NA, NA, 0.10, -0.05),
    subsample_stability = c(NA, NA, 0.95, 0.10))
  out <- ascClassifyConfidence(changes, evidence)
  expect_equal(out$confidence[1], "High")   # expansion is easy
  expect_equal(out$confidence[2], "High")   # new ASC is easy
  expect_equal(out$confidence[3], "High")   # strong reproducible split
  expect_equal(out$confidence[4], "Low")    # weak split
})

# --- end-to-end driver smoke test ---------------------------------------------

test_that("runIncrementalASCExperiment returns a populated result", {
  data(HVGERM)
  ref <- HVGERM[1:30]
  new <- HVGERM[31:38]
  res <- runIncrementalASCExperiment(ref, new, resolution = 0.1,
                                     n_leiden = 3, n_subsample = 2)

  expect_s3_class(res, "IncrementalASCResult")
  expect_equal(length(res$new_alleles), 8L)
  # reconciled table covers every allele in the expanded set
  expect_equal(nrow(res$reconciled), 38L)
  expect_true(all(c("summary_counts", "retention", "changes", "reconciled") %in%
                    names(res$diagnostics)))
  # every reconciled allele receives an ASC label
  expect_false(any(is.na(res$reconciled$reconciled_asc)))
  expect_output(print(res), "IncrementalASCResult")
})
