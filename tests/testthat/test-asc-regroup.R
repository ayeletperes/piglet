# regroupASCByLabel (R/asc_regroup.R)

make_tb <- function() data.frame(
  allele         = paste0("a", 1:5),
  allele_cluster = c("G1", "G1", "G2", "G3", "G3"),
  imgt_allele    = c("IGHD1-1*01", "IGHD1-2*01", "IGHD1-1*02", "IGHD2-2*01", "IGHD2-8*01"),
  duplicated     = c("", "", "", "IGHD2-8*01", ""),   # a4 is shared between 2-2 and 2-8
  stringsAsFactors = FALSE)

test_that("action='both' re-groups by label: splits mixed groups, reunites split labels, merges shared", {
  b <- regroupASCByLabel(make_tb(), action = "both")
  g <- function(al) b$regrouped[b$allele == al]
  expect_true("regrouped" %in% names(b))
  expect_false(g("a1") == g("a2"))   # IGHD1-1 and IGHD1-2 split out of the merged group G1
  expect_equal(g("a1"), g("a3"))     # IGHD1-1 split across G1/G2 is reunited
  expect_equal(g("a4"), g("a5"))     # 2-2 and 2-8 merged via the shared allele
  expect_true(grepl("IGHD2-2", g("a4")) && grepl("IGHD2-8", g("a4")))  # both names in the tag
  expect_length(unique(b$regrouped), 3L)  # {1-1}, {1-2}, {2-2/2-8}
})

test_that("action='split' only breaks mixed groups and does not reunite across groups", {
  s <- regroupASCByLabel(make_tb(), action = "split")
  g <- function(al) s$regrouped[s$allele == al]
  expect_false(g("a1") == g("a2"))   # G1 split by gene
  expect_false(g("a1") == g("a3"))   # IGHD1-1 in G1 and G2 stays split
})

test_that("action='merge' keeps mixed groups but merges shared-allele-linked genes", {
  m <- regroupASCByLabel(make_tb(), action = "merge")
  g <- function(al) m$regrouped[m$allele == al]
  expect_equal(g("a1"), g("a2"))     # G1 kept intact (not split)
  expect_equal(g("a4"), g("a5"))     # 2-2 and 2-8 still merged via shared allele
})

test_that("missing duplicated column degrades to label-only grouping", {
  tb <- make_tb(); tb$duplicated <- NULL
  b <- regroupASCByLabel(tb, action = "both")
  g <- function(al) b$regrouped[b$allele == al]
  expect_false(g("a4") == g("a5"))   # no shared-allele info -> 2-2 and 2-8 stay separate
  expect_equal(g("a1"), g("a3"))     # label reunion still works
})

test_that("a shared allele listed inside the label (comma / slash) forces a merge", {
  tb <- data.frame(
    asc_allele  = c("IGHD2-A*01", "IGHD2-B*01", "IGHD4-C*01"),
    iuis_allele = c("IGHD2-2*01,IGHD2-8*01_g1a", "IGHD2-8*02", "IGHD4-4*01/IGHD4-11*01"),
    stringsAsFactors = FALSE)
  b <- regroupASCByLabel(tb, group_col = "asc_allele", label_col = "iuis_allele",
                         duplicated_col = NULL, action = "both")
  g <- function(al) b$regrouped[b$asc_allele == al]
  expect_equal(g("IGHD2-A*01"), g("IGHD2-B*01"))   # 2-2/2-8 label links them to 2-8*02
  t4 <- b$regrouped[b$asc_allele == "IGHD4-C*01"]  # slash-shared allele
  expect_true(grepl("IGHD4-4", t4) && grepl("IGHD4-11", t4))  # both gene names kept in the tag
})

test_that("group column with an allele suffix collapses to the ASC gene", {
  tb <- data.frame(asc_allele = c("IGHD1-X*01", "IGHD1-X*02"),
                   iuis_allele = c("IGHD1-1*01", "IGHD1-1*02"), stringsAsFactors = FALSE)
  b <- regroupASCByLabel(tb, group_col = "asc_allele", label_col = "iuis_allele",
                         duplicated_col = NULL, action = "split")
  expect_equal(length(unique(b$regrouped)), 1L)    # same group, same gene -> one bin
})

test_that("errors on missing required columns", {
  expect_error(regroupASCByLabel(data.frame(x = 1), group_col = "allele_cluster"))
})
