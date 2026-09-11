## Guards the regression that motivated fix/leiden-modal-consensus.
##
## igraph::cluster_leiden is randomised. Before this fix the package called it
## once, with no seed, so the same germline set could return different ASCs on
## consecutive runs: 10 seeds gave 3 distinct partitions on the human IGHV set.
## The original infer_asc_usofa.R took the modal partition over 100 seeds; that
## step was lost when the script moved into the package. These tests fail if it
## is lost again.

test_that("detect_communities_leiden is reproducible for a fixed seed", {
  set.seed(11)
  g <- igraph::sample_gnm(120, 900)
  igraph::E(g)$weight <- 1
  a <- igraph::membership(detect_communities_leiden(g, resolution = 0.05, seed = 3, n_runs = 5))
  b <- igraph::membership(detect_communities_leiden(g, resolution = 0.05, seed = 3, n_runs = 5))
  expect_identical(as.integer(a), as.integer(b))
})

test_that("the modal partition is stable across different starting seeds", {
  set.seed(12)
  g <- igraph::sample_sbm(150, matrix(c(.25, .02, .02, .02, .25, .02, .02, .02, .25), 3), rep(50, 3))
  igraph::E(g)$weight <- 1
  key <- function(s) paste(as.integer(igraph::membership(
    detect_communities_leiden(g, resolution = 0.05, seed = s, n_runs = 25))), collapse = ",")
  expect_equal(length(unique(vapply(c(1L, 5L, 17L), key, character(1)))), 1L)
})

test_that("inferAlleleClusters returns the same table twice on the leiden path", {
  skip_if_not(exists("HVGERM"))
  data(HVGERM, envir = environment())
  ref <- HVGERM[order(names(HVGERM))][1:120]
  a <- inferAlleleClusters(ref, locus = "IGHV", clustering_method = "leiden",
                           distance_method = "lv", n_runs = 10L, quiet = TRUE)
  b <- inferAlleleClusters(ref, locus = "IGHV", clustering_method = "leiden",
                           distance_method = "lv", n_runs = 10L, quiet = TRUE)
  ta <- a$alleleClusterTable; tb <- b$alleleClusterTable
  expect_identical(ta[order(ta[[1]]), ], tb[order(tb[[1]]), ])
})
