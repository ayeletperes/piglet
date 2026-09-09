# =============================================================================
# Experimental incremental / hysteresis Allele Similarity Cluster (ASC) framework
# -----------------------------------------------------------------------------
# When new alleles are added to a reference - a large batch, or a new species
# with unknown orthology - re-running Leiden community detection from scratch can
# reshuffle established ASCs even when the new evidence is weak. This module adds
# an experimental framework that:
#   * runs the "add a batch, re-cluster, compare" experiment,
#   * scores how much of each original ASC is retained,
#   * classifies changes (expansion / new ASC / candidate split / merge / reassign),
#   * subjects changes to established ASCs to strong-evidence tests,
#   * classifies each change High / Intermediate / Low confidence, and
#   * emits a reconciled ASC table that applies changes only where confidence is
#     High (soft stability / hysteresis).
#
# Design note - dataset independence: the production distance_to_graph() rescales
# distances by max(distance) and similarities by their max, so adding sequences
# changes the edge weights between two *existing* alleles. Here the transform is
# isolated and fixed (S = 1 - d, or a fixed-floor -log), so the weight between two
# alleles depends only on their pairwise distance and is invariant to what else is
# in the set. igDistance(method = "decipher") returns a pairwise fraction in
# [0, 1], which is the assumption this module relies on.
# =============================================================================

# ---- 1. Dataset-independent distance -> similarity -> graph ------------------

#' Fixed, dataset-independent distance-to-similarity transform
#'
#' Converts a pairwise distance matrix to a similarity matrix without any
#' dataset-wide normalisation, so the similarity between two alleles depends only
#' on their pairwise distance and never changes when other alleles are added or
#' removed. This is the key property that lets ASCs stay stable as the reference
#' grows.
#'
#' @param distance_matrix A square distance matrix (or a \code{dist} object) with
#'   values in \[0, 1\] (as returned by \code{igDistance(method = "decipher")}).
#' @param method \code{"linear"} (default) gives \code{S = 1 - d}; \code{"log"}
#'   gives \code{S = min(-log(max(d, floor)), -log(floor)) / -log(floor)}.
#' @param floor Constant lower bound used by the \code{"log"} method (default
#'   \code{1e-3}). It is a fixed constant, never derived from the data.
#'
#' @return A square similarity matrix with a zero diagonal and the same dimnames
#'   as the input.
#'
#' @keywords internal
.asc_fixed_similarity <- function(distance_matrix,
                                  method = c("linear", "log"),
                                  floor = 1e-3) {
  method <- match.arg(method)
  if (inherits(distance_matrix, "dist")) distance_matrix <- as.matrix(distance_matrix)
  stopifnot(is.matrix(distance_matrix),
            nrow(distance_matrix) == ncol(distance_matrix))

  D <- distance_matrix
  if (method == "linear") {
    S <- 1 - D
  } else {
    denom <- -log(floor)
    S <- pmin(-log(pmax(D, floor)), denom) / denom
  }
  S[!is.finite(S)] <- 0
  S[S < 0] <- 0
  diag(S) <- 0
  dimnames(S) <- dimnames(distance_matrix)
  S
}

#' Build a weighted graph from a fixed similarity matrix
#'
#' @param similarity A square similarity matrix (e.g. from
#'   \code{.asc_fixed_similarity}).
#' @param min_similarity Fixed pruning threshold; edges with similarity below
#'   this constant are dropped (default \code{0}, i.e. keep all). The threshold is
#'   a constant, not derived from the data.
#'
#' @return A weighted undirected \code{igraph} graph.
#'
#' @keywords internal
.asc_build_graph <- function(similarity, min_similarity = 0) {
  S <- similarity
  if (min_similarity > 0) S[S < min_similarity] <- 0
  diag(S) <- 0
  igraph::graph_from_adjacency_matrix(S, mode = "undirected",
                                      weighted = TRUE, diag = FALSE)
}

#' Allele labels of a distance matrix / dist object
#' @keywords internal
.asc_labels <- function(distance_matrix) {
  if (inherits(distance_matrix, "dist")) return(attr(distance_matrix, "Labels"))
  rownames(distance_matrix)
}

# ---- 2. Clustering + stability harness --------------------------------------

#' Single Leiden clustering with the fixed transform
#'
#' Runs one Leiden (CPM) community detection on the fixed, dataset-independent
#' similarity graph. Because \code{igraph}'s Leiden is randomised, \code{seed}
#' makes a run reproducible and lets callers sample the run-to-run distribution.
#'
#' @param distance_matrix A square distance matrix (or \code{dist}).
#' @param resolution Leiden CPM resolution.
#' @param seed Optional integer seed (\code{set.seed} before clustering).
#' @param sim_method,floor Passed to \code{.asc_fixed_similarity}.
#' @param min_similarity Passed to \code{.asc_build_graph}.
#'
#' @return A named integer membership vector (names = allele labels).
#'
#' @examples
#' data(HVGERM)
#' d <- igDistance(HVGERM[1:20], method = "decipher")
#' m <- ascClusterOnce(d, resolution = 0.1, seed = 1)
#'
#' @export
ascClusterOnce <- function(distance_matrix, resolution, seed = NULL,
                           sim_method = "linear", floor = 1e-3,
                           min_similarity = 0) {
  S <- .asc_fixed_similarity(distance_matrix, method = sim_method, floor = floor)
  g <- .asc_build_graph(S, min_similarity = min_similarity)
  if (!is.null(seed)) set.seed(seed)
  comm <- detect_communities_leiden(g, resolution = resolution)
  memb <- as.integer(igraph::membership(comm))
  names(memb) <- igraph::V(g)$name
  labs <- .asc_labels(distance_matrix)
  if (!is.null(labs)) memb <- memb[labs]
  memb
}

#' Sweep resolutions and score partitions by silhouette
#'
#' Module-local resolution sweep using the fixed transform (the production
#' \code{optimize_resolution} hardcodes the dataset-dependent transform, so it is
#' not reused here). Returns per-resolution cluster counts and mean silhouette.
#'
#' @keywords internal
.asc_resolution_sweep <- function(distance_matrix, resolutions, seed = 1,
                                  sim_method = "linear", floor = 1e-3,
                                  min_similarity = 0) {
  D <- if (inherits(distance_matrix, "dist")) distance_matrix else stats::as.dist(distance_matrix)
  out <- data.frame(resolution = resolutions, n_clusters = NA_integer_,
                    silhouette = NA_real_)
  parts <- vector("list", length(resolutions))
  for (i in seq_along(resolutions)) {
    m <- ascClusterOnce(distance_matrix, resolutions[i], seed = seed,
                        sim_method = sim_method, floor = floor,
                        min_similarity = min_similarity)
    parts[[i]] <- m
    k <- length(unique(m))
    out$n_clusters[i] <- k
    if (k > 1 && k < length(m)) {
      si <- cluster::silhouette(as.integer(factor(m)), D)
      if (!any(is.na(si))) out$silhouette[i] <- mean(si[, 3])
    }
  }
  best <- which.max(replace(out$silhouette, is.na(out$silhouette), -Inf))
  list(results = out, partitions = parts,
       best_resolution = resolutions[best], best_partition = parts[[best]])
}

#' Leiden stability across seeded repetitions
#'
#' Runs Leiden \code{n_runs} times with different seeds and summarises how often
#' each pair of alleles is placed in the same community (the consensus /
#' co-assignment matrix). This is the harness behind the "run Leiden 100x" part
#' of the strong-evidence tests.
#'
#' @param distance_matrix A square distance matrix (or \code{dist}).
#' @param resolution Leiden CPM resolution.
#' @param n_runs Number of seeded repetitions (default 100).
#' @param seeds Optional integer seeds; defaults to \code{seq_len(n_runs)}.
#' @param ... Passed to \code{ascClusterOnce}.
#'
#' @return A list with \code{memberships} (list of membership vectors),
#'   \code{coassign} (an allele x allele matrix of co-clustering fractions in
#'   \[0, 1\], unit diagonal), and \code{n_runs}.
#'
#' @export
ascLeidenStability <- function(distance_matrix, resolution, n_runs = 100,
                               seeds = NULL, ...) {
  if (is.null(seeds)) seeds <- seq_len(n_runs)
  labs <- .asc_labels(distance_matrix)
  n <- length(labs)
  C <- matrix(0, n, n, dimnames = list(labs, labs))
  memberships <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    m <- ascClusterOnce(distance_matrix, resolution, seed = seeds[i], ...)
    m <- m[labs]
    memberships[[i]] <- m
    C <- C + outer(m, m, FUN = "==")
  }
  C <- C / length(seeds)
  list(memberships = memberships, coassign = C, n_runs = length(seeds))
}

# ---- 3. Retention + change classification -----------------------------------

#' Per-ASC retention / continuity score
#'
#' For each original ASC, measures how much of it stays together in the expanded
#' clustering: \code{retention(k) = max_j |k_old intersect j_new| / |k_old|}.
#'
#' @param old_assign Named vector mapping allele -> original ASC label.
#' @param new_assign Named vector mapping allele -> new cluster label (must cover
#'   all original alleles plus any added ones).
#' @param min_share Minimum fraction of an original ASC that a new cluster must
#'   receive to count as "spanned" (default 0.1).
#'
#' @return A data.frame with columns \code{old_asc}, \code{size},
#'   \code{best_new_cluster}, \code{retention}, \code{n_new_spanned},
#'   \code{n_moved}.
#'
#' @export
ascRetention <- function(old_assign, new_assign, min_share = 0.1) {
  old_assign <- .asc_named_chr(old_assign)
  new_assign <- .asc_named_chr(new_assign)
  ks <- unique(old_assign)
  rows <- lapply(ks, function(k) {
    members <- names(old_assign)[old_assign == k]
    members <- members[members %in% names(new_assign)]
    if (length(members) == 0) {
      return(data.frame(old_asc = k, size = 0L, best_new_cluster = NA_character_,
                        retention = NA_real_, n_new_spanned = 0L, n_moved = 0L,
                        stringsAsFactors = FALSE))
    }
    tab <- table(new_assign[members])
    shares <- as.numeric(tab) / length(members)
    best <- names(tab)[which.max(tab)]
    data.frame(old_asc = k, size = length(members), best_new_cluster = best,
               retention = max(shares),
               n_new_spanned = sum(shares >= min_share),
               n_moved = length(members) - max(tab),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Classify how ASCs change after adding new alleles
#'
#' Labels each original ASC and each new cluster as one of: simple_expansion,
#' new_asc, candidate_split, candidate_merge, or reassignment.
#'
#' @param old_assign Named vector allele -> original ASC label (original alleles).
#' @param new_assign Named vector allele -> new cluster label (original + new).
#' @param new_alleles Character vector of the newly added allele names.
#' @param tau_keep Retention at/above which an ASC is considered "kept together"
#'   (default 0.9).
#' @param p_new Fraction of a new cluster that must be newly added alleles for it
#'   to count as a brand-new ASC (default 0.8).
#' @param p_share Fraction of an original ASC (or of a new cluster's old members)
#'   that a group must hold to count toward a split / merge (default 0.2).
#'
#' @return A data.frame, one row per detected change, with columns
#'   \code{change_type}, \code{old_asc}, \code{new_cluster}, \code{n_old},
#'   \code{n_new_alleles}, \code{retention}, \code{support}, \code{members}.
#'   \code{support} is the number of distinct new alleles associated with the
#'   change (strong-evidence item F).
#'
#' @export
ascPartitionChanges <- function(old_assign, new_assign, new_alleles,
                                tau_keep = 0.9, p_new = 0.8, p_share = 0.2) {
  old_assign <- .asc_named_chr(old_assign)
  new_assign <- .asc_named_chr(new_assign)
  new_alleles <- intersect(new_alleles, names(new_assign))
  changes <- list()

  ## ---- per original ASC: expansion / split / merge / reassignment ----------
  for (k in unique(old_assign)) {
    members <- names(old_assign)[old_assign == k]
    members <- members[members %in% names(new_assign)]
    if (length(members) == 0) next
    tab <- table(new_assign[members])
    shares <- as.numeric(tab) / length(members)
    names(shares) <- names(tab)
    retention <- max(shares)
    span <- names(shares)[shares >= p_share]

    if (length(span) >= 2) {
      new_in_span <- new_alleles[new_assign[new_alleles] %in% span]
      changes[[length(changes) + 1L]] <- data.frame(
        change_type = "candidate_split", old_asc = k,
        new_cluster = paste(span, collapse = "|"), n_old = length(members),
        n_new_alleles = length(new_in_span), retention = retention,
        support = length(new_in_span),
        members = paste(members, collapse = ","), stringsAsFactors = FALSE)
      next
    }

    majority <- names(shares)[which.max(shares)]
    j_members <- names(new_assign)[new_assign == majority]
    j_old <- setdiff(j_members, new_alleles)
    other_old <- unique(old_assign[j_old])
    other_old <- other_old[other_old != k & !is.na(other_old)]
    ## other original ASCs that put a substantial share into this same cluster
    merged_partners <- character(0)
    for (o in other_old) {
      o_members <- names(old_assign)[old_assign == o]
      if (mean(new_assign[o_members] == majority, na.rm = TRUE) >= p_share)
        merged_partners <- c(merged_partners, o)
    }
    new_in_j <- new_alleles[new_assign[new_alleles] == majority]

    if (length(merged_partners) >= 1) {
      changes[[length(changes) + 1L]] <- data.frame(
        change_type = "candidate_merge", old_asc = k,
        new_cluster = majority, n_old = length(members),
        n_new_alleles = length(new_in_j), retention = retention,
        support = length(new_in_j),
        members = paste(sort(unique(c(k, merged_partners))), collapse = "+"),
        stringsAsFactors = FALSE)
    } else if (retention >= tau_keep) {
      changes[[length(changes) + 1L]] <- data.frame(
        change_type = "simple_expansion", old_asc = k,
        new_cluster = majority, n_old = length(members),
        n_new_alleles = length(new_in_j), retention = retention,
        support = length(new_in_j),
        members = paste(members, collapse = ","), stringsAsFactors = FALSE)
    } else {
      moved <- members[new_assign[members] != majority]
      changes[[length(changes) + 1L]] <- data.frame(
        change_type = "reassignment", old_asc = k,
        new_cluster = majority, n_old = length(members),
        n_new_alleles = length(new_in_j), retention = retention,
        support = length(moved),
        members = paste(moved, collapse = ","), stringsAsFactors = FALSE)
    }
  }

  ## ---- per new cluster: brand-new ASCs -------------------------------------
  for (j in unique(new_assign)) {
    members <- names(new_assign)[new_assign == j]
    new_in_j <- intersect(members, new_alleles)
    if (length(members) == 0) next
    frac_new <- length(new_in_j) / length(members)
    if (frac_new >= p_new) {
      changes[[length(changes) + 1L]] <- data.frame(
        change_type = "new_asc", old_asc = NA_character_, new_cluster = j,
        n_old = length(members) - length(new_in_j),
        n_new_alleles = length(new_in_j), retention = NA_real_,
        support = length(new_in_j),
        members = paste(new_in_j, collapse = ","), stringsAsFactors = FALSE)
    }
  }

  if (length(changes) == 0) {
    return(data.frame(change_type = character(0), old_asc = character(0),
                      new_cluster = character(0), n_old = integer(0),
                      n_new_alleles = integer(0), retention = numeric(0),
                      support = integer(0), members = character(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, changes)
}

# ---- 4. Strong-evidence tests -----------------------------------------------

#' Within- vs between-group similarity (evidence C)
#'
#' @param similarity A square similarity matrix (labelled).
#' @param groups A list of character vectors of allele names (the proposed
#'   daughter groups).
#'
#' @return A list with \code{within}, \code{between}, and \code{margin}
#'   (\code{within - between}).
#'
#' @export
ascWithinBetween <- function(similarity, groups) {
  groups <- lapply(groups, function(g) g[g %in% rownames(similarity)])
  groups <- groups[vapply(groups, length, 1L) > 0]
  within_vals <- c(); between_vals <- c()
  for (a in seq_along(groups)) {
    ga <- groups[[a]]
    if (length(ga) >= 2) {
      sub <- similarity[ga, ga, drop = FALSE]
      within_vals <- c(within_vals, sub[upper.tri(sub)])
    }
    for (b in seq_along(groups)) {
      if (b <= a) next
      between_vals <- c(between_vals, as.vector(similarity[ga, groups[[b]], drop = FALSE]))
    }
  }
  w <- if (length(within_vals)) mean(within_vals) else NA_real_
  bt <- if (length(between_vals)) mean(between_vals) else NA_real_
  list(within = w, between = bt, margin = w - bt)
}

#' Whether proposed groups separate under a membership vector
#' @keywords internal
.asc_groups_separate <- function(membership, groups) {
  doms <- vapply(groups, function(g) {
    g <- g[g %in% names(membership)]
    if (length(g) == 0) return(NA_character_)
    tb <- table(membership[g])
    names(tb)[which.max(tb)]
  }, character(1))
  doms <- doms[!is.na(doms)]
  length(doms) >= 2 && length(unique(doms)) == length(doms)
}

#' Leiden stability of a proposed split (evidence A)
#'
#' Fraction of seeded Leiden runs in which the proposed groups land in distinct
#' communities.
#'
#' @param distance_matrix A square distance matrix (or \code{dist}).
#' @param resolution Leiden CPM resolution.
#' @param groups A list of character vectors of allele names.
#' @param n_runs Number of seeded runs (default 100).
#' @param ... Passed to \code{ascLeidenStability}.
#'
#' @return A list with \code{separation} (fraction of runs the groups separate)
#'   and \code{n_runs}.
#'
#' @export
ascEvidenceLeiden <- function(distance_matrix, resolution, groups,
                              n_runs = 100, ...) {
  stab <- ascLeidenStability(distance_matrix, resolution, n_runs = n_runs, ...)
  sep <- vapply(stab$memberships, .asc_groups_separate, logical(1), groups = groups)
  list(separation = mean(sep), n_runs = stab$n_runs)
}

#' Resolution stability of a proposed split (evidence B)
#'
#' Fraction of neighbouring resolutions at which the proposed groups separate.
#'
#' @param distance_matrix A square distance matrix (or \code{dist}).
#' @param resolutions Numeric vector of resolutions to test.
#' @param groups A list of character vectors of allele names.
#' @param seed Seed for the single run at each resolution (default 1).
#' @param ... Passed to \code{ascClusterOnce}.
#'
#' @return A list with \code{persistence} (fraction of resolutions the groups
#'   separate) and \code{detail} (per-resolution logical).
#'
#' @export
ascEvidenceResolution <- function(distance_matrix, resolutions, groups,
                                  seed = 1, ...) {
  sep <- vapply(resolutions, function(r) {
    m <- ascClusterOnce(distance_matrix, r, seed = seed, ...)
    .asc_groups_separate(m, groups)
  }, logical(1))
  list(persistence = mean(sep), detail = stats::setNames(sep, resolutions))
}

#' Silhouette improvement from a split (evidence D)
#'
#' Mean silhouette width of the new partition minus that of the partition in
#' which the original ASC is kept intact, optionally on a subset of alleles.
#'
#' @param distance_matrix A square distance matrix (or \code{dist}).
#' @param partition_new Named membership vector (split applied).
#' @param partition_old_intact Named membership vector (original ASC kept whole).
#' @param subset Optional character vector of alleles to restrict the comparison.
#'
#' @return A list with \code{silhouette_new}, \code{silhouette_old}, and
#'   \code{delta}.
#'
#' @export
ascSilhouetteDelta <- function(distance_matrix, partition_new,
                               partition_old_intact, subset = NULL) {
  M <- if (inherits(distance_matrix, "dist")) as.matrix(distance_matrix) else distance_matrix
  labs <- rownames(M)
  if (is.null(subset)) subset <- labs
  subset <- intersect(subset, labs)
  subset <- intersect(subset, names(partition_new))
  subset <- intersect(subset, names(partition_old_intact))
  sil <- function(part) {
    p <- as.integer(factor(part[subset]))
    if (length(unique(p)) < 2 || length(unique(p)) >= length(p)) return(NA_real_)
    d <- stats::as.dist(M[subset, subset, drop = FALSE])
    si <- cluster::silhouette(p, d)
    if (any(is.na(si))) NA_real_ else mean(si[, 3])
  }
  s_new <- sil(partition_new)
  s_old <- sil(partition_old_intact)
  list(silhouette_new = s_new, silhouette_old = s_old, delta = s_new - s_old)
}

#' Subsampling robustness of a proposed split (evidence E)
#'
#' Repeatedly subsamples the newly added alleles, re-clusters the reference plus
#' the subsample, and reports how often the proposed groups still separate.
#'
#' @param full_distance_matrix A square distance matrix over reference + all new
#'   alleles.
#' @param reference_alleles Character vector of the reference allele names.
#' @param new_alleles Character vector of the newly added allele names.
#' @param groups A list of character vectors of allele names (the proposed split).
#' @param resolution Leiden CPM resolution.
#' @param frac Fraction of new alleles kept per subsample (default 0.8).
#' @param n Number of subsamples (default 20).
#' @param seed_base Base seed; run \code{i} uses \code{seed_base + i}.
#' @param ... Passed to \code{ascClusterOnce}.
#'
#' @return A list with \code{stability} (fraction of subsamples the groups
#'   separate) and \code{n}.
#'
#' @export
ascSubsampleStability <- function(full_distance_matrix, reference_alleles,
                                  new_alleles, groups, resolution,
                                  frac = 0.8, n = 20, seed_base = 1000, ...) {
  M <- if (inherits(full_distance_matrix, "dist")) as.matrix(full_distance_matrix) else full_distance_matrix
  keepn <- max(1L, floor(length(new_alleles) * frac))
  sep <- logical(n)
  for (i in seq_len(n)) {
    set.seed(seed_base + i)
    sampled <- sample(new_alleles, keepn)
    keep <- c(reference_alleles, sampled)
    sub <- M[keep, keep, drop = FALSE]
    m <- ascClusterOnce(sub, resolution, seed = seed_base + i, ...)
    g_present <- lapply(groups, function(g) g[g %in% keep])
    sep[i] <- .asc_groups_separate(m, g_present)
  }
  list(stability = mean(sep), n = n)
}

# ---- 5. Confidence + reconciliation (hysteresis) ----------------------------

#' Classify each change High / Intermediate / Low confidence
#'
#' Applies the asymmetric hysteresis rule: assigning new alleles and creating
#' brand-new ASCs are easy (auto High); splitting, merging, or reassigning
#' established ASCs requires reproducible evidence.
#'
#' @param changes A data.frame from \code{ascPartitionChanges}.
#' @param evidence Optional data.frame joined to \code{changes} by row order with
#'   any of the columns \code{leiden_stability}, \code{resolution_stability},
#'   \code{margin}, \code{silhouette_delta}, \code{subsample_stability}. Missing
#'   columns are treated as failing.
#' @param thresholds Named list of evidence thresholds (see defaults).
#'
#' @return \code{changes} with added \code{n_pass} and \code{confidence} columns.
#'
#' @export
ascClassifyConfidence <- function(changes, evidence = NULL,
                                  thresholds = list(leiden = 0.9,
                                                    resolution = 0.6,
                                                    margin = 0,
                                                    silhouette = 0,
                                                    subsample = 0.8,
                                                    support = 2)) {
  n <- nrow(changes)
  if (n == 0) { changes$n_pass <- integer(0); changes$confidence <- character(0); return(changes) }
  get_col <- function(nm) if (!is.null(evidence) && nm %in% names(evidence)) evidence[[nm]] else rep(NA_real_, n)
  leiden <- get_col("leiden_stability")
  reso   <- get_col("resolution_stability")
  margin <- get_col("margin")
  sildel <- get_col("silhouette_delta")
  subs   <- get_col("subsample_stability")
  support <- changes$support

  pass_mat <- cbind(
    leiden  = !is.na(leiden) & leiden  >= thresholds$leiden,
    reso    = !is.na(reso)   & reso    >= thresholds$resolution,
    margin  = !is.na(margin) & margin  >  thresholds$margin,
    sil     = !is.na(sildel) & sildel  >  thresholds$silhouette,
    subs    = !is.na(subs)   & subs    >= thresholds$subsample,
    support = !is.na(support) & support >= thresholds$support
  )
  n_pass <- rowSums(pass_mat)
  ncrit <- ncol(pass_mat)

  conf <- character(n)
  for (i in seq_len(n)) {
    type <- changes$change_type[i]
    if (type %in% c("simple_expansion", "new_asc")) {
      conf[i] <- "High"                       # easy: absorb / create
    } else if (n_pass[i] == ncrit) {
      conf[i] <- "High"                        # all criteria met
    } else if (n_pass[i] >= ceiling(ncrit / 2)) {
      conf[i] <- "Intermediate"
    } else {
      conf[i] <- "Low"
    }
  }
  changes$n_pass <- n_pass
  changes$confidence <- conf
  changes
}

#' Reconcile old and new ASC assignments with hysteresis
#'
#' Produces a reconciled allele -> ASC assignment that absorbs new alleles into
#' their best original ASC, creates brand-new ASCs, but only applies changes to
#' established ASCs (splits, merges, reassignments) where confidence is High.
#' Intermediate changes are retained as the old assignment and flagged for review.
#'
#' @param old_assign Named vector allele -> original ASC label (original alleles).
#' @param new_assign Named vector allele -> new cluster label (original + new).
#' @param changes A data.frame from \code{ascClassifyConfidence} (must have a
#'   \code{confidence} column).
#' @param new_alleles Character vector of newly added allele names.
#'
#' @return A data.frame with columns \code{allele}, \code{old_asc},
#'   \code{new_cluster}, \code{reconciled_asc}, \code{change_type},
#'   \code{confidence}, \code{review_flag}.
#'
#' @export
ascReconcile <- function(old_assign, new_assign, changes, new_alleles) {
  old_assign <- .asc_named_chr(old_assign)
  new_assign <- .asc_named_chr(new_assign)
  new_alleles <- intersect(new_alleles, names(new_assign))
  all_alleles <- names(new_assign)

  reconciled <- stats::setNames(rep(NA_character_, length(all_alleles)), all_alleles)
  change_type <- stats::setNames(rep(NA_character_, length(all_alleles)), all_alleles)
  confidence <- stats::setNames(rep(NA_character_, length(all_alleles)), all_alleles)
  review <- stats::setNames(rep(FALSE, length(all_alleles)), all_alleles)

  ## start: original alleles keep their old ASC (hysteresis default)
  reconciled[names(old_assign)] <- old_assign

  ## dominant original ASC of a new cluster (for absorbing new alleles)
  dom_old <- function(j) {
    members <- names(new_assign)[new_assign == j]
    j_old <- setdiff(members, new_alleles)
    if (length(j_old) == 0) return(NA_character_)
    tb <- table(old_assign[j_old]); names(tb)[which.max(tb)]
  }

  ## which new clusters are High-confidence brand-new ASCs
  new_asc_clusters <- changes$new_cluster[changes$change_type == "new_asc" &
                                            changes$confidence == "High"]

  ## absorb each new allele
  for (a in new_alleles) {
    j <- new_assign[a]
    if (j %in% new_asc_clusters) {
      reconciled[a] <- paste0("NEW_", j)
      change_type[a] <- "new_asc"; confidence[a] <- "High"
    } else {
      d <- dom_old(j)
      reconciled[a] <- if (is.na(d)) paste0("NEW_", j) else d
      change_type[a] <- if (is.na(d)) "new_asc" else "simple_expansion"
      confidence[a] <- "High"
    }
  }

  ## apply established-ASC changes only where High
  for (i in seq_len(nrow(changes))) {
    ct <- changes$change_type[i]; cf <- changes$confidence[i]
    if (ct %in% c("simple_expansion", "new_asc")) next
    k <- changes$old_asc[i]
    k_members <- names(old_assign)[old_assign == k]
    if (cf == "High") {
      if (ct == "candidate_split") {
        ## relabel each daughter by the new cluster it fell into
        for (a in k_members) {
          reconciled[a] <- paste0(k, ".", new_assign[a])
          change_type[a] <- ct; confidence[a] <- cf
        }
      } else if (ct == "candidate_merge") {
        merged_label <- paste0("MERGE_", changes$members[i])
        partners <- strsplit(changes$members[i], "\\+")[[1]]
        for (a in names(old_assign)[old_assign %in% partners]) {
          reconciled[a] <- merged_label
          change_type[a] <- ct; confidence[a] <- cf
        }
      } else if (ct == "reassignment") {
        moved <- strsplit(changes$members[i], ",")[[1]]
        tgt <- dom_old(changes$new_cluster[i])
        for (a in moved) {
          if (a %in% names(reconciled) && !is.na(tgt)) reconciled[a] <- tgt
          change_type[a] <- ct; confidence[a] <- cf
        }
      }
    } else {
      ## Intermediate / Low: keep old assignment, flag Intermediate for review
      for (a in k_members) {
        if (is.na(change_type[a])) { change_type[a] <- ct; confidence[a] <- cf }
        if (cf == "Intermediate") review[a] <- TRUE
      }
    }
  }

  data.frame(
    allele = all_alleles,
    old_asc = ifelse(all_alleles %in% names(old_assign), old_assign[all_alleles], NA_character_),
    new_cluster = new_assign[all_alleles],
    reconciled_asc = reconciled[all_alleles],
    change_type = change_type[all_alleles],
    confidence = confidence[all_alleles],
    review_flag = review[all_alleles],
    row.names = NULL, stringsAsFactors = FALSE)
}

# ---- 6. Driver + diagnostics ------------------------------------------------

#' Run the incremental-ASC experiment end to end
#'
#' Starts from a reference and its ASC assignments, adds a batch of new sequences,
#' re-runs Leiden (CPM) with the fixed, dataset-independent transform on the
#' expanded set, and produces retention / change / evidence / confidence
#' diagnostics plus a reconciled ASC table.
#'
#' @param reference Named character vector of reference germline sequences
#'   (IMGT-gapped), or a \code{GermlineCluster} whose \code{alleleClusterTable}
#'   provides the baseline assignment.
#' @param new_seqs Named character vector of new germline sequences to add.
#' @param baseline_assign Optional named vector allele -> original ASC. If
#'   \code{NULL}, it is computed by clustering \code{reference}.
#' @param distance_method Distance metric passed to \code{igDistance} (default
#'   \code{"decipher"}).
#' @param trim_3prime_side Fixed trim position for a stable distance frame
#'   (default 318; \code{NULL} to disable).
#' @param resolution Leiden CPM resolution; if \code{NULL} it is chosen by a
#'   silhouette sweep over \code{resolution_grid}.
#' @param resolution_grid Resolutions swept when \code{resolution} is \code{NULL}.
#' @param sim_method Fixed similarity transform, \code{"linear"} (\code{S = 1 - d})
#'   or \code{"log"} (see \code{.asc_fixed_similarity}). \code{"log"} spreads
#'   separation across a gentler resolution band and is preferable for larger,
#'   more diverse references.
#' @param floor Constant floor for the \code{"log"} transform (default 1e-3).
#' @param n_leiden Seeded Leiden repetitions for the stability evidence (default
#'   100).
#' @param n_subsample Subsamples for the robustness evidence (default 20).
#' @param quiet Suppress progress messages (default TRUE).
#'
#' @return An \code{IncrementalASCResult} S3 object (a list) with elements
#'   \code{old_assign}, \code{new_assign}, \code{new_alleles}, \code{resolution},
#'   \code{retention}, \code{changes} (with confidence), \code{reconciled},
#'   \code{diagnostics} (a list of tables), and \code{distance_matrix}.
#'
#' @examples
#' \donttest{
#' data(HVGERM)
#' ref <- HVGERM[1:60]
#' new <- HVGERM[61:80]
#' res <- runIncrementalASCExperiment(ref, new, n_leiden = 20, n_subsample = 5)
#' res$diagnostics$summary_counts
#' }
#'
#' @export
runIncrementalASCExperiment <- function(reference, new_seqs,
                                        baseline_assign = NULL,
                                        distance_method = "decipher",
                                        trim_3prime_side = 318,
                                        resolution = NULL,
                                        resolution_grid = seq(0.02, 0.30, by = 0.02),
                                        sim_method = "linear", floor = 1e-3,
                                        n_leiden = 100, n_subsample = 20,
                                        quiet = TRUE) {
  say <- function(...) if (!quiet) message(...)

  ## ---- baseline assignment -------------------------------------------------
  if (inherits(reference, "GermlineCluster")) {
    tab <- reference$alleleClusterTable
    baseline_assign <- stats::setNames(as.character(tab$allele_cluster), tab$iuis_allele)
    reference <- reference$germlineSet
  }
  ref_names <- names(reference)
  new_names <- names(new_seqs)
  new_names <- setdiff(new_names, ref_names)          # only genuinely new alleles
  new_seqs <- new_seqs[new_names]

  ## ---- fixed distance frame over the union --------------------------------
  prep <- function(x) if (!is.null(trim_3prime_side)) substr(x, 1, trim_3prime_side) else x
  union_seqs <- c(prep(reference), prep(new_seqs))
  say("Computing distances over ", length(union_seqs), " alleles ...")
  full_dist <- as.matrix(igDistance(union_seqs, method = distance_method))

  if (is.null(baseline_assign)) {
    say("Clustering baseline reference ...")
    ref_dist <- full_dist[ref_names, ref_names, drop = FALSE]
    sweep0 <- .asc_resolution_sweep(ref_dist, resolution_grid,
                                    sim_method = sim_method, floor = floor)
    baseline_assign <- as.character(sweep0$best_partition)
    names(baseline_assign) <- names(sweep0$best_partition)
    if (is.null(resolution)) resolution <- sweep0$best_resolution
  }
  old_assign <- .asc_named_chr(baseline_assign)

  ## ---- resolution for the expanded set ------------------------------------
  if (is.null(resolution)) {
    say("Choosing resolution on the expanded set ...")
    sweep1 <- .asc_resolution_sweep(full_dist, resolution_grid,
                                    sim_method = sim_method, floor = floor)
    resolution <- sweep1$best_resolution
  }

  ## ---- new partition on the expanded set ----------------------------------
  new_assign <- as.character(ascClusterOnce(full_dist, resolution, seed = 1,
                                            sim_method = sim_method, floor = floor))
  names(new_assign) <- rownames(full_dist)

  ## ---- retention + changes ------------------------------------------------
  retention <- ascRetention(old_assign, new_assign)
  changes <- ascPartitionChanges(old_assign, new_assign, new_names)

  ## ---- evidence for changes touching established ASCs ---------------------
  S <- .asc_fixed_similarity(full_dist, method = sim_method, floor = floor)
  ev_rows <- vector("list", nrow(changes))
  for (i in seq_len(nrow(changes))) {
    ct <- changes$change_type[i]
    if (ct %in% c("simple_expansion", "new_asc")) {
      ev_rows[[i]] <- data.frame(leiden_stability = NA_real_, resolution_stability = NA_real_,
                                 within = NA_real_, between = NA_real_, margin = NA_real_,
                                 silhouette_delta = NA_real_, subsample_stability = NA_real_)
      next
    }
    ## groups = the split of the old ASC's members by their new cluster
    k <- changes$old_asc[i]
    k_members <- names(old_assign)[old_assign == k]
    grp <- split(k_members, new_assign[k_members])
    grp <- grp[vapply(grp, length, 1L) > 0]
    if (length(grp) < 2) grp <- list(k_members, setdiff(new_names, k_members))
    res_grid <- sort(unique(pmax(resolution * c(0.8, 0.9, 1.0, 1.1, 1.2), 1e-3)))
    lei <- ascEvidenceLeiden(full_dist, resolution, grp, n_runs = n_leiden,
                             sim_method = sim_method, floor = floor)
    rez <- ascEvidenceResolution(full_dist, res_grid, grp,
                                 sim_method = sim_method, floor = floor)
    wb  <- ascWithinBetween(S, grp)
    ## counterfactual: keep the original ASC intact (one label), compare global
    ## silhouette so both partitions have >= 2 clusters.
    old_intact <- new_assign; old_intact[k_members] <- paste0("INTACT_", k)
    sil <- ascSilhouetteDelta(full_dist, new_assign, old_intact)
    sub <- ascSubsampleStability(full_dist, ref_names, new_names, grp, resolution,
                                 n = n_subsample, sim_method = sim_method, floor = floor)
    ev_rows[[i]] <- data.frame(leiden_stability = lei$separation,
                               resolution_stability = rez$persistence,
                               within = wb$within, between = wb$between, margin = wb$margin,
                               silhouette_delta = sil$delta,
                               subsample_stability = sub$stability)
  }
  evidence <- if (length(ev_rows)) do.call(rbind, ev_rows) else NULL
  changes <- ascClassifyConfidence(changes, evidence)
  if (!is.null(evidence)) changes <- cbind(changes, evidence[, c("within", "between",
                                           "silhouette_delta", "leiden_stability",
                                           "resolution_stability", "subsample_stability")])

  ## ---- reconciliation ------------------------------------------------------
  reconciled <- ascReconcile(old_assign, new_assign, changes, new_names)

  ## ---- diagnostics ---------------------------------------------------------
  n_old_changed <- sum(reconciled$allele %in% names(old_assign) &
                       reconciled$reconciled_asc != reconciled$old_asc, na.rm = TRUE)
  summary_counts <- data.frame(
    n_reference = length(ref_names),
    n_new = length(new_names),
    n_old_changed = n_old_changed,
    n_new_asc = sum(changes$change_type == "new_asc"),
    n_candidate_split = sum(changes$change_type == "candidate_split"),
    n_candidate_merge = sum(changes$change_type == "candidate_merge"),
    n_reassignment = sum(changes$change_type == "reassignment"),
    stringsAsFactors = FALSE)

  diagnostics <- list(
    summary_counts = summary_counts,
    retention = retention,
    changes = changes,
    reconciled = reconciled)

  structure(list(
    old_assign = old_assign, new_assign = new_assign, new_alleles = new_names,
    resolution = resolution, retention = retention, changes = changes,
    reconciled = reconciled, diagnostics = diagnostics,
    distance_matrix = full_dist),
    class = "IncrementalASCResult")
}

#' @export
print.IncrementalASCResult <- function(x, ...) {
  cat("IncrementalASCResult\n")
  cat("  reference alleles :", length(x$old_assign), "\n")
  cat("  new alleles       :", length(x$new_alleles), "\n")
  cat("  resolution        :", signif(x$resolution, 3), "\n")
  sc <- x$diagnostics$summary_counts
  cat("  old alleles changed ASC:", sc$n_old_changed, "\n")
  cat("  new ASCs / splits / merges / reassignments:",
      sc$n_new_asc, "/", sc$n_candidate_split, "/",
      sc$n_candidate_merge, "/", sc$n_reassignment, "\n")
  if (nrow(x$changes)) {
    cat("  confidence of changes:\n")
    print(table(x$changes$change_type, x$changes$confidence))
  }
  invisible(x)
}

# ---- helpers ----------------------------------------------------------------

#' Coerce a membership object to a named character vector
#' @keywords internal
.asc_named_chr <- function(x) {
  if (is.null(names(x))) stop("assignment vectors must be named by allele")
  stats::setNames(as.character(x), names(x))
}
