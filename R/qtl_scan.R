# Scan primitives for gene-usage QTL analysis.
#
# Moved from husa_manuscript/lib_manuscript/analysis/igqtl.R, where they were developed.
# The arithmetic is unchanged: these functions reproduce that script's published output
# row for row (see the package tests and vignette).

#' Univariate QTL scan of every phenotype on every variant
#'
#' Ordinary least squares of each phenotype column on each variant, for the whole matrix
#' at once. The slope of a simple regression is a ratio of cross products, so the entire
#' scan is two matrix products per missingness group rather than one \code{lm()} call per
#' variant-phenotype pair.
#'
#' Variants are grouped by which subjects have no call, so every fit inside a group uses
#' the same complete-case subject set and \code{n} is exact for each row rather than an
#' upper bound. Groups smaller than \code{min_n} are skipped entirely.
#'
#' A p-value from this scan must not be read without the smallest genotype class behind
#' the fit: the extreme tail is anti-conservative, and a non-significant result means
#' either no effect or no power. Use \code{\link{qtlSmallestGenotypeClass}} on the same
#' dosage matrix, as \code{\link{runGeneUsageQTL}} does. Filtering associations on the
#' p-value alone preferentially removes the underpowered rows, which leaves what survives
#' looking better powered than it is.
#'
#' @param pheno Numeric matrix of phenotypes, subjects in rows (row names are subject
#'   identifiers) and phenotypes in columns (column names name the phenotype). No missing
#'   values; drop or impute before calling.
#' @param dosage Numeric matrix of genotype dosages, variants in rows (row names are
#'   variant identifiers) and subjects in columns, in the same subject order as
#'   \code{pheno}. \code{NA} marks no call.
#' @param min_n Minimum number of complete-case subjects a missingness group needs before
#'   it is scanned. Default 60.
#'
#' @return A \code{data.table} with one row per variant-phenotype pair that produced a
#'   finite p-value: \code{variant}, \code{phenotype}, \code{n} (subjects in the fit),
#'   \code{beta}, \code{se}, \code{t_stat}, \code{p_value}. Effect direction is relative
#'   to the coding of \code{dosage}, which is not necessarily the minor allele, so
#'   \code{beta} must not be described in terms of a named allele without resolving
#'   polarity first.
#'
#' @seealso \code{\link{qtlScanMultivariate}} for a multivariate response,
#'   \code{\link{qtlSmallestGenotypeClass}}, \code{\link{runGeneUsageQTL}}
#' @export
qtlScanUnivariate <- function(pheno, dosage, min_n = 60L) {
  .qtl_check_matrices(pheno, dosage)
  groups <- .qtl_missingness_groups(dosage)

  out <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    rows <- .qtl_group_rows(names(groups)[k], nrow(pheno))
    n <- length(rows)
    if (n < min_n) next

    y <- sweep(pheno[rows, , drop = FALSE], 2L, colMeans(pheno[rows, , drop = FALSE]))
    ss <- colSums(y^2)
    g <- dosage[idx, rows, drop = FALSE]
    g <- g - rowMeans(g)
    den <- rowSums(g^2)

    num <- g %*% y
    beta <- num / den
    se <- sqrt(pmax(sweep(-beta * num, 2L, ss, "+"), 0) / (n - 2) / den)
    t_stat <- beta / se
    p <- 2 * stats::pt(-abs(t_stat), n - 2)

    cells <- which(is.finite(p), arr.ind = TRUE)
    if (!nrow(cells)) next
    out[[k]] <- data.table::data.table(
      variant = rownames(g)[cells[, 1L]], phenotype = colnames(y)[cells[, 2L]],
      n = n, beta = beta[cells], se = se[cells],
      t_stat = t_stat[cells], p_value = p[cells])
  }
  data.table::rbindlist(out, use.names = TRUE)
}

#' Multivariate QTL scan by Pillai's trace
#'
#' Scores every variant against a multivariate response, one MANOVA per variant. For a
#' one degree of freedom predictor Pillai's trace has a closed form,
#' \deqn{V = (u' S^{-1} u) / (g'g), \quad u = Y_c' g_c, \quad S = Y_c' Y_c}
#' and the F approximation is exact, so a whole missingness group is scored by one matrix
#' product instead of a \code{stats::manova()} call per variant. The package tests check
#' this against \code{stats::manova()}.
#'
#' Response columns with zero variance are dropped within each missingness group, and a
#' group is skipped when fewer than two response columns survive or when \code{n} is too
#' small for the covariance to be estimable (\code{n < p + 3}).
#'
#' @inheritParams qtlScanUnivariate
#' @param pheno Numeric matrix of the multivariate response, subjects in rows, response
#'   variables in columns. Complete cases only.
#'
#' @return A \code{data.table} with \code{variant}, \code{n}, \code{n_response} (response
#'   columns retained), \code{pillai}, \code{f_stat} and \code{p_value}.
#'
#' @seealso \code{\link{pairingScan}}, which uses this for the conditional pairing scan
#' @export
qtlScanMultivariate <- function(pheno, dosage, min_n = 60L) {
  .qtl_check_matrices(pheno, dosage)
  groups <- .qtl_missingness_groups(dosage)

  out <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    rows <- .qtl_group_rows(names(groups)[k], nrow(pheno))
    n <- length(rows)
    if (n < min_n) next

    y <- sweep(pheno[rows, , drop = FALSE], 2L, colMeans(pheno[rows, , drop = FALSE]))
    y <- y[, apply(y, 2L, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0), drop = FALSE]
    p <- ncol(y)
    if (p < 2L || n < p + 3L) next
    s_inv <- tryCatch(solve(crossprod(y)), error = function(e) NULL)
    if (is.null(s_inv)) next

    g <- dosage[idx, rows, drop = FALSE]
    g <- g - rowMeans(g)
    den <- rowSums(g^2)
    u <- crossprod(y, t(g))
    v <- colSums(u * (s_inv %*% u)) / den
    v[!is.finite(v) | v <= 0 | v >= 1] <- NA_real_
    f <- ((n - p - 1) / p) * (v / (1 - v))
    out[[k]] <- data.table::data.table(
      variant = rownames(g), n = n, n_response = p, pillai = v, f_stat = f,
      p_value = stats::pf(f, p, n - p - 1, lower.tail = FALSE))
  }
  data.table::rbindlist(out, use.names = TRUE)
}

#' Collapse variants carrying identical genotype information
#'
#' Two variants carry the same information when their standardised dosage vectors match
#' up to sign, so the mirrored case is folded in. This is exact linkage disequilibrium
#' rather than a correlation threshold; use \code{\link{qtlClump}} for the latter.
#'
#' The number of distinct groups is the number of independent tests, which is what a
#' Bonferroni threshold should be set from -- the number of variants tested overstates it,
#' often by half. The contig prefix stops identical patterns at different loci from
#' merging, so positions from different loci are never pooled.
#'
#' @param dosage Numeric matrix of genotype dosages, variants in rows, subjects in
#'   columns. Missing calls are mean-imputed within the variant for the comparison only.
#' @param contig Character vector, one contig per row of \code{dosage}.
#'
#' @return A character vector, one group key per row of \code{dosage}. The key is opaque:
#'   use it to group and count, not to read. Variants with no variance each get their own
#'   group rather than collapsing together.
#' @export
qtlLDGroups <- function(dosage, contig) {
  if (length(contig) != nrow(dosage)) {
    stop("`contig` must have one entry per row of `dosage`: ",
         length(contig), " against ", nrow(dosage), ".", call. = FALSE)
  }
  z <- dosage
  mu <- rowMeans(z, na.rm = TRUE)
  na_idx <- which(is.na(z), arr.ind = TRUE)
  if (nrow(na_idx)) z[na_idx] <- mu[na_idx[, 1L]]
  sdev <- apply(z, 1L, stats::sd)
  flat <- !is.finite(sdev) | sdev == 0
  z <- round((z - mu) / pmax(sdev, .Machine$double.eps), 8L)
  key <- pmin(apply(z, 1L, paste, collapse = ","), apply(-z, 1L, paste, collapse = ","))
  key[flat] <- paste0("const_", seq_len(sum(flat)))
  paste(contig, key, sep = "::")
}

#' Independent lead variants by greedy clumping
#'
#' Takes the strongest signal, absorbs everything correlated with it above \code{r2},
#' and repeats. What comes back is one row per independent signal rather than one row per
#' significant variant, which is the number to quote.
#'
#' @param assoc A \code{data.table} of significant associations carrying at least
#'   \code{variant} and \code{p_value}. Reordered by \code{p_value} in place.
#' @param dosage Numeric matrix of genotype dosages, variants in rows, subjects in
#'   columns, covering every variant in \code{assoc}. Missing calls are mean-imputed
#'   within the variant for the correlation only.
#' @param r2 Squared-correlation threshold above which a variant is absorbed into an
#'   already-taken lead.
#'
#' @return The subset of \code{assoc} that are leads, in increasing p-value order.
#' @export
qtlClump <- function(assoc, dosage, r2 = 0.8) {
  if (!nrow(assoc)) return(assoc)
  data.table::setorderv(assoc, "p_value")
  sub <- dosage[unique(assoc$variant), , drop = FALSE]
  mu <- rowMeans(sub, na.rm = TRUE)
  na_idx <- which(is.na(sub), arr.ind = TRUE)
  if (nrow(na_idx)) sub[na_idx] <- mu[na_idx[, 1L]]
  r2m <- stats::cor(t(sub), use = "pairwise.complete.obs")^2
  taken <- character(0)
  leads <- integer(0)
  for (i in seq_len(nrow(assoc))) {
    v <- assoc$variant[i]
    if (v %in% taken) next
    leads <- c(leads, i)
    taken <- union(taken, colnames(r2m)[is.finite(r2m[v, ]) & r2m[v, ] >= r2])
  }
  assoc[leads]
}

#' Smallest observed genotype class per variant
#'
#' How many subjects sit in the least populated genotype class at each variant. This is
#' what says whether a fit is worth reading: an association driven by two homozygotes is
#' not the same claim as one driven by forty, and the scan's extreme tail is
#' anti-conservative precisely where this number is small.
#'
#' Classes with no subject are not counted, so a variant seen only as 0 and 1 reports the
#' smaller of those two rather than zero. Computed for the whole matrix at once; index
#' the result by variant name rather than calling this per variant.
#'
#' @param dosage Numeric matrix of genotype dosages, variants in rows, subjects in
#'   columns. Dosages are rounded to 0, 1 or 2.
#'
#' @return A named numeric vector, one entry per row of \code{dosage}, \code{NA} where no
#'   class had a subject.
#' @export
qtlSmallestGenotypeClass <- function(dosage) {
  g <- round(dosage)
  n <- vapply(0:2, function(k) rowSums(g == k, na.rm = TRUE), numeric(nrow(g)))
  n[n == 0L] <- NA_real_
  out <- suppressWarnings(apply(n, 1L, min, na.rm = TRUE))
  out[!is.finite(out)] <- NA_real_
  stats::setNames(out, rownames(dosage))
}

# Phenotype means in the lowest and highest observed genotype class, for every variant by
# response pair at once. The per-allele slope is fitted on the model scale; this is the
# raw difference to read it against. An indicator matrix per class turns the whole
# calculation into three matrix products.
.qtl_class_means <- function(pheno, dosage) {
  g <- round(dosage)
  ind <- lapply(0:2, function(k) {
    x <- (g == k) & !is.na(g)
    storage.mode(x) <- "numeric"
    x
  })
  n <- vapply(ind, rowSums, numeric(nrow(g)))
  mean_k <- lapply(seq_along(ind), function(k) (ind[[k]] %*% pheno) / pmax(n[, k], 1))
  present <- n > 0L
  lo <- max.col(present, ties.method = "first")
  hi <- 4L - max.col(present[, 3:1, drop = FALSE], ties.method = "first")
  pick <- function(idx) {
    m <- mean_k[[1L]]
    for (k in 2:3) m[idx == k, ] <- mean_k[[k]][idx == k, , drop = FALSE]
    m
  }
  ml <- pick(lo)
  mh <- pick(hi)
  data.table::data.table(
    variant = rep(rownames(g), times = ncol(pheno)),
    response_gene = rep(colnames(pheno), each = nrow(g)),
    n_low = n[cbind(seq_along(lo), lo)], n_high = n[cbind(seq_along(hi), hi)],
    mean_low = as.vector(ml), mean_high = as.vector(mh),
    delta_mean = as.vector(mh - ml))
}

# Group variants by which subjects are missing, so each group shares one subject subset
# and every fit inside it is a complete-case fit with the same n.
.qtl_missingness_groups <- function(dosage) {
  key <- apply(is.na(dosage), 1L, function(x) paste(which(x), collapse = ","))
  split(seq_len(nrow(dosage)), key)
}

.qtl_group_rows <- function(name, n_total) {
  drop_i <- suppressWarnings(as.integer(strsplit(name, ",", fixed = TRUE)[[1]]))
  drop_i <- drop_i[!is.na(drop_i)]
  if (length(drop_i)) setdiff(seq_len(n_total), drop_i) else seq_len(n_total)
}

# Fail before the arithmetic rather than after: a subject-order mismatch between the two
# matrices produces numbers that look plausible and are wrong.
.qtl_check_matrices <- function(pheno, dosage) {
  if (!is.matrix(pheno) || !is.numeric(pheno)) {
    stop("`pheno` must be a numeric matrix (subjects in rows).", call. = FALSE)
  }
  if (!is.matrix(dosage) || !is.numeric(dosage)) {
    stop("`dosage` must be a numeric matrix (variants in rows, subjects in columns).",
         call. = FALSE)
  }
  if (ncol(dosage) != nrow(pheno)) {
    stop("`dosage` has ", ncol(dosage), " subject columns but `pheno` has ",
         nrow(pheno), " subject rows.", call. = FALSE)
  }
  if (anyNA(pheno)) {
    stop("`pheno` has missing values; drop or impute before scanning.", call. = FALSE)
  }
  if (!is.null(rownames(pheno)) && !is.null(colnames(dosage)) &&
      !identical(rownames(pheno), colnames(dosage))) {
    stop("Subject order differs between `pheno` rows and `dosage` columns.", call. = FALSE)
  }
  invisible(TRUE)
}
