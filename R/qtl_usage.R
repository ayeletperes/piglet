# Gene-usage phenotype and the two in-memory drivers.
#
# Moved from husa_manuscript/lib_manuscript/analysis/igqtl.R (build_usage,
# usage_one_locus, step_pairing). The drivers take matrices and tables rather than paths:
# reading a genotype matrix, a gene BED and a repertoire, and recording which pinned
# release each came from, stays with whatever runs the analysis.

#' Per-subject usage of each gene group
#'
#' A group's share of its own segment's repertoire in each subject, with a pseudocount so
#' the logit is finite. This is the phenotype the usage scan explains.
#'
#' Zero counts are kept deliberately. A variant that deletes a gene drives its usage to
#' zero, and that is the canonical signal in this family, so dropping zeros would drop the
#' effect. Every subject gets a row for every retained group, whether or not they used it.
#'
#' Groups seen in fewer than \code{min_subject_fraction} of subjects are dropped before
#' the shares are taken, and the count kept is reported. \code{n_asc} is how many groups
#' the share was taken over and is part of the pseudocount, so \code{usage} cannot be
#' recomputed from \code{count} and \code{total} without it.
#'
#' @param data A \code{data.frame} of rearrangements, one row per sequence.
#' @param segments Named character vector of the columns holding each segment's gene
#'   call, e.g. \code{c(V = "v_gene_iuis", D = "d_gene_iuis", J = "j_gene_iuis")}. Names
#'   become the \code{segment} column.
#' @param subject Column naming the subject. Default \code{"subject"}.
#' @param min_subject_fraction Minimum fraction of subjects a group must appear in to be
#'   retained. Default 0.1.
#' @param pseudocount Added to the count, and \code{pseudocount * n_asc} to the total.
#'   Default 0.5.
#'
#' @return A \code{data.table} with \code{subject}, \code{asc}, \code{count},
#'   \code{total}, \code{segment}, \code{n_asc}, \code{usage} and \code{logit_usage}.
#' @seealso \code{\link{runGeneUsageQTL}}
#' @examples
#' set.seed(1)
#' rep_dt <- data.frame(
#'   subject = rep(sprintf("s%02d", 1:70), each = 20),
#'   v_gene  = sample(c("V1-2", "V3-23", "V4-34"), 1400, replace = TRUE),
#'   j_gene  = sample(c("J4", "J6"), 1400, replace = TRUE))
#'
#' usage <- ascUsagePhenotype(rep_dt, segments = c(V = "v_gene", J = "j_gene"))
#' head(usage)
#'
#' # Zero counts are kept on purpose: a variant that deletes a gene drives its
#' # usage to zero, and that is the signal rather than missing data.
#' sum(usage$count == 0)
#' @export
ascUsagePhenotype <- function(data, segments, subject = "subject",
                              min_subject_fraction = 0.1, pseudocount = 0.5) {
  x <- data.table::as.data.table(data)
  missing <- setdiff(c(subject, unname(segments)), names(x))
  if (length(missing)) {
    stop("`data` is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(names(segments)) || any(!nzchar(names(segments)))) {
    stop("`segments` must be named, e.g. c(V = \"v_gene_iuis\").", call. = FALSE)
  }

  out <- lapply(names(segments), function(seg) {
    calls <- x[[segments[[seg]]]]
    keep <- !is.na(calls) & nzchar(calls)
    if (!any(keep)) return(NULL)
    d <- data.table::data.table(subject = as.character(x[[subject]])[keep], asc = calls[keep])

    seen <- d[, list(n_subjects = data.table::uniqueN(get("subject"))), by = "asc"]
    n_sub <- data.table::uniqueN(d$subject)
    retained <- seen[get("n_subjects") >= min_subject_fraction * n_sub][["asc"]]
    message(sprintf("ascUsagePhenotype [%s]: %d of %d groups seen in >= %.0f%% of %d subjects.",
                    seg, length(retained), nrow(seen), 100 * min_subject_fraction, n_sub))
    if (!length(retained)) return(NULL)

    counts <- d[, list(count = .N), by = c("subject", "asc")]
    totals <- counts[, list(total = sum(get("count"))), by = "subject"]
    full <- data.table::CJ(subject = unique(d$subject), asc = retained, unique = TRUE)
    full <- merge(full, counts, by = c("subject", "asc"), all.x = TRUE)
    full[is.na(full$count), "count" := 0L]
    full <- merge(full, totals, by = "subject")
    full[, c("segment", "n_asc") := list(seg, length(retained))]
    full[]
  })

  usage <- data.table::rbindlist(Filter(Negate(is.null), out), use.names = TRUE)
  if (!nrow(usage)) return(usage)
  usage[, "usage" := (get("count") + pseudocount) /
          (get("total") + pseudocount * get("n_asc"))]
  usage[, "logit_usage" := stats::qlogis(get("usage"))]
  usage[]
}

# Subjects by group, of logit usage, in the dosage matrix's subject order. Groups with no
# variance are dropped: there is nothing to regress.
.qtl_usage_matrix <- function(usage, subjects) {
  wide <- data.table::dcast(usage, subject ~ asc, value.var = "logit_usage")
  m <- as.matrix(wide[, -1L])
  rownames(m) <- wide$subject
  m <- m[subjects, , drop = FALSE]
  m[, apply(m, 2L, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0), drop = FALSE]
}

# Every member midpoint of a slash-separated group label. Members after the first drop the
# leading segment letter. All midpoints, not their mean: a group can merge genes far
# apart (at IGK the proximal and distal copies of a duplicated gene sit about 800 kb
# apart, so their mean is a location no gene occupies).
.qtl_member_positions <- function(label, positions) {
  seg <- substr(label, 1L, 1L)
  members <- sub("^[VDJ]", "", strsplit(label, "/", fixed = TRUE)[[1]])
  positions$mid[positions$segment == seg & positions$short %in% members]
}

#' Gene-usage QTL scan for one locus
#'
#' Runs the whole usage analysis in memory: builds the phenotype, scans every variant
#' against every gene group, sets the significance threshold from the number of
#' independent variants, and clumps the significant hits into leads.
#'
#' The threshold is \code{alpha} over the number of \emph{independent} variants after
#' exact-LD collapse (\code{\link{qtlLDGroups}}), not over the number tested. Those differ
#' by about half in this data, and using the number tested is the more conservative but
#' wrong denominator.
#'
#' Power columns are attached to the leads, where \code{igqtl.R} put them, and the counts
#' behind every step are messaged: a filter that drops variants without saying so is a bug
#' even when it returns cleanly.
#'
#' @inheritParams ascUsagePhenotype
#' @param dosage Numeric matrix of genotype dosages, variants in rows, subjects in
#'   columns. Subjects are matched by name against the phenotype; filter on MAF and
#'   missingness before calling.
#' @param variants A \code{data.frame} describing the rows of \code{dosage}:
#'   \code{variant}, \code{contig}, \code{pos}, \code{maf}. Not derived from the variant
#'   names, which are a convention of one genotype matrix rather than a fact.
#' @param positions Optional gene coordinates for the locus, with columns \code{segment},
#'   \code{short} and \code{mid}, used for \code{asc_position}, \code{asc_span},
#'   \code{n_member}, \code{distance_to_asc}, \code{is_cis} and \code{nearest_gene}. Those
#'   columns are omitted when it is \code{NULL}.
#' @param locus Locus label written onto the output. Optional.
#' @param min_subjects Minimum complete-case subjects for a missingness group to be
#'   scanned. Default 60.
#' @param alpha Family-wise significance before the independent-variant correction.
#'   Default 0.05.
#' @param r2_clump Squared correlation above which a variant is absorbed into a lead.
#'   Default 0.8.
#' @param min_genotype_group Smallest genotype class at or above which a lead is called
#'   well powered. Default 5.
#' @param cis_window Distance within which a lead is called cis. Default 50000.
#'
#' @return A list of \code{data.table}s: \code{phenotype}, \code{associations},
#'   \code{leads}, \code{per_asc} and \code{thresholds}.
#'
#' @section Ceiling:
#' The scan is unweighted OLS on usage fractions, so its standard errors assume constant
#' variance while binomial sampling variance across subjects spans roughly twentyfold. Any
#' R-squared read off it is diluted and is an underestimate. Do not report an effect size
#' from this as exact.
#'
#' @seealso \code{\link{qtlScanUnivariate}}, \code{\link{runPairingQTL}}
#' @examples
#' set.seed(1)
#' n <- 70L
#' subjects <- sprintf("s%02d", seq_len(n))
#' dosage <- rbind(v_hit  = rep(c(0, 1, 2), length.out = n),
#'                 v_null = rep(c(0, 0, 1, 2), length.out = n))
#' colnames(dosage) <- subjects
#'
#' genes <- c("V1-2", "V3-23", "V4-34")
#' rep_dt <- do.call(rbind, lapply(seq_len(n), function(i) {
#'   w <- c(1 + dosage["v_hit", i], 1, 1)   # dosage raises V1-2 usage
#'   data.frame(subject = subjects[i],
#'              v_gene = sample(genes, 40, replace = TRUE, prob = w / sum(w)))
#' }))
#'
#' variants <- data.frame(variant = rownames(dosage), contig = "igh",
#'                        pos = c(1000L, 50000L), maf = rowMeans(dosage) / 2)
#'
#' res <- runGeneUsageQTL(rep_dt, dosage, variants,
#'                        segments = c(V = "v_gene"), min_subjects = 60)
#' res$associations
#' res$thresholds
#' @export
runGeneUsageQTL <- function(data, dosage, variants, segments, positions = NULL,
                            locus = NA_character_, subject = "subject",
                            min_subject_fraction = 0.1, pseudocount = 0.5,
                            min_subjects = 60L, alpha = 0.05, r2_clump = 0.8,
                            min_genotype_group = 5L, cis_window = 50000) {
  usage <- ascUsagePhenotype(data, segments, subject, min_subject_fraction, pseudocount)
  if (!nrow(usage)) stop("No gene group survived `min_subject_fraction`.", call. = FALSE)

  info <- data.table::as.data.table(variants)
  if (!all(c("variant", "contig", "pos", "maf") %in% names(info))) {
    stop("`variants` needs columns variant, contig, pos, maf.", call. = FALSE)
  }
  subjects <- colnames(dosage)
  pheno <- .qtl_usage_matrix(usage, subjects)
  message(sprintf("runGeneUsageQTL: %d subjects, %d groups scanned, %d variants.",
                  nrow(pheno), ncol(pheno), nrow(dosage)))

  ld <- qtlLDGroups(dosage, info$contig[match(rownames(dosage), info$variant)])
  n_indep <- data.table::uniqueN(ld)
  thr <- alpha / n_indep
  message(sprintf("runGeneUsageQTL: %d variants collapse to %d independent; threshold %.3g.",
                  nrow(dosage), n_indep, thr))

  assoc <- qtlScanUnivariate(pheno, dosage, min_subjects)
  data.table::setnames(assoc, "phenotype", "asc")
  assoc <- merge(assoc, info, by = "variant", sort = FALSE)

  asc_meta <- unique(usage[, list(asc = get("asc"), segment = get("segment"))])
  asc_meta[, "locus" := locus]
  members <- NULL
  if (!is.null(positions)) {
    members <- data.table::rbindlist(lapply(asc_meta$asc, function(a) {
      data.table::data.table(asc = a, mid = .qtl_member_positions(a, positions))
    }))
    asc_meta <- merge(asc_meta,
                      members[, list(asc_position = mean(get("mid")),
                                     asc_span = diff(range(get("mid"))),
                                     n_member = .N), by = "asc"],
                      by = "asc", all.x = TRUE)
  }
  assoc <- merge(assoc, asc_meta, by = "asc", sort = FALSE)

  if (!is.null(members)) {
    # Distance to the nearest member gene, not to the mean of the members.
    assoc[, "distance_to_asc" := {
      m <- members$mid[members$asc == get("asc")[1L]]
      if (!length(m)) NA_real_ else vapply(get("pos"), function(p) min(abs(p - m)), 1)
    }, by = "asc"]
  }
  assoc[, c("threshold", "significant") := list(thr, get("p_value") < thr)]
  data.table::setorderv(assoc, "p_value")

  sig <- assoc[assoc$significant %in% TRUE]
  leads <- if (nrow(sig)) {
    data.table::rbindlist(lapply(split(sig, by = "asc"), qtlClump,
                                 dosage = dosage, r2 = r2_clump))
  } else {
    sig
  }
  if (nrow(leads)) {
    mgc <- qtlSmallestGenotypeClass(dosage)
    leads[, c("min_genotype_group", "well_powered") :=
            list(mgc[get("variant")], mgc[get("variant")] >= min_genotype_group)]
    if (!is.null(positions)) {
      leads[, "is_cis" := get("distance_to_asc") <= cis_window]
      leads[, "nearest_gene" := vapply(get("pos"), function(p) {
        positions$gene[order(pmin(abs(positions$start - p), abs(positions$end - p)))][1L]
      }, character(1))]
    }
    data.table::setorderv(leads, "p_value")
  }
  message(sprintf("runGeneUsageQTL: %d associations, %d significant, %d leads over %d groups.",
                  nrow(assoc), nrow(sig), nrow(leads),
                  data.table::uniqueN(leads$asc)))

  per_asc <- merge(asc_meta,
                   assoc[, list(min_p = min(get("p_value")),
                                n_significant = sum(get("significant")),
                                n_significant_groups = data.table::uniqueN(
                                  ld[get("variant")[get("significant")]])),
                         by = "asc"],
                   by = "asc", all.x = TRUE)
  per_asc <- merge(per_asc,
                   usage[, list(median_usage = stats::median(get("count") / get("total")),
                                subjects_nonzero = sum(get("count") > 0L)), by = "asc"],
                   by = "asc", all.x = TRUE)

  sig_ld <- ld[sig$variant]
  list(phenotype = usage,
       associations = assoc[],
       leads = leads[],
       per_asc = per_asc[],
       thresholds = data.table::data.table(
         locus = locus, analysis = "usage", n_subjects = nrow(pheno),
         n_variants = nrow(dosage), n_independent = n_indep, n_asc = nrow(asc_meta),
         threshold = thr, n_significant_variants = data.table::uniqueN(sig$variant),
         n_independent_significant = data.table::uniqueN(sig_ld)))
}

#' Conditional gene-pairing QTL scan for one anchoring
#'
#' Runs one side of the pairing analysis: builds the enrichment table, scans every variant
#' against each anchor's partner distribution, sets the threshold from the number of
#' independent variants, and clumps the hits into leads.
#'
#' One anchoring per call, on purpose. Anchoring on each of the two segments gives two
#' different multivariate tests over the same variants, and a function that ran both and
#' returned them stacked would invite exactly the pooling that double-reports. Call it
#' twice, with the \code{conditional} each side estimates.
#'
#' @inheritParams pairingScan
#' @inheritParams runGeneUsageQTL
#' @param data A \code{data.frame} of rearrangements, one row per sequence.
#' @param anchor,partner Columns holding the gene to condition on and the partner gene.
#'
#' @return A list of \code{data.table}s: \code{pairs}, \code{associations}, \code{leads}
#'   and \code{thresholds}.
#' @seealso \code{\link{pairingTable}}, \code{\link{pairingScan}},
#'   \code{\link{pairingCellTests}}
#' @examples
#' set.seed(1)
#' subjects <- sprintf("s%02d", 1:70)
#' rep_dt <- data.frame(
#'   subject = rep(subjects, each = 20),
#'   d_gene  = sample(c("D1", "D2", "D3"), 1400, replace = TRUE),
#'   j_gene  = sample(c("J1", "J2"), 1400, replace = TRUE))
#'
#' dosage <- matrix(rep(c(0, 1, 2), length.out = 70), nrow = 1,
#'                  dimnames = list("v1", subjects))
#' variants <- data.frame(variant = "v1", contig = "igh", pos = 1000L, maf = 0.33)
#'
#' # One anchoring per call. Running both and stacking them would double-report,
#' # so the caller states which conditional each side estimates.
#' jd <- runPairingQTL(rep_dt, dosage, variants,
#'                     anchor = "j_gene", partner = "d_gene",
#'                     conditional = "P(J|D)", min_subjects = 60)
#' jd$associations
#' @export
runPairingQTL <- function(data, dosage, variants, anchor, partner, conditional,
                          subject = "subject", locus = NA_character_, pseudocount = 0.5,
                          min_subjects = 60L, min_complete_fraction = 0.9, alpha = 0.05,
                          r2_clump = 0.8, min_genotype_group = 5L) {
  info <- data.table::as.data.table(variants)
  if (!all(c("variant", "contig", "pos", "maf") %in% names(info))) {
    stop("`variants` needs columns variant, contig, pos, maf.", call. = FALSE)
  }

  pairs <- pairingTable(data, anchor = anchor, partner = partner, subject = subject,
                        pseudocount = pseudocount)

  ld <- qtlLDGroups(dosage, info$contig[match(rownames(dosage), info$variant)])
  n_indep <- data.table::uniqueN(ld)
  thr <- alpha / n_indep
  message(sprintf("runPairingQTL [%s]: %d variants collapse to %d independent; threshold %.3g.",
                  conditional, nrow(dosage), n_indep, thr))

  assoc <- pairingScan(pairs, dosage, conditional, min_subjects = min_subjects,
                       min_complete_fraction = min_complete_fraction,
                       min_genotype_group = min_genotype_group)
  if (!nrow(assoc)) {
    return(list(pairs = pairs, associations = assoc, leads = assoc,
                thresholds = data.table::data.table()))
  }
  assoc <- merge(assoc, info, by = "variant", sort = FALSE)
  assoc[, c("threshold", "significant") := list(thr, get("p_value") < thr)]
  data.table::setorderv(assoc, "p_value")

  sig <- assoc[assoc$significant %in% TRUE]
  leads <- if (nrow(sig)) {
    data.table::rbindlist(lapply(split(sig, by = "anchor_gene"), qtlClump,
                                 dosage = dosage, r2 = r2_clump))
  } else {
    sig
  }
  sig_ld <- ld[sig$variant]
  message(sprintf("runPairingQTL [%s]: %d anchors, %d associations, %d significant over %d independent groups, %d leads.",
                  conditional, data.table::uniqueN(assoc$anchor_gene), nrow(assoc),
                  nrow(sig), data.table::uniqueN(sig_ld), nrow(leads)))

  list(pairs = pairs,
       associations = assoc[],
       leads = leads[],
       thresholds = data.table::data.table(
         locus = locus, analysis = "pairing", grouped_by = anchor,
         conditional = conditional, n_subjects = length(unique(pairs$subject)),
         n_variants = nrow(dosage), n_independent = n_indep,
         n_asc = data.table::uniqueN(assoc$anchor_gene), threshold = thr,
         n_significant_variants = data.table::uniqueN(sig$variant),
         n_independent_significant = data.table::uniqueN(sig_ld)))
}
