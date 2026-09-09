# Conditional gene-pairing QTL: does a variant change which partner a gene pairs with?
#
# Moved from husa_manuscript/lib_manuscript/analysis/igqtl.R (build_pairs,
# pairing_phenotype, run_pairing_side, step_cell_tests), generalised from the D/J
# hardcoding to any anchor and partner segment.

#' Observed against expected gene pairing, per subject
#'
#' Builds the anchor-by-partner contingency of one subject's repertoire and compares it
#' to what independence of the two marginals would give. This is the phenotype the
#' conditional pairing scan explains.
#'
#' \strong{Both marginals are divided out.} \code{expected} is the product of the two
#' margins over the subject's depth, so \code{enrichment} measures pairing preference
#' rather than how much either gene is used. That distinction is the whole analysis: in
#' the heavy chain the individual-level signal in the raw conditional
#' \code{p_partner_given_anchor} was found to be almost entirely marginal usage, and it
#' went to roughly zero once both marginals were removed. The raw conditionals are
#' returned for inspection, but \code{enrichment} is what should be scanned.
#'
#' Zero cells are created only where both genes occur somewhere in that subject, so a
#' gene the subject does not carry stays a structural absence rather than a measured
#' zero. This matters: a deletion is not the same observation as a pairing the subject
#' never makes.
#'
#' The enrichment cell is symmetric -- swapping which margin conditions leaves the number
#' unchanged -- so to scan the other side, call this again with \code{anchor} and
#' \code{partner} swapped. The omnibus test is \emph{not} symmetric, because grouping the
#' cells by one segment gives each anchor a response vector of a different width, which
#' is a different multivariate test on a different covariance structure. Both sides are
#' therefore separate scans and their results must not be pooled.
#'
#' @section Which combinations are established:
#' IGH D-anchored and J-anchored pairing is the analysis this was developed for and the
#' one that has been checked. IGK/IGL V--J, and IGH V--D and V--J, run through the same
#' code and are \strong{exploratory}: nobody has shown that pairing preference in those
#' segments survives removing the marginals, and the heavy-chain result above is a reason
#' to expect trouble rather than reassurance. Treat their output as a hypothesis.
#'
#' @param data A \code{data.frame} of rearrangements, one row per sequence.
#' @param subject Column naming the subject. Default \code{"subject"}.
#' @param anchor Column holding the gene to condition on.
#' @param partner Column holding the partner gene whose distribution is the response.
#' @param pseudocount Added to both count and expectation before the log ratio, so a zero
#'   cell gives a finite enrichment. Default 0.5.
#'
#' @return A \code{data.table} with one row per subject, anchor gene and partner gene:
#'   \code{subject}, \code{anchor_gene}, \code{partner_gene}, \code{count}, \code{depth},
#'   \code{anchor_total}, \code{partner_total}, \code{expected}, \code{enrichment},
#'   \code{p_anchor}, \code{p_partner}, \code{p_partner_given_anchor},
#'   \code{p_anchor_given_partner}. The \code{anchor} and \code{partner} column names are
#'   kept on the result as attributes.
#'
#' @seealso \code{\link{pairingScan}}, \code{\link{pairingCellTests}}
#' @export
pairingTable <- function(data, anchor, partner, subject = "subject", pseudocount = 0.5) {
  x <- data.table::as.data.table(data)
  missing <- setdiff(c(subject, anchor, partner), names(x))
  if (length(missing)) {
    stop("`data` is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  # Pulled out by name before any data.table expression: `anchor` and `partner` are
  # column names in a repertoire too, and inside `[` a bare name resolves to the column.
  keep <- data.table::data.table(subject = as.character(x[[subject]]),
                                 anchor_gene = x[[anchor]],
                                 partner_gene = x[[partner]])
  n_in <- nrow(keep)
  keep <- keep[!is.na(keep$anchor_gene) & nzchar(keep$anchor_gene) &
                 !is.na(keep$partner_gene) & nzchar(keep$partner_gene)]
  message(sprintf("pairingTable: %d of %d rearrangements carry both a %s and a %s call.",
                  nrow(keep), n_in, anchor, partner))
  if (!nrow(keep)) return(keep[0])

  counts <- keep[, list(count = .N), by = c("subject", "anchor_gene", "partner_gene")]
  grid <- merge(unique(counts[, c("subject", "anchor_gene"), with = FALSE]),
                unique(counts[, c("subject", "partner_gene"), with = FALSE]),
                by = "subject", allow.cartesian = TRUE)
  p <- merge(grid, counts, by = c("subject", "anchor_gene", "partner_gene"), all.x = TRUE)
  p[is.na(p$count), "count" := 0L]
  p[, "depth" := sum(get("count")), by = "subject"]
  p[, "anchor_total" := sum(get("count")), by = c("subject", "anchor_gene")]
  p[, "partner_total" := sum(get("count")), by = c("subject", "partner_gene")]
  p[, "expected" := get("anchor_total") * get("partner_total") / get("depth")]
  p[, "enrichment" := log2((get("count") + pseudocount) / (get("expected") + pseudocount))]
  p[, c("p_anchor", "p_partner", "p_partner_given_anchor", "p_anchor_given_partner") :=
      list(get("anchor_total") / get("depth"), get("partner_total") / get("depth"),
           get("count") / get("anchor_total"), get("count") / get("partner_total"))]

  data.table::setattr(p, "anchor_column", anchor)
  data.table::setattr(p, "partner_column", partner)
  p[]
}

# The response matrix for one anchor: subjects by partner gene, of enrichment. Partner
# genes seen in too few subjects are dropped before complete cases are taken, so one
# rare partner cannot cost the whole anchor its subjects.
.qtl_pairing_phenotype <- function(pairs, anchor, subjects, min_complete_fraction, min_subjects) {
  sub <- pairs[pairs$anchor_gene == anchor]
  w <- data.table::dcast(sub, subject ~ partner_gene, value.var = "enrichment")
  m <- as.matrix(w[, -1L])
  rownames(m) <- w$subject
  m <- m[, colMeans(is.finite(m)) >= min_complete_fraction, drop = FALSE]
  m <- m[rownames(m) %in% subjects & stats::complete.cases(m), , drop = FALSE]
  if (ncol(m) < 2L || nrow(m) < min_subjects) NULL else m
}

#' Scan variants against a gene's partner distribution
#'
#' One anchor gene at a time, tests every variant against the whole vector of partner
#' enrichments by MANOVA (Pillai's trace, via \code{\link{qtlScanMultivariate}}). A
#' variant that merely raises or lowers how much the anchor is used shifts every partner
#' in the same direction and produces no signal; one that repartitions the anchor's
#' partners does.
#'
#' \code{conditional} has no default and must be stated. The two anchorings are separate
#' scans over the same variants, so a result that does not say which one it is cannot be
#' pooled, counted or plotted without double-reporting. The direction is easy to invert
#' by mistake and has been inverted once: at fixed partner,
#' \deqn{enrichment = \log_2 P(\mathrm{partner} \mid \mathrm{anchor}) - \log_2 P(\mathrm{partner})}
#' and the second term is constant across the anchors, so the anchored response vector is
#' \eqn{\log_2 P(\mathrm{partner} \mid \mathrm{anchor})} shifted by a constant. Anchoring
#' on J therefore gives the P(J|D) scan, and anchoring on D gives P(D|J).
#'
#' \code{min_genotype_group} is counted on each anchor's own subject set, since complete
#' cases differ between anchors and a variant can be well powered for one anchor and not
#' another. It is returned on every row, not only on the leads: a p-value from this scan
#' cannot be read without it.
#'
#' @param pairs A \code{data.table} from \code{\link{pairingTable}}.
#' @param dosage Numeric matrix of genotype dosages, variants in rows, subjects in
#'   columns.
#' @param conditional A string naming what this scan estimates, e.g. \code{"P(J|D)"}.
#'   Required.
#' @param anchors Anchor genes to scan. Default all of them.
#' @param min_subjects Minimum complete-case subjects for an anchor to be scanned, and
#'   for a missingness group inside it. Default 60.
#' @param min_complete_fraction Minimum fraction of subjects in which a partner gene must
#'   be observed for it to stay in the response. Default 0.9.
#' @param min_genotype_group Smallest genotype class at or above which a fit is called
#'   well powered. Default 5.
#'
#' @return A \code{data.table} with \code{conditional}, \code{variant},
#'   \code{anchor_gene}, \code{n}, \code{n_response}, \code{pillai}, \code{f_stat},
#'   \code{p_value}, \code{min_genotype_group} and \code{well_powered}. No significance
#'   threshold is applied: set it from the number of independent variants
#'   (\code{\link{qtlLDGroups}}), not the number tested.
#'
#' @seealso \code{\link{pairingTable}}, \code{\link{pairingCellTests}}
#' @export
pairingScan <- function(pairs, dosage, conditional, anchors = NULL,
                        min_subjects = 60L, min_complete_fraction = 0.9,
                        min_genotype_group = 5L) {
  if (missing(conditional) || !is.character(conditional) || length(conditional) != 1L ||
      is.na(conditional) || !nzchar(conditional)) {
    stop("`conditional` must be given, e.g. \"P(J|D)\". The two anchorings are separate ",
         "scans over the same variants, so a result that does not say which it is will ",
         "be double-counted.", call. = FALSE)
  }
  subjects <- colnames(dosage)
  if (is.null(anchors)) anchors <- sort(unique(pairs$anchor_gene))

  out <- lapply(anchors, function(a) {
    pheno <- .qtl_pairing_phenotype(pairs, a, subjects, min_complete_fraction, min_subjects)
    if (is.null(pheno)) return(NULL)
    d <- dosage[, match(rownames(pheno), subjects), drop = FALSE]
    r <- qtlScanMultivariate(pheno, d, min_subjects)
    if (!nrow(r)) return(NULL)
    mgc <- qtlSmallestGenotypeClass(d)
    r[, c("min_genotype_group", "well_powered") :=
        list(mgc[get("variant")], mgc[get("variant")] >= min_genotype_group)]
    r[, "anchor_gene" := a]
    r
  })
  out <- data.table::rbindlist(Filter(Negate(is.null), out), use.names = TRUE)
  if (!nrow(out)) return(out)
  out[, "conditional" := conditional]
  data.table::setcolorder(out, c("conditional", "variant", "anchor_gene"))
  data.table::setorderv(out, "p_value")
  out[]
}

#' Per-cell follow-up under a significant pairing hit
#'
#' Regresses each single anchor-partner cell on the variant, so a significant omnibus can
#' be read as which partners moved and in which direction. Returns the fitted slope
#' alongside the raw phenotype means in the lowest and highest observed genotype class,
#' and how many subjects sat either side, which is what says whether a difference in
#' means is worth anything.
#'
#' A cell is only meaningful inside a row the omnibus already called significant, so the
#' scan is restricted to those variants and \code{omnibus_significant} travels with every
#' row. \code{marked} applies nominal alpha inside a cleared row; \code{marked_strict}
#' divides alpha by the row width, the family being that row's partner genes. Both are
#' returned so a downstream table and a figure make the same call rather than each
#' re-deriving it.
#'
#' @inheritParams pairingScan
#' @param omnibus A \code{\link{pairingScan}} result carrying a logical
#'   \code{significant} column. It may hold both anchorings stacked together, and
#'   normally should: a variant that cleared the omnibus on either side is worth reading
#'   cells for on both, and the per-row \code{omnibus_p_value} and
#'   \code{omnibus_significant} are matched on \code{variant} and \code{anchor_gene},
#'   so the other side's rows simply do not match this side's anchors.
#' @param alpha Nominal significance for a cell inside a cleared row. Default 0.05.
#'
#' @return A \code{data.table} with \code{conditional}, \code{variant},
#'   \code{anchor_gene}, \code{partner_gene}, \code{n}, \code{beta}, \code{p_value},
#'   \code{n_low}, \code{n_high}, \code{mean_low}, \code{mean_high}, \code{delta_mean},
#'   \code{omnibus_p_value}, \code{omnibus_significant}, \code{min_genotype_group},
#'   \code{marked} and \code{marked_strict}.
#' @seealso \code{\link{pairingScan}}
#' @export
pairingCellTests <- function(pairs, dosage, omnibus, conditional, alpha = 0.05,
                             min_subjects = 60L, min_complete_fraction = 0.9) {
  if (!"significant" %in% names(omnibus)) {
    stop("`omnibus` needs a `significant` column; a cell test is only meaningful inside ",
         "a row the omnibus already cleared.", call. = FALSE)
  }
  keep_variants <- unique(omnibus[omnibus$significant %in% TRUE][["variant"]])
  message(sprintf("pairingCellTests: %d variants with at least one significant %s row.",
                  length(keep_variants), conditional))
  if (!length(keep_variants)) return(data.table::data.table())

  subjects <- colnames(dosage)
  anchors <- sort(unique(pairs$anchor_gene))
  fits <- lapply(anchors, function(a) {
    pheno <- .qtl_pairing_phenotype(pairs, a, subjects, min_complete_fraction, min_subjects)
    if (is.null(pheno)) return(NULL)
    d <- dosage[rownames(dosage) %in% keep_variants,
                match(rownames(pheno), subjects), drop = FALSE]
    if (!nrow(d)) return(NULL)
    r <- qtlScanUnivariate(pheno, d, min_subjects)
    if (!nrow(r)) return(NULL)
    data.table::setnames(r, "phenotype", "partner_gene")
    r[, "anchor_gene" := a]
    merge(r, .qtl_class_means(pheno, d), by.x = c("variant", "partner_gene"),
          by.y = c("variant", "response_gene"), all.x = TRUE)
  })
  out <- data.table::rbindlist(Filter(Negate(is.null), fits), use.names = TRUE)
  if (!nrow(out)) return(out)

  om <- data.table::as.data.table(omnibus)[, list(
    variant = get("variant"), anchor_gene = get("anchor_gene"),
    omnibus_p_value = get("p_value"), omnibus_significant = get("significant"),
    min_genotype_group = get("min_genotype_group"))]
  out <- merge(out, om, by = c("variant", "anchor_gene"), all.x = TRUE)

  # Protected follow-up: nominal alpha inside a row the omnibus already cleared. The
  # family is that row's partner genes, so the strict column divides by the row width.
  out[, "row_width" := data.table::uniqueN(get("partner_gene"))]
  out[, c("conditional", "marked", "marked_strict") := list(
    conditional,
    get("omnibus_significant") %in% TRUE & get("p_value") < alpha,
    get("omnibus_significant") %in% TRUE & get("p_value") < alpha / get("row_width"))]

  out <- out[, list(conditional = get("conditional"), variant = get("variant"),
                    anchor_gene = get("anchor_gene"), partner_gene = get("partner_gene"),
                    n = get("n"), beta = get("beta"), p_value = get("p_value"),
                    n_low = get("n_low"), n_high = get("n_high"),
                    mean_low = get("mean_low"), mean_high = get("mean_high"),
                    delta_mean = get("delta_mean"),
                    omnibus_p_value = get("omnibus_p_value"),
                    omnibus_significant = get("omnibus_significant"),
                    min_genotype_group = get("min_genotype_group"),
                    marked = get("marked"), marked_strict = get("marked_strict"))]
  data.table::setorderv(out, c("variant", "anchor_gene", "p_value"))
  out[]
}
