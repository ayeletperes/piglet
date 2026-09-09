# Adjudication of pairing leads: is a significant omnibus a real reallocation of partners,
# or a rigid shift, one subject, or marginal usage in disguise?
#
# Moved from husa_manuscript/lib_manuscript/analysis/igqtl.R (characterise_leads,
# overlap_coefficient), generalised from the D/J hardcoding to any anchor and partner.

# Overlap of two normal densities separated by d pooled standard deviations. Reported
# beside Cohen's d because "0.8" means little to most readers and "these two genotype
# groups share 69% of their range" means quite a lot.
.qtl_overlap <- function(d) 2 * stats::pnorm(-abs(d) / 2)

#' Characterise pairing leads: what kind of signal is this?
#'
#' A significant omnibus says a variant is associated with an anchor gene's partner
#' distribution. It does not say the partners were reallocated: the same p-value is
#' produced by a rigid shift of the whole profile, by one unusual subject, and by a
#' variant that simply changes how much the anchor or a partner is used. This runs the
#' follow-ups that separate those, per lead, and returns a verdict.
#'
#' @section What is tested:
#' \describe{
#'   \item{Rigid shift against reallocation}{The mean across partners regressed on dosage
#'     is the shift every partner shares (\code{p_common}). Removing that mean and testing
#'     what is left asks whether partners moved by \emph{different} amounts
#'     (\code{p_heterogeneous}). The second is the reallocation claim; the first is not.
#'     One deviation column is dropped before the test because deviations from a row mean
#'     sum to zero and the covariance would be singular.}
#'   \item{Separation}{\code{mahalanobis}, the distance between the highest and lowest
#'     genotype class centroids in units of the pooled within-genotype covariance, after
#'     residualising on the per-partner slopes. With \code{overlap_joint} beside it.}
#'   \item{One subject}{\code{p_worst_loo} refits the omnibus dropping each subject in turn
#'     and keeps the worst p-value. A signal that rests on one person shows up here and
#'     nowhere else, and a small genotype class is exactly when it happens.}
#'   \item{Both margins}{The enrichment divides both marginals out, so a lead is not
#'     pairing-specific until neither margin moved. Logit marginal usage is regressed on
#'     dosage for the anchor (\code{anchor_usage_p}) and for every partner
#'     (\code{partner_usage_p_min}, Bonferroni across partners).}
#'   \item{Ancestry}{\code{ancestry_p}, whether the dosage itself is stratified, and
#'     \code{p_replication}, whether the omnibus holds inside one ancestry group. A cis
#'     variant is expected to track ancestry, so neither is a rejection on its own.}
#'   \item{Position}{\code{rho_position}, Spearman of per-partner slope against partner
#'     position. A monotone trend along the locus is a different claim from a scattered
#'     one. Needs \code{positions}.}
#' }
#'
#' @section The verdict:
#' \code{verdict} is assigned in order: \code{single_subject_driven},
#' \code{overlapping_small_effect}, \code{pairing_plus_both_margins},
#' \code{pairing_plus_anchor_usage}, \code{pairing_plus_partner_usage}, and otherwise
#' \code{pairing_only}. Its cutoffs are arguments rather than constants, and their defaults
#' are the ones the heavy-chain D/J analysis used. They are that analysis's policy, not
#' facts about the method: a different cohort or a different pair of segments should set
#' them deliberately rather than inherit them.
#'
#' @param leads Significant leads for one anchoring, from \code{\link{runPairingQTL}},
#'   carrying \code{variant}, \code{anchor_gene}, \code{p_value} and \code{threshold}.
#' @param pairs The \code{\link{pairingTable}} the scan used.
#' @param dosage Numeric matrix of genotype dosages, variants in rows, subjects in columns.
#' @param positions Optional \code{data.frame} of \code{partner_gene} and
#'   \code{partner_position}. Without it \code{rho_position} is \code{NA} and the
#'   per-partner table is left in column order.
#' @param ancestry Optional \code{data.frame} of \code{subject} and \code{ancestry}.
#'   Without it \code{ancestry_p} and \code{p_replication} are \code{NA}.
#' @param replication_ancestry Ancestry group to re-test within. Default \code{"EUR"}.
#' @param min_replication_n Minimum subjects in that group before it is attempted.
#'   Default 40.
#' @param min_genotype_group Smallest genotype class at or above which a lead is called
#'   well powered. Default 5.
#' @param alpha Nominal significance for the per-partner and margin tests. Default 0.05.
#' @param min_sep_d,max_overlap A partner counts as separated at or above this Cohen's d
#'   and at or below this overlap. Defaults 0.8 and 0.70.
#' @param min_mahalanobis Joint separation a lead needs before it is called separated at
#'   all. Default 1.
#' @inheritParams pairingScan
#'
#' @return A list of two \code{data.table}s: \code{leads}, the input with the summary
#'   columns and \code{verdict} merged on, and \code{per_partner}, one row per lead and
#'   partner with slopes and effect sizes.
#' @seealso \code{\link{runPairingQTL}}, \code{\link{pairingCellTests}}
#' @export
pairingLeadCharacter <- function(leads, pairs, dosage, positions = NULL, ancestry = NULL,
                                 replication_ancestry = "EUR", min_replication_n = 40L,
                                 min_genotype_group = 5L, alpha = 0.05,
                                 min_sep_d = 0.8, max_overlap = 0.70, min_mahalanobis = 1,
                                 min_subjects = 60L, min_complete_fraction = 0.9) {
  leads <- data.table::as.data.table(leads)
  needed <- c("variant", "anchor_gene", "threshold")
  if (!all(needed %in% names(leads))) {
    stop("`leads` needs columns: ", paste(needed, collapse = ", "), ".", call. = FALSE)
  }
  if (!nrow(leads)) return(list(leads = leads, per_partner = data.table::data.table()))

  subjects <- colnames(dosage)
  anc <- if (is.null(ancestry)) {
    data.table::data.table(subject = character(0), ancestry = character(0))
  } else {
    data.table::as.data.table(ancestry)[, list(subject = as.character(get("subject")),
                                               ancestry = get("ancestry"))]
  }

  # Both marginals, straight off the pairing table: p_anchor and p_partner already are
  # anchor_total/depth and partner_total/depth.
  anchor_usage <- unique(pairs[, list(subject = get("subject"), anchor_gene = get("anchor_gene"),
                                      own = get("p_anchor"))])
  partner_wide <- data.table::dcast(
    unique(pairs[, list(subject = get("subject"), partner_gene = get("partner_gene"),
                        own = get("p_partner"))]),
    subject ~ partner_gene, value.var = "own")

  one <- function(i) {
    ag <- leads$anchor_gene[i]
    v <- leads$variant[i]
    pheno <- pairingPhenotype(pairs, ag, subjects, min_complete_fraction, min_subjects)
    if (is.null(pheno) || !v %in% rownames(dosage)) return(NULL)

    g <- dosage[v, match(rownames(pheno), subjects)]
    ok <- is.finite(g)
    y <- pheno[ok, , drop = FALSE]
    g <- g[ok]
    ids <- rownames(y)
    cls <- round(g)
    tab <- table(cls)
    lo <- which(cls == as.numeric(names(tab)[1L]))
    hi <- which(cls == as.numeric(names(tab)[length(tab)]))
    gc <- g - mean(g)

    per_partner <- data.table::rbindlist(lapply(colnames(y), function(pg) {
      val <- y[, pg]
      beta <- sum(gc * (val - mean(val))) / sum(gc^2)
      resid <- val - mean(val) - beta * gc
      se <- sqrt(sum(resid^2) / (length(val) - 2) / sum(gc^2))
      a <- val[lo]
      b <- val[hi]
      pooled <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
                       (length(a) + length(b) - 2))
      sep <- if (is.finite(pooled) && pooled > 0) (mean(b) - mean(a)) / pooled else NA_real_
      wins <- sum(outer(b, a, ">")) + 0.5 * sum(outer(b, a, "=="))
      data.table::data.table(
        anchor_gene = ag, variant = v, partner_gene = pg,
        partner_position = .qtl_partner_position(pg, positions),
        beta_per_allele = beta, se = se, t_stat = beta / se,
        p_value = 2 * stats::pt(-abs(beta / se), length(val) - 2),
        mean_low = mean(a), mean_high = mean(b), sep_d = sep,
        auc = wins / (length(a) * length(b)), overlap = .qtl_overlap(sep))
    }))
    if (!all(is.na(per_partner$partner_position))) {
      data.table::setorderv(per_partner, "partner_position")
    }

    # A shift shared by every partner, and what is left after removing it. The second is
    # the test of whether the partners move by different amounts, which is what separates
    # a reallocation from a rigid shift of the whole profile.
    beta_by_col <- per_partner$beta_per_allele[match(colnames(y), per_partner$partner_gene)]
    row_mean <- rowMeans(y)
    p_common <- summary(stats::lm(row_mean ~ gc))$coefficients["gc", "Pr(>|t|)"]
    dev <- y - row_mean
    p_hetero <- .qtl_pillai(dev[, -ncol(dev), drop = FALSE], g)[["p"]]

    # Joint separation, in units of the pooled within-genotype covariance.
    resid <- y - outer(gc, beta_by_col)
    delta <- colMeans(y[hi, , drop = FALSE]) - colMeans(y[lo, , drop = FALSE])
    mahal <- tryCatch(sqrt(max(drop(delta %*% solve(stats::cov(resid)) %*% delta), 0)),
                      error = function(e) NA_real_)

    # A tiny genotype class only matters if the signal rests on one person, so test that.
    p_loo <- max(vapply(seq_len(nrow(y)),
                        function(k) .qtl_pillai(y[-k, , drop = FALSE], g[-k])[["p"]], 1),
                 na.rm = TRUE)

    # Both margins. The enrichment removes each of them, so both must be tested before a
    # lead can be called usage neutral.
    logit_fit <- function(u) {
      x <- stats::qlogis(pmin(pmax(u, 1e-6), 1 - 1e-6))
      if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) return(c(NA_real_, NA_real_))
      cf <- summary(stats::lm(x ~ g))$coefficients
      c(cf[2L, "Estimate"], cf[2L, "Pr(>|t|)"])
    }
    a_own <- anchor_usage[anchor_usage$anchor_gene == ag][match(ids, get("subject")), get("own")]
    a_fit <- logit_fit(a_own)
    p_own <- as.matrix(partner_wide[match(ids, partner_wide$subject), -1L])
    p_p <- apply(p_own, 2L, function(u) logit_fit(u)[2L])
    p_moved <- names(sort(p_p[is.finite(p_p) & p_p < alpha / length(p_p)]))

    # Within-ancestry replication, since a cis variant is expected to track ancestry.
    rep_ids <- intersect(ids, anc$subject[anc$ancestry == replication_ancestry])
    rep_p <- if (length(rep_ids) >= min_replication_n) {
      gr <- dosage[v, match(rep_ids, subjects)]
      keep <- is.finite(gr)
      .qtl_pillai(pheno[rep_ids[keep], , drop = FALSE], gr[keep])[["p"]]
    } else {
      NA_real_
    }
    anc_p <- tryCatch({
      a <- merge(data.table::data.table(subject = ids, dose = g), anc, by = "subject")
      stats::anova(stats::lm(dose ~ ancestry, a))$`Pr(>F)`[1L]
    }, error = function(e) NA_real_)

    best <- per_partner[which.max(abs(per_partner$sep_d))]
    rho <- suppressWarnings(stats::cor(per_partner$partner_position,
                                       per_partner$beta_per_allele, method = "spearman"))
    list(summary = data.table::data.table(
           anchor_gene = ag, variant = v, min_genotype_group = min(tab),
           n_class = length(tab), p_worst_loo = p_loo, mahalanobis = mahal,
           overlap_joint = .qtl_overlap(mahal),
           best_partner_gene = best$partner_gene, best_sep_d = best$sep_d, best_auc = best$auc,
           n_partner_up = sum(per_partner$beta_per_allele > 0),
           n_partner_down = sum(per_partner$beta_per_allele < 0),
           n_partner_nominal = sum(per_partner$p_value < alpha, na.rm = TRUE),
           n_partner_separated = sum(abs(per_partner$sep_d) >= min_sep_d &
                                       per_partner$overlap <= max_overlap, na.rm = TRUE),
           p_common = p_common, p_heterogeneous = p_hetero, rho_position = rho,
           anchor_usage_beta = a_fit[1L], anchor_usage_p = a_fit[2L],
           partner_usage_p_min = suppressWarnings(min(p_p, na.rm = TRUE)),
           partner_usage_n_moved = length(p_moved),
           partner_usage_genes = paste(p_moved, collapse = "; "),
           ancestry_p = anc_p, n_replication = length(rep_ids), p_replication = rep_p),
         per_partner = per_partner)
  }

  parts <- Filter(Negate(is.null), lapply(seq_len(nrow(leads)), one))
  if (!length(parts)) return(list(leads = leads[0], per_partner = data.table::data.table()))
  detail <- data.table::rbindlist(lapply(parts, `[[`, "summary"), use.names = TRUE)

  # min_genotype_group and well_powered are recomputed here on the lead's own subject set,
  # so the scan's copies are dropped to keep the merge from suffixing duplicates.
  dup <- intersect(c(setdiff(names(detail), c("anchor_gene", "variant")), "well_powered"),
                   names(leads))
  out <- merge(leads[, !dup, with = FALSE], detail, by = c("anchor_gene", "variant"),
               sort = FALSE)
  # Held in a differently-named local first: `min_genotype_group` is also a column of
  # `out`, and inside `[` the column shadows the argument, so the comparison would be the
  # column against itself and every lead would come out well powered.
  well_powered_at <- min_genotype_group
  out[, "well_powered" := get("min_genotype_group") >= well_powered_at]
  out[, "separated" := is.finite(get("mahalanobis")) &
        get("mahalanobis") >= min_mahalanobis & get("n_partner_separated") >= 1L]
  out[, "verdict" := data.table::fcase(
    get("p_worst_loo") > get("threshold"), "single_subject_driven",
    !get("separated"), "overlapping_small_effect",
    get("anchor_usage_p") < alpha & get("partner_usage_n_moved") > 0L, "pairing_plus_both_margins",
    get("anchor_usage_p") < alpha, "pairing_plus_anchor_usage",
    get("partner_usage_n_moved") > 0L, "pairing_plus_partner_usage",
    default = "pairing_only")]
  data.table::setorderv(out, "p_value")

  tally <- out[, .N, by = "verdict"]
  message("pairingLeadCharacter: ",
          paste(sprintf("%s %d", tally$verdict, tally$N), collapse = ", "))

  list(leads = out[],
       per_partner = data.table::rbindlist(lapply(parts, `[[`, "per_partner"), use.names = TRUE))
}

.qtl_partner_position <- function(gene, positions) {
  if (is.null(positions)) return(NA_real_)
  p <- data.table::as.data.table(positions)
  if (!all(c("partner_gene", "partner_position") %in% names(p))) {
    stop("`positions` needs columns partner_gene and partner_position.", call. = FALSE)
  }
  p$partner_position[match(gene, p$partner_gene)]
}
