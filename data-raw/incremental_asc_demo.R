# =============================================================================
# Demo: experimental incremental / hysteresis ASC framework
# -----------------------------------------------------------------------------
# Runs two experiments with runIncrementalASCExperiment() and writes the
# diagnostic outputs the framework produces:
#   * ASC retention before vs. after adding sequences
#   * number of old alleles that changed ASC
#   * number of new ASCs created
#   * candidate splits and merges
#   * Leiden stability of each proposed change
#   * resolution stability
#   * within- vs. between-group similarity for candidate splits
#   * silhouette change
#   * subsampling stability
#
# Experiment 1 (always runnable): HVGERM hold-out. A batch of human IGHV alleles
#   -- including whole subgroups, to mimic "previously unrepresented structure"
#   -- is held out, the remainder is clustered as the baseline, then the held-out
#   alleles are added back and the framework measures what changes.
#
# Experiment 2 (cross-species, if the FASTA is present): add rhesus (macaque)
#   IGHV onto the human reference -- the genuine "new species, unknown orthology"
#   case.
#
# Run with:  Rscript data-raw/incremental_asc_demo.R
# =============================================================================

if (Sys.getenv("ASC_DEV") == "1") {
  suppressMessages(devtools::load_all(Sys.getenv("ASC_PKG", ".")))
} else {
  suppressMessages(library(piglet))
}

outdir <- Sys.getenv("ASC_DEMO_OUTDIR", "data-raw/incremental_asc_demo_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- paths to the "real" data -----------------------------------------------
## Human "husa larger" reference: the piglet all-loci ASC table. It carries both
## the germline sequences (seq_gapped) and the established ASC assignment
## (asc_allele); the ASC *cluster* is asc_allele with the *NN allele number
## stripped (614 IGHV alleles in 56 clusters).
husa_tsv <- "/home/ayelet/dropbox/husa_manuscript/results/asc_final/current/VDJ_piglet_all_loci.tsv"
## Rhesus: OGRDB/IMGT macaque IGHV. This IMGT file is usable now:
rhesus_fasta <- "/home/ayelet/Documents/Macaca-mulatta_Ig_Heavy_V_1-1.fasta"

## Load the human IGHV reference + established ASC clusters from the husa table.
load_husa_ighv <- function(path) {
  t <- utils::read.delim(path, stringsAsFactors = FALSE)
  v <- t[t$tag == "V" & grepl("^IGHV", t$allele), ]
  seqs <- stats::setNames(v$seq_gapped, v$allele)
  cluster <- sub("[*].*", "", v$asc_allele)              # ASC cluster (multi-allele)
  baseline <- stats::setNames(cluster, v$allele)
  list(seqs = seqs, baseline = baseline)
}

write_diagnostics <- function(res, tag) {
  f <- function(name) file.path(outdir, paste0(tag, "_", name, ".tsv"))
  utils::write.table(res$diagnostics$summary_counts, f("summary_counts"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(res$retention, f("retention"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(res$changes, f("changes"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(res$reconciled, f("reconciled"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  message("[", tag, "] wrote diagnostics to ", outdir)
  print(res)
  ## quick retention plot (before vs after): retention == 1 means fully preserved
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(res$retention)) {
    p <- ggplot2::ggplot(res$retention,
                         ggplot2::aes(x = stats::reorder(old_asc, retention),
                                      y = retention)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "original ASC", y = "retention (max overlap)",
                    title = paste0(tag, ": ASC retention after adding alleles")) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(outdir, paste0(tag, "_retention.png")), p,
                    width = 6, height = max(3, nrow(res$retention) * 0.15))
  }
}

# =============================================================================
# Experiment 1: HVGERM hold-out (self-contained)
# =============================================================================
data(HVGERM)

## Hold out two whole subgroups (unrepresented structure) plus a random batch.
subgroup <- sub("[-*].*", "", names(HVGERM))          # e.g. IGHV3
held_subgroups <- c("IGHV6", "IGHV7")                  # small subgroups
set.seed(42)
holdout <- unique(c(
  names(HVGERM)[subgroup %in% held_subgroups],
  sample(names(HVGERM), 80)
))
reference <- HVGERM[setdiff(names(HVGERM), holdout)]
new_seqs  <- HVGERM[holdout]

message(sprintf("Experiment 1: reference=%d, new=%d alleles",
                length(reference), length(new_seqs)))
res1 <- runIncrementalASCExperiment(
  reference, new_seqs,
  sim_method = "log",          # gentler, more resolution-stable than linear
  resolution = NULL,           # baseline is self-clustered, so a sweep is fair here
  resolution_grid = seq(0.15, 0.55, by = 0.02),
  n_leiden = 100, n_subsample = 20, quiet = FALSE)
write_diagnostics(res1, "holdout")

# =============================================================================
# Experiment 2: rhesus (macaque) added onto the husa human reference
# =============================================================================
if (file.exists(husa_tsv) && file.exists(rhesus_fasta)) {
  husa <- load_husa_ighv(husa_tsv)
  message(sprintf("husa human IGHV reference: %d alleles in %d ASC clusters",
                  length(husa$seqs), length(unique(husa$baseline))))

  rhesus <- tigger::readIgFasta(rhesus_fasta)          # named character vector
  set.seed(7)
  rhesus_batch <- rhesus[sample(seq_along(rhesus), min(150, length(rhesus)))]
  # prefix so macaque names can never collide with human allele names (labels
  # only - clustering is on sequence, never on the species tag)
  names(rhesus_batch) <- paste0("rh_", names(rhesus_batch))

  message(sprintf("Experiment 2: human reference=%d, rhesus new=%d alleles",
                  length(husa$seqs), length(rhesus_batch)))
  res2 <- runIncrementalASCExperiment(
    husa$seqs, rhesus_batch,
    baseline_assign = husa$baseline,   # start from the established husa ASCs
    sim_method = "log",
    ## resolution calibrated to reproduce ~56 clusters (the husa ASC count) so
    ## retention reflects the *added data*, not a resolution mismatch against the
    ## established structure. Re-calibrate if the reference changes.
    resolution = 0.46,
    n_leiden = 100, n_subsample = 20, quiet = FALSE)
  write_diagnostics(res2, "rhesus_on_husa")
} else {
  message("husa TSV or rhesus FASTA not found -- skipping the cross-species ",
          "experiment. Confirm `husa_tsv` and `rhesus_fasta`.")
}

message("Done. Outputs in: ", normalizePath(outdir))
