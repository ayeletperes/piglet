#' Regroup Allele Similarity Clusters by a per-allele label
#'
#' Post-processes an ASC assignment (e.g. the \code{alleleClusterTable} produced by
#' \code{\link{inferAlleleClusters}}) by re-grouping alleles according to a per-allele
#' label -- typically the IMGT gene held in \code{imgt_allele}. It corrects the two
#' ways an ASC grouping and a gene-level view disagree: an ASC group that merges
#' several distinct genes, and a gene whose alleles are split across ASC groups. The
#' clustering itself is not re-run; only the output is relabelled, so this drops onto
#' the end of the pipeline without changing anything upstream.
#'
#' The label is not produced by the clustering; supply it on the table (the IMGT gene
#' in \code{imgt_allele}, and any other genes that share the allele in
#' \code{duplicated}). A gene is taken as the part of a label before \code{*}, so
#' allele suffixes and mutation tags (\code{IGHD1-20*01_a7g}) are ignored.
#'
#' Genes are grouped within a locus. A gene name is unique only inside its locus, so a
#' label stripped of its prefix (\code{"J1"} rather than \code{"IGHJ1"}) would otherwise
#' put IGHJ1, IGKJ1 and IGLJ1 in one bin. The locus is read from the label where it still
#' carries one and from \code{group_col} otherwise; if neither supplies it, the grouping
#' falls back to plain gene names and warns. \code{action = "split"} is unaffected, since
#' its output is prefixed by the ASC group.
#'
#' @param x A \code{data.frame} / \code{AlleleClusterTable} carrying at least the ASC
#'   group column and a label column.
#' @param group_col Name of the ASC group column. Default \code{"allele_cluster"}.
#'   An ASC-allele column such as \code{"asc_allele"} works too: an allele suffix
#'   (\code{*NN}) is stripped, so \code{IGHD1-CO5H*01} collapses to the ASC gene
#'   \code{IGHD1-CO5H}.
#' @param label_col Name of the label column (e.g. \code{"imgt_allele"} or
#'   \code{"iuis_allele"}). A single value names one gene; a value listing several
#'   genes separated by \code{,} or \code{/} (a shared / ambiguous allele) names them
#'   all and will force those genes together.
#' @param duplicated_col Optional column naming \emph{other} genes that share the
#'   same allele, separated by \code{,} or \code{/}. Default \code{"duplicated"};
#'   ignored if the column is absent or is a logical flag (which carries no gene
#'   names). Not needed when the shared genes are already listed in \code{label_col}.
#' @param action One of:
#'   \describe{
#'     \item{\code{"both"}}{(default) re-group so each label becomes one bin; a shared
#'       allele (a label present under two genes) forces those genes together.}
#'     \item{\code{"split"}}{only break ASC groups that mix labels; a label split
#'       across two ASC groups stays split.}
#'     \item{\code{"merge"}}{only merge ASC groups that share a label or a shared
#'       allele; groups that mix labels are left intact.}
#'   }
#'
#' @return \code{x} with an added character column \code{regrouped}. A group made of
#'   several genes (a merge, e.g. a shared allele) is tagged with all of them joined
#'   by \code{/}, e.g. \code{"IGHD2-2/IGHD2-8"}, so no gene name is dropped.
#'
#' @examples
#' tb <- data.frame(
#'   allele         = paste0("a", 1:5),
#'   allele_cluster = c("G1", "G1", "G2", "G3", "G3"),
#'   imgt_allele    = c("IGHD1-1*01", "IGHD1-2*01", "IGHD1-1*02",
#'                      "IGHD2-2*01", "IGHD2-8*01"),
#'   duplicated     = c("", "", "", "IGHD2-8*01", ""),   # a4 shared: 2-2 and 2-8
#'   stringsAsFactors = FALSE)
#' regroupASCByLabel(tb, action = "both")$regrouped
#'
#' @seealso \code{\link{inferAlleleClusters}}
#' @export
regroupASCByLabel <- function(x,
                              group_col = "allele_cluster",
                              label_col = "imgt_allele",
                              duplicated_col = "duplicated",
                              action = c("both", "split", "merge")) {
  action <- match.arg(action)
  if (!group_col %in% names(x)) stop("group_col '", group_col, "' not found in x")
  if (!label_col %in% names(x)) stop("label_col '", label_col, "' not found in x")

  ## the ASC group is the group column with any allele suffix (*NN) removed, so an
  ## ASC-allele column (IGHD1-CO5H*01) collapses to its ASC gene (IGHD1-CO5H); a
  ## plain group id (a number, or IGHDF1-G1) is left unchanged.
  grp <- sub("\\*.*", "", as.character(x[[group_col]]))
  ## other genes sharing the allele; a logical flag (TRUE/FALSE) carries no gene
  ## names, so it is ignored -- supply gene names (comma/slash separated) to use it.
  dup <- if (!is.null(duplicated_col) && duplicated_col %in% names(x) &&
             !is.logical(x[[duplicated_col]]))
    as.character(x[[duplicated_col]]) else rep("", nrow(x))

  ## per-allele set of genes: from label_col plus duplicated_col; gene = text before '*'
  gene_of <- function(v) {
    v <- trimws(v[!is.na(v) & nzchar(trimws(v))])
    unique(sub("\\*.*", "", v))
  }
  labels_of <- Map(function(a, b)
                     gene_of(c(unlist(strsplit(a, "[,/]")), unlist(strsplit(b, "[,/]")))),
                   as.character(x[[label_col]]), dup)
  prim <- vapply(labels_of, function(z) if (length(z)) z[1] else NA_character_,
                 character(1))

  if (action == "split") {
    x$regrouped <- paste0(grp, " | ", prim)
    return(x)
  }

  ## A gene name is unique only within a locus, and a label that has lost its locus
  ## prefix ("J1" rather than "IGHJ1") is one string for all three loci, so the
  ## union-find below would fuse IGHJ1, IGKJ1 and IGLJ1 into a single bin. Each row is
  ## given a locus, from its own label where that still carries one and otherwise from
  ## the ASC group, and the union-find runs on locus-qualified keys.
  locus_of <- function(v) {
    p <- substr(as.character(v), 1L, 3L)
    ifelse(!is.na(p) & grepl("^IG[HKL]$", p), p, NA_character_)
  }
  row_locus <- locus_of(prim)
  from_grp <- locus_of(grp)
  row_locus[is.na(row_locus)] <- from_grp[is.na(row_locus)]
  if (anyNA(row_locus)) {
    warning("no locus could be read from ", label_col, " or ", group_col, " for ",
            sum(is.na(row_locus)), " row(s); genes of the same name in different loci ",
            "will be grouped together")
    row_locus[is.na(row_locus)] <- ""
  }

  ## union-find over the distinct locus-qualified labels
  keys_of <- Map(function(ll, lo) if (length(ll)) paste0(lo, "\r", ll) else character(),
                 labels_of, row_locus)
  prim_key <- ifelse(is.na(prim), NA_character_, paste0(row_locus, "\r", prim))
  labs <- unique(unlist(keys_of))
  labs <- labs[!is.na(labs)]
  parent <- stats::setNames(labs, labs)
  find <- function(z) { while (parent[[z]] != z) z <- parent[[z]]; z }
  link <- function(a, b) { ra <- find(a); rb <- find(b); if (ra != rb) parent[[rb]] <<- ra }

  ## shared / ambiguous alleles link the genes they name (used by both & merge)
  for (ll in keys_of) if (length(ll) > 1)
    for (j in 2:length(ll)) link(ll[[1L]], ll[[j]])

  ## merge also links labels that co-occur in the same ASC group
  if (action == "merge") for (g in unique(grp)) {
    gl <- unique(stats::na.omit(prim_key[grp == g]))
    if (length(gl) > 1) for (j in 2:length(gl)) link(gl[[1L]], gl[[j]])
  }

  ## tag each merged group with ALL its genes, so a shared allele keeps both names
  ## (e.g. a merge of IGHD2-2 and IGHD2-8 is labelled "IGHD2-2/IGHD2-8", not one of them)
  gene_part <- function(z) sub("^[^\r]*\r", "", z)
  root <- vapply(labs, find, character(1))
  tag  <- vapply(split(gene_part(labs), root),
                 function(m) paste(sort(unique(m)), collapse = "/"), character(1))

  ## Two groups kept apart by locus still carry the same text when the label had no
  ## prefix to begin with ("J1" for both IGHJ1 and IGKJ1), which would fuse them again
  ## the moment anything groups on the result. Those are prefixed with their locus.
  root_locus <- vapply(split(sub("\r.*$", "", labs), root),
                       function(m) { u <- unique(m[nzchar(m)])
                                     if (length(u) == 1L) u else NA_character_ },
                       character(1))
  ambiguous <- tag %in% tag[duplicated(tag)] & !is.na(root_locus)
  tag[ambiguous] <- paste0(root_locus[ambiguous], ":", tag[ambiguous])

  spans <- vapply(split(gene_part(labs), root),
                  function(m) length(unique(stats::na.omit(locus_of(m)))) > 1L, logical(1))
  if (any(spans)) {
    warning(sum(spans), " group(s) name genes from more than one locus: ",
            paste(utils::head(tag[spans], 3), collapse = ", "))
  }
  x$regrouped <- vapply(prim_key, function(p) if (is.na(p)) NA_character_ else tag[[find(p)]],
                        character(1))
  x
}
