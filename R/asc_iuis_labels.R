# ASC subgroup to IUIS gene-group labelling.
#
# Moved from husa_manuscript/lib_manuscript/analysis/repertoire_prep.R, which built this
# vocabulary inline. Owning it here means the labels a downstream analysis joins on come
# from one checked function rather than from a copy of a script.

# "D5-12/5-18" expands to D5-12, D5-18: parts after the first drop the segment letter.
.asc_expand_iuis_group <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  parts <- trimws(unlist(strsplit(x, "/", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  full <- parts[grepl("^[A-Za-z]", parts)]
  if (length(full)) {
    tag <- sub("^([A-Za-z]+).*", "\\1", full[1])
    parts <- ifelse(grepl("^[A-Za-z]", parts), parts, paste0(tag, parts))
  }
  unique(parts)
}

# The inverse: D5-12, D5-18, D5-5 compacts back to D5-12/5-18/5-5.
.asc_format_iuis_group <- function(parts) {
  parts <- sort(unique(parts))
  if (!length(parts)) return(NA_character_)
  tags <- unique(sub("^([A-Za-z]+).*", "\\1", parts))
  if (length(tags) > 1L) return(paste(parts, collapse = "/"))
  if (length(parts) == 1L) return(parts)
  paste0(parts[1], "/", paste0(sub(paste0("^", tags[1]), "", parts[-1]), collapse = "/"))
}

# A call can be a comma-separated list; map each part and keep the result only when the
# parts agree. A call whose alleles name different IUIS genes is not a single-gene
# assignment and is left unlabelled.
.asc_map_one <- function(voc, calls) {
  if (is.na(calls) || !nzchar(calls)) return(NA_character_)
  m <- unique(voc[unlist(strsplit(calls, ",", fixed = TRUE))])
  if (length(m) == 1L && !is.na(m)) m else NA_character_
}

#' Build the ASC subgroup to IUIS gene-group vocabulary
#'
#' Derives the display labels that translate an Allele Similarity Cluster assignment into
#' IUIS gene names, and returns the allele set behind each label.
#'
#' Two ASC subgroups that share any IUIS gene must carry the same label, so subgroups and
#' IUIS parts are treated as a bipartite graph and every connected component gets one
#' merged label. Where a component holds more than one subgroup the label alone does not
#' identify the subgroup, so \code{asc_and_iuis_group} suffixes them A, B, C for display.
#'
#' D is handled differently from V and J, and that asymmetry is the reason ASC labels look
#' inconsistent downstream. D is relabelled at IUIS gene resolution by
#' \code{\link{regroupASCByLabel}}: an ASC group holding alleles of more than one IUIS gene
#' is broken apart and groups naming the same gene are joined, so the D label is a property
#' of the allele rather than of the ASC gene. V and J keep the ASC grouping. This is why
#' the resulting labels carry the locus for D (\code{IGHD5-12}) but not for V and J
#' (\code{V1-18}, \code{J4}).
#'
#' \strong{The label is a label.} It is a display string, and it moves when the reference
#' set behind it moves. What identifies a group is the set of member allele sequences,
#' returned in \code{identity}; join on that, or compare it to detect that an upstream
#' reference changed underneath a downstream table.
#'
#' @param husa A \code{data.frame} of the ASC reference (the \code{husa.tsv} produced by
#'   the ASC pipeline), with columns \code{allele} (the ASC allele name), \code{chain},
#'   \code{gene_type}, \code{seq} and \code{iuis_allele}. \code{iuis_allele} may list
#'   several IUIS alleles separated by \code{,} or \code{/}.
#'
#' @return A list of class \code{ascIUISVocabulary}:
#'   \describe{
#'     \item{\code{labels}}{one row per ASC subgroup: \code{chain}, \code{subgroup},
#'       \code{iuis_group}, \code{asc_and_iuis_group}, \code{iuis_group_size}.}
#'     \item{\code{asc_to_iuis}}{named character vector, ASC subgroup to IUIS label, for V
#'       and J.}
#'     \item{\code{d_allele_to_iuis}}{named character vector, ASC \emph{allele} to IUIS
#'       label, for D.}
#'     \item{\code{identity}}{one row per \code{chain} and \code{iuis_label}:
#'       \code{asc_genes}, \code{n_alleles}, \code{member_alleles} and
#'       \code{member_seq_md5}. This is what the label means; the label itself is not.}
#'   }
#'
#' @seealso \code{\link{annotateRepertoireIUIS}}, \code{\link{regroupASCByLabel}}
#' @examples
#' husa <- data.frame(
#'   allele      = c("IGHV1-2*01", "IGHV1-2*02", "IGHD1-CO5H*01"),
#'   chain       = "IGH",
#'   gene_type   = c("IGHV", "IGHV", "IGHD"),
#'   seq         = c("ACGTACGT", "ACGTACGA", "GGTTGGTT"),
#'   iuis_allele = c("IGHV1-2*01", "IGHV1-2*02", "IGHD1-20*01"))
#'
#' voc <- ascIUISVocabulary(husa)
#' voc$labels
#'
#' # What each label stands for. The label is a display string that moves when the
#' # reference set moves; these member alleles are what identify the group.
#' voc$identity
#' @export
ascIUISVocabulary <- function(husa) {
  husa <- data.table::as.data.table(husa)

  # The ASC reference writes the IUIS assignment as `iuis_allele` in newer runs and as
  # `husa` in older ones. Accept either and say which was taken, rather than defaulting a
  # missing column to NA and labelling nothing.
  if (!"iuis_allele" %in% names(husa)) {
    if (!"husa" %in% names(husa)) {
      stop("`husa` carries neither an `iuis_allele` nor a `husa` column, so no IUIS ",
           "assignment is available.", call. = FALSE)
    }
    message("ascIUISVocabulary: no `iuis_allele` column; using `husa` as the IUIS assignment.")
    # Pulled out before the assignment on purpose: there is a column called `husa`, and
    # inside `[` that name shadows the table itself, so `husa[["husa"]]` would index the
    # column by name and fail.
    iuis <- husa[["husa"]]
    husa[, "iuis_allele" := iuis]
  }

  required <- c("allele", "chain", "gene_type", "seq", "iuis_allele")
  missing <- setdiff(required, names(husa))
  if (length(missing)) {
    stop("`husa` is missing required column(s): ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }

  # ---- D: relabelled per allele at IUIS gene resolution ----------------------------
  asc_iuis_alleles <- unique(husa[, list(asc_allele = get("allele"),
                                         iuis_allele = get("iuis_allele"),
                                         seq = get("seq"))])
  asc_iuis_alleles[, "iuis_allele_clean" := paste0(
    unique(unlist(strsplit(unlist(strsplit(get("iuis_allele"), "/")), ","))), collapse = ","),
    by = "iuis_allele"]

  d_rows <- substr(asc_iuis_alleles$asc_allele, 4L, 4L) == "D"
  d_allele_to_iuis <- character(0)
  if (any(d_rows)) {
    d_iuis <- regroupASCByLabel(as.data.frame(asc_iuis_alleles[d_rows]),
                                group_col = "asc_allele", label_col = "iuis_allele_clean",
                                action = "both")
    d_allele_to_iuis <- stats::setNames(d_iuis$regrouped, d_iuis$asc_allele)
  }

  # ---- V and J: the ASC subgroup keeps its grouping --------------------------------
  husa_set <- unique(husa[, list(allele = get("allele"), chain = get("chain"),
                                 gene_type = get("gene_type"),
                                 iuis_allele = get("iuis_allele"))])
  husa_set[, "subgroup" := alakazam::getGene(get("allele"), strip_d = FALSE,
                                             omit_nl = FALSE, first = FALSE)]
  husa_set[, "iuis_group" := paste0(gsub("IG[KLH]", "", unique(alakazam::getGene(
    unlist(strsplit(get("iuis_allele"), ",", fixed = TRUE)),
    strip_d = FALSE, collapse = TRUE, first = FALSE))), collapse = "/"), by = "subgroup"]

  labels <- unique(husa_set[, list(subgroup = get("subgroup"), iuis_group = get("iuis_group"))])
  labels[, "chain" := substr(get("subgroup"), 1L, 3L)]

  label_parts <- labels[!is.na(get("iuis_group")) & nzchar(get("iuis_group")),
                        list(iuis_part = .asc_expand_iuis_group(get("iuis_group"))),
                        by = c("chain", "subgroup", "iuis_group")]

  g <- igraph::graph_from_data_frame(
    label_parts[, list(from = paste(get("chain"), get("subgroup"), sep = "::subgroup::"),
                       to = paste(get("chain"), get("iuis_part"), sep = "::part::"))],
    directed = FALSE)
  memb <- igraph::components(g)$membership
  comp <- data.table::data.table(node = names(memb), component = as.integer(memb))

  subgroup_comp <- comp[grepl("::subgroup::", get("node")),
                        list(chain = sub("::subgroup::.*$", "", get("node")),
                             subgroup = sub("^.*::subgroup::", "", get("node")),
                             component = get("component"))]
  part_comp <- merge(
    label_parts[, list(chain = get("chain"), iuis_part = get("iuis_part"),
                       part_node = paste(get("chain"), get("iuis_part"), sep = "::part::"))],
    comp[, list(part_node = get("node"), component = get("component"))],
    by = "part_node", allow.cartesian = TRUE)
  comp_groups <- part_comp[, list(iuis_group_merged = .asc_format_iuis_group(get("iuis_part"))),
                           by = "component"]

  labels <- merge(labels, subgroup_comp[, list(chain = get("chain"), subgroup = get("subgroup"),
                                               component = get("component"))],
                  by = c("chain", "subgroup"), all.x = TRUE)
  labels <- merge(labels, comp_groups, by = "component", all.x = TRUE)
  labels[!is.na(get("iuis_group_merged")), "iuis_group" := get("iuis_group_merged")]
  labels[, c("component", "iuis_group_merged") := NULL]

  labels[, "iuis_group_size" := data.table::uniqueN(get("subgroup")),
         by = c("chain", "iuis_group")]
  labels[, "asc_and_iuis_group" := get("iuis_group")]
  labels[get("iuis_group_size") > 1L,
         "asc_and_iuis_group" := sprintf("%s<sub>%s</sub>", get("iuis_group"),
                                         LETTERS[seq_len(.N)]),
         by = c("chain", "iuis_group")]
  asc_to_iuis <- stats::setNames(labels$iuis_group, labels$subgroup)

  structure(list(labels = labels[],
                 asc_to_iuis = asc_to_iuis,
                 d_allele_to_iuis = d_allele_to_iuis,
                 identity = .asc_iuis_identity(husa, asc_to_iuis, d_allele_to_iuis)),
            class = c("ascIUISVocabulary", "list"))
}

# What each label actually stands for: the alleles, and a digest of their ungapped
# sequences. Reported per chain because the labels are written without their locus
# prefix for V and J, so `V1-2` is not unique across chains on its own.
.asc_iuis_identity <- function(husa, asc_to_iuis, d_allele_to_iuis) {
  x <- data.table::as.data.table(husa)[, list(allele = get("allele"), chain = get("chain"),
                                              gene_type = get("gene_type"), seq = get("seq"))]
  x[, "asc_gene" := alakazam::getGene(get("allele"), strip_d = FALSE, omit_nl = FALSE,
                                      first = FALSE)]
  is_d <- substr(x$allele, 4L, 4L) == "D"
  x[, "iuis_label" := asc_to_iuis[get("asc_gene")]]
  if (length(d_allele_to_iuis)) {
    x[is_d, "iuis_label" := d_allele_to_iuis[get("allele")]]
  }

  dropped <- sum(is.na(x$iuis_label))
  if (dropped) {
    message(sprintf("ascIUISVocabulary: %d of %d alleles carry no IUIS label and are not in `identity`.",
                    dropped, nrow(x)))
  }

  x <- x[!is.na(get("iuis_label"))]
  x[, list(asc_genes = paste(sort(unique(get("asc_gene"))), collapse = ","),
           n_alleles = data.table::uniqueN(get("allele")),
           member_alleles = paste(sort(unique(get("allele"))), collapse = ","),
           member_seq_md5 = .asc_seq_digest(get("seq"))),
    by = c("chain", "gene_type", "iuis_label")]
}

# A fixed-width stand-in for the member sequence set, for cheap comparison. The member
# alleles are already returned, but an allele can be redefined upstream while keeping its
# name, and that is precisely the change a name list cannot show. Gaps are stripped so a
# gapped and an ungapped copy of the same allele hash alike.
.asc_seq_digest <- function(seq) {
  rlang::hash(paste(sort(unique(gsub("[-.]", "", seq))), collapse = "\n"))
}

#' Add IUIS gene-group columns to a genotype-corrected repertoire
#'
#' Translates the ASC gene calls on a repertoire into the IUIS labels built by
#' \code{\link{ascIUISVocabulary}}, adding \code{v_gene_iuis}, \code{j_gene_iuis} and,
#' where the chain has one, \code{d_gene_iuis}.
#'
#' For IGH the D label is per allele, so it is mapped from the allele-level call rather
#' than from the collapsed gene, and a call spanning more than one ASC group names no
#' single gene and stays \code{NA}. V and J map straight from their ASC gene.
#'
#' Rows in and rows out are reported: a call that maps to nothing is left \code{NA} rather
#' than dropped, and the count is messaged so a labelling gap cannot pass unnoticed.
#'
#' @param data A \code{data.frame} of rearrangements carrying the ASC gene calls.
#' @param vocabulary An \code{ascIUISVocabulary} from \code{\link{ascIUISVocabulary}}.
#' @param chain The locus of \code{data}, one of \code{"IGH"}, \code{"IGK"}, \code{"IGL"}.
#' @param v_gene,j_gene,d_gene Column names holding the ASC gene calls. Defaults
#'   \code{"v_gene"}, \code{"j_gene"}, \code{"d_gene"}.
#' @param d_call Column name holding the ASC \emph{allele} call for D, used for IGH only.
#'   Default \code{"d_call_new"}.
#'
#' @return \code{data} as a \code{data.table} with the \code{*_gene_iuis} columns added.
#' @seealso \code{\link{ascIUISVocabulary}}
#' @examples
#' husa <- data.frame(
#'   allele      = c("IGHV1-2*01", "IGHV1-2*02", "IGHD1-CO5H*01"),
#'   chain       = "IGH",
#'   gene_type   = c("IGHV", "IGHV", "IGHD"),
#'   seq         = c("ACGTACGT", "ACGTACGA", "GGTTGGTT"),
#'   iuis_allele = c("IGHV1-2*01", "IGHV1-2*02", "IGHD1-20*01"))
#' voc <- ascIUISVocabulary(husa)
#'
#' rep_dt <- data.frame(
#'   v_gene     = c("IGHV1-2", "IGHV1-2"),
#'   d_gene     = c("IGHD1-CO5H", "IGHD1-CO5H"),
#'   j_gene     = c("IGHJ4", "IGHJ4"),
#'   d_call_new = c("IGHD1-CO5H*01", "IGHD1-CO5H*01"))
#'
#' annotateRepertoireIUIS(rep_dt, voc, chain = "IGH")
#' @export
annotateRepertoireIUIS <- function(data, vocabulary, chain,
                                   v_gene = "v_gene", j_gene = "j_gene", d_gene = "d_gene",
                                   d_call = "d_call_new") {
  if (!inherits(vocabulary, "ascIUISVocabulary")) {
    stop("`vocabulary` must come from ascIUISVocabulary().", call. = FALSE)
  }
  if (!chain %in% c("IGH", "IGK", "IGL")) {
    stop("`chain` must be one of IGH, IGK, IGL; got '", chain, "'.", call. = FALSE)
  }
  x <- data.table::as.data.table(data)
  voc <- vocabulary$asc_to_iuis

  needed <- c(v_gene, j_gene, d_gene, if (chain == "IGH") d_call)
  missing <- setdiff(needed, names(x))
  if (length(missing)) {
    stop("`data` is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }

  # Columns are pulled out by name before any assignment. The argument names here are
  # also column names in a repertoire, and inside `[` a bare name resolves to the column,
  # so `get(v_gene)` would look up the contents of `v_gene` rather than the argument.
  vg <- x[[v_gene]]
  jg <- x[[j_gene]]
  dg <- x[[d_gene]]

  data.table::set(x, j = "v_gene_iuis", value = unname(voc[vg]))
  data.table::set(x, j = "j_gene_iuis", value = unname(voc[jg]))

  if (chain == "IGH") {
    # The D label is a property of the allele, not of the ASC gene, so it is mapped from
    # the allele-level call. A call spanning more than one ASC group names no single gene
    # and stays unlabelled. Mapped once per distinct call rather than once per row.
    dvoc <- vocabulary$d_allele_to_iuis
    dcall <- x[[d_call]]
    lab <- rep(NA_character_, nrow(x))
    single <- !grepl(",", dg) & !is.na(dcall)
    calls <- unique(dcall[single])
    if (length(calls)) {
      mapped <- vapply(calls, function(cc) .asc_map_one(dvoc, cc), character(1))
      lab[single] <- unname(mapped[dcall[single]])
    }
    data.table::set(x, j = "d_gene_iuis", value = lab)
  } else {
    data.table::set(x, j = "d_gene_iuis", value = unname(voc[dg]))
  }

  for (col in c("v_gene_iuis", "j_gene_iuis", "d_gene_iuis")) {
    n_na <- sum(is.na(x[[col]]))
    message(sprintf("annotateRepertoireIUIS [%s]: %s labelled on %d of %d rows (%d unlabelled).",
                    chain, col, nrow(x) - n_na, nrow(x), n_na))
  }
  x[]
}
