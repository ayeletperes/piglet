library(markr)
library(piglet)

# Directories
pkg_path <- "."
doc_path <- "./docs"

# Regenerate the API reference (topics), the vignettes and the release notes
# from the package sources, WITHOUT touching:
#   * docs/index.md  -- the hand-maintained "choose your path" home page
#   * mkdocs.yml     -- the hand-maintained navigation
#
# build_mkdocs() would also call build_md_index() (overwriting index.md from the
# README) and, with yaml=TRUE, build_yaml() (overwriting the nav). We therefore
# call the individual builders directly and deliberately skip those two steps.
#
# When a new exported function is added, its topic page is generated here but
# must be added to the "API Reference" section of mkdocs.yml by hand.
pkg <- as.sd_package(pkg_path, site_path = doc_path, mathjax = FALSE)
pkg$topics    <- build_md_topics(pkg, style = "mkdocs")
pkg$vignettes <- build_md_vignettes(pkg)
pkg$news      <- build_md_news(pkg, style = "mkdocs")

# markr renders \email{x} as a bare Markdown link "[x](x)"; mkdocs --strict then
# treats the address as a broken internal link. Rewrite such email links to use
# the mailto: scheme across the generated topic pages.
for (tf in list.files(file.path(doc_path, "topics"), pattern = "\\.md$", full.names = TRUE)) {
  txt <- readLines(tf, warn = FALSE)
  fixed <- gsub("\\]\\((?!mailto:|https?://|\\./|/)([^()\\s]+@[^()\\s]+)\\)",
                "](mailto:\\1)", txt, perl = TRUE)
  if (!identical(txt, fixed)) writeLines(fixed, tf)
}

# markr can fail to strip the YAML front matter from a vignette whose Rmd has a
# complex header (e.g. a multi-line `vignette: >` block); it then leaves the raw
# YAML as visible text right after the H1 title ("Error finding yaml block").
# Remove such a leftover block: when line 2 of a generated vignette page is a
# "---" fence, drop everything up to and including the matching closing fence.
for (vf in list.files(file.path(doc_path, "vignettes"), pattern = "\\.md$", full.names = TRUE)) {
  lines <- readLines(vf, warn = FALSE)
  if (length(lines) >= 2L && grepl("^# ", lines[1]) && lines[2] == "---") {
    fences <- which(lines == "---")
    close_idx <- fences[fences > 2L][1]
    if (!is.na(close_idx)) {
      writeLines(c(lines[1], lines[(close_idx + 1):length(lines)]), vf)
    }
  }
}
