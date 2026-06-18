library(markr)
library(piglet)

# Directories
pkg_path <- "."
doc_path <- "./docs"

# Build
# NOTE: yaml=FALSE so the hand-maintained mkdocs.yml nav (Getting Started /
# Lessons / Core Concepts / Task Guides / API Reference) is preserved. markr only
# regenerates the topic and vignette markdown under docs/, not the navigation.
build_mkdocs(pkg_path, doc_path=doc_path, yaml=F)
