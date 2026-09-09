# Working on piglet

## Branches

`dev` is the integration branch. **Start new branches from `dev`, not `master`.**

```bash
git switch dev && git pull
git switch -c feature/my-thing
```

`master` is the release branch: `dev` merges into it when a version is ready.
Merge feature branches into `dev` with `--no-ff`, so a change that turns out wrong
can be reverted as one commit rather than unpicked.

## Before you push

`devtools::check()` runs the vignettes and every example, which takes ten minutes or
so. That is the gate before a release, not the one to run while iterating. For a fast
pass that still catches the things that actually break:

```r
devtools::document()                       # regenerate NAMESPACE and man/ first
devtools::test()                           # the test suite, about a minute
devtools::check(args = c("--no-manual"), build_args = "--no-build-vignettes")
```

Run the full `devtools::check()` before merging `dev` into `master`. It must be
0 errors, 0 warnings, 0 notes: this package goes to CRAN, and a note that is
tolerated locally is a rejection there.

Regenerate documentation whenever roxygen comments change. A stale `.Rd` or a
missing `NAMESPACE` export is the most common reason a branch fails CI having passed
locally, because `load_all()` sees functions that an installed package would not.

## CI

`.github/workflows/R-CMD-check.yaml` runs `R CMD check` on macOS, Windows and three
Linux R versions, on every push to `master` or `dev` and on every pull request into
either. So a branch is checked when it merges into `dev`, not only at the release
merge. If you would rather not spend a full matrix on routine `dev` pushes, trim the
`push:` branch list and leave `pull_request:` as it is; the pull request is the gate
that matters.

## Changes that move published numbers

Some of this package reproduces analyses that have already been published. If a change
alters a result rather than the code around it, say so in the commit message with the
size of the difference, and check it against the run it is meant to reproduce rather
than against your expectation. "The tests still pass" is not the same claim.
