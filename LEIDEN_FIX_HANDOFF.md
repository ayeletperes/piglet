# Handoff: Leiden reproducibility fix

Branch `fix/leiden-modal-consensus`, two commits, **not pushed**. Pushing a change to a
CRAN package is your call, not mine.

## The defect

`igraph::cluster_leiden` is randomised. The package called it **once, with no seed anywhere
in the clustering path**, so the same germline set could return different ASCs on consecutive
runs.

Measured on `piglet::HVGERM` at resolution 0.955, n = 521:

| | |
|---|---|
| distinct partitions over 10 seeds | **3** |
| cluster count | 128 in every case |
| alleles differing between a single run and the modal answer | **0 to 2 of 521** |
| ARI between a single run and the modal answer | 0.995 to 1.000 |

`scripts/legacy_import/human_asc/infer_asc_usofa.R` ran Leiden 100 times and took the most
frequent membership. That step was lost when the script moved into the package. None of the
three current HUSA runners restore it; all call `inferAlleleClusters`.

## Severity

**Reproducibility, not correctness.** Any released ASC set is within a couple of alleles of
the modal answer and has the right cluster count. Two consecutive HUSA releases
(`results_v3/asc_final/2026-09-09` and `2026-09-10`) show **0 of 1454 shared alleles changing
ASC**, so the wobble has not visibly bitten a published result. It was luck, not design.

## What changed

- `detect_communities_leiden()` gains `seed` (default 1) and `n_runs` (default 1, so direct
  callers are unaffected) and returns the modal partition over `n_runs` seeded runs.
- `.igClust_leiden()`, `igClust()` and `inferAlleleClusters()` thread both, defaulting to
  `n_runs = 100` to match the original procedure.
- `optimize_resolution()` gains `seed` and `n_runs` (default 25) and routes through
  `detect_communities_leiden()`, so the silhouette sweep scores modal partitions rather than
  single draws.
- `.getNClusters()` gains `seed`; both its branches are now seeded.
- `tests/testthat/test-leiden-determinism.R` guards all three levels.

## Verification

- 213 tests pass, 0 failures, including the three new guards.
- `R CMD check`: 0 errors, 0 notes. Two warnings, both from running with
  `--no-build-vignettes`; unrelated to the change.
- Same seed twice through `inferAlleleClusters` now returns an identical table.
- On HVGERM the mode holds **44 of 100** runs across 5 distinct partitions. `n_runs = 25` is
  already enough for stability across starting seeds (2.4 s); 100 costs 8.2 s.

## Decisions left for you

1. **Default `n_runs`.** 100 matches the legacy script and costs ~8 s on 521 alleles. 25 gives
   the same answer here for a third of the time. I defaulted to 100 on the grounds that
   matching the original procedure is worth more than 6 seconds.
2. **Whether to push.** The change is on a branch. A CRAN package's release is your decision.
3. **Whether to regenerate released ASC sets.** Not required on the evidence above, but the
   next release will differ from the last by up to a couple of alleles purely from adopting
   the modal answer, and that is worth a line in the release notes rather than a surprise.

## One thing I did not do

I did not check which path produced each historical release. If you want certainty that a
specific published table came from the legacy script rather than the package, that is a
separate check against the run logs.
