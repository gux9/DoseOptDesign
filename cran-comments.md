# CRAN Submission Comments — dosopt 0.1.0

## Test environments

- Local: macOS 14 (Sonoma), R 4.4.1
- win-builder: R-devel, R-release
- rhub: ubuntu-release, windows-x86_64-release, macos-arm64

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

  New submission.

## Downstream dependencies

None. This is a new package.

## Notes for CRAN reviewers

- The Shiny application in `inst/shiny/` requires several packages listed under `Suggests`.
  `run_dosopt_app()` checks for their availability and provides an informative error if missing.
- `simulate_dose_opt()` uses `future.apply` for optional parallelism (also in `Suggests`).
  The function falls back gracefully to sequential computation if unavailable.
- The `\donttest{}` examples in `simulate_dose_opt()` and exact sample size functions are
  wrapped because they can take > 5 seconds with default parameters.
