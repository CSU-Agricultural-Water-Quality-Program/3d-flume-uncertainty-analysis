# Release history

## 1.0.0 — prepared for the first Zenodo release

- Established `code/bayes_flume.R` as the main analysis script, with one active
  Bayesian model and automatic fitting when run.
- Calibrated and reported in meters and L/s using the October 2026 ten-point
  observations for each of three flumes; retained the original input CSV and
  metric-conversion provenance.
- Kept editable priors in the main script, labeled by device, source, units,
  and modeling assumptions. The release's proportional-SD prior median is
  1.5% (0.015); the manual's 0.2% component remains a separately cited
  specification, not the selected prior median.
- Organized reporting functions into `code/supporting_code/` and optional
  preparation/validation/regeneration commands into `code/tools/`.
- Preserved earlier R/Python scripts and historical calibration outputs.
- Added citation metadata, reproducibility documentation, and Zenodo release
  instructions. Excluded local scratch files and R session history from Git
  tracking and release archives without deleting local copies.

The included native metric fit began on 2026-10-13. Its saved diagnostics report
623 divergences, 7,201 maximum-treedepth hits, and E-BFMI below 0.3 in three
chains. This release preserves those results with their diagnostics; it does
not certify convergence. No R code was executed during release preparation.
