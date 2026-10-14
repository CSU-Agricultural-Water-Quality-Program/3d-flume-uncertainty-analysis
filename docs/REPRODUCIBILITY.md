# Reproducing the analysis

## Entry point and data

Run `code/bayes_flume.R`, which imports the four reporting/unit modules and the
path helper in `code/supporting_code/`. The code guide identifies optional tools
and historical scripts. All editable calibration settings are in the main
script's `dat` list. Run that entire section after changing settings, including
the stage-bound calculations, before fitting again.

The active metric CSV contains 30 observations, ten each for `3D1`, `3D2`, and
`FIBER`. Its MD5 is `d94cd13bde4b36f88c84a034cd2126c1`, matching the current
saved fit's run metadata. Its adjacent provenance CSV records conversion from
the original imperial observations using 0.3048 m/ft and 28.316846592 L/s per cfs.
No new measurements are inferred by conversion.

Input paths saved in provenance/run metadata identify the original workstation.
Use the repository-relative paths in the scripts on another computer.

## Software environment

The saved native metric run records `rethinking` 2.42, CmdStanR 0.9.0.9000, and
CmdStan 2.36.0. The table supplements those records with installed-package
DESCRIPTION files inspected during release preparation, without loading R.

| Component | Version | Source of version information |
| --- | --- | --- |
| `rethinking` | 2.42 | Saved run metadata and installed DESCRIPTION |
| `cmdstanr` | 0.9.0.9000 | Saved run metadata and installed DESCRIPTION |
| CmdStan | 2.36.0 | Saved run metadata |
| `posterior` | 1.6.1 | Installed DESCRIPTION |
| `ggplot2` | 3.5.1 | Installed DESCRIPTION |
| `dplyr` | 1.1.4 | Installed DESCRIPTION |
| `tidyr` | 1.3.1 | Installed DESCRIPTION |

The installed GitHub package sources are `rmcelreath/rethinking` at commit
`ac1b3b2cda83f3e14096e2d997a6e30ad109eeee` and `stan-dev/cmdstanr` at commit
`5f2a8c850890e54e6fccc8108e5a31cb08010ab8` (from DESCRIPTION `RemoteSha` fields).
These distinguish the development CmdStanR build from a numbered release.
This is an environment record, not a complete dependency lock or a verified
clean-machine installation. The saved metadata do not record the R version.

Install R and the listed packages plus their dependencies. CmdStan requires a
working C++ toolchain; use the
[official CmdStanR installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
and the [rethinking installation instructions](https://github.com/rmcelreath/rethinking).
No packages, compiler, or model were installed or executed during release preparation.

## Saved run and diagnostic status

`output/bayes/oct2026_metric/` contains the current native metric fit. The run
started at `2026-10-13 22:18:38 UTC` (16:18:38 MDT), with seed 1815795276, four
chains, 2,000 warmup and 2,000 retained iterations per chain, `adapt_delta=0.95`,
and maximum treedepth 10, as recorded in its metadata.

Its chain diagnostics record 623 divergences, 7,201 treedepth-limit hits, and
E-BFMI below 0.3 in three chains. These warnings require resolution before
treating posterior summaries as reliable manuscript results; see
[Stan's diagnostic guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
The release preserves the fit and diagnostics; it does not certify convergence.
We did not rerun sampling, compilation, or reporting for release preparation.

All numerical model settings were left unchanged during this preparation.
The current empirical proportional-SD prior median is 0.015 (1.5%), not the
manual's 0.002 component. The latter remains cited as an instrument specification.
The additive prior median is 0.56633693184 L/s, the stage bound is 0.0015 m,
and the rating prior medians are 941.2 and 2.58.

## Optional commands

From the repository root:

```powershell
# Prepare the metric CSV from the original observations.
Rscript code/tools/prepare-flume-metric-data.R

# Check data, compile without sampling, and preview reporting with historical draws.
Rscript code/tools/validate-flume-metric.R

# Regenerate figures from the included native metric fit; no sampling.
Rscript code/tools/regenerate-flume-figures.R
```

The validator's previews use transformed historical draws as reporting fixtures,
not inference from a newly fitted model. Figure regeneration uses the matched
fit/input hashes and metadata. Both tools write files; copy existing products
first if you need to preserve them exactly. Historical fits and archived code
are not substitutes for the current native metric likelihood and priors.
