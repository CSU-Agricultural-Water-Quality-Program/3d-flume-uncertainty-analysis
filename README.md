# 3D-Printed Trapezoidal Flume Calibration and Bayesian Uncertainty Analysis

**Authors (manuscript order):** Emmanuel Deleon, Perry Cabot, A. J. Brown,
Jakob F. Ladow, Daniel C. White, Karen Schlatter, and Erik Wardle.
All authors are affiliated with Colorado State University; specific affiliations
are recorded in `CITATION.cff`.

This repository contains the R analysis, calibration inputs, and saved outputs
for estimating stage-discharge relationships for two 3D-printed trapezoidal
flumes and a fiberglass reference flume. The current calibration uses ten
October 2026 observations per flume, with stage in **meters** and discharge in
**liters per second** throughout.

## Run the analysis

The main script is **[code/bayes_flume.R](code/bayes_flume.R)**. With R,
`rethinking`, CmdStanR/CmdStan, and the plotting/data packages installed, open
it in RStudio and click **Source**, or run from the repository root:

```r
source("code/bayes_flume.R")
```

```powershell
Rscript code/bayes_flume.R
```

Running the script fits the model and exports results to
`output/bayes/oct2026_metric/`. A rerun uses the same output filenames; copy
that folder elsewhere first if you need to preserve a previous run. No fitting
switch or environment variable is required. Supporting functions are imported
automatically. See the [code guide](code/README.md) and
[reproducibility notes](docs/REPRODUCIBILITY.md) for dependencies and saved-fit
regeneration.

Edit priors and uncertainty settings in the labeled `dat` list in the main
script. Its comments identify the point gauge, ModMAG M2000 flow meter, CSU lab
infiltration allowance, rating-curve reference, and process scatter, with units,
sources, and modeling assumptions. The current proportional-SD prior median is
**0.015 (1.5%)**. The meter manual's **0.2%-of-reading component plus 1 mm/s**
is documented separately; the selected empirical prior does not reproduce that
full accuracy envelope or equate it with a manufacturer-reported SD.

The stage bound is +/-0.0015 m, represented by latent stage uniform within the
reported interval. Rating-coefficient prior medians are 941.2 and 2.58, matching
Openchannelflow's published metric equation. The lab's +/-0.02-cfs allowance is
used as an additive-SD prior median of 0.56633693184 L/s, rather than a mean
loss or hard bound. The chosen lognormal prior widths are 0.5.

## Inputs and outputs

The active input is
[`flume_data_oct2026_metric.csv`](<data/New Data Oct 2026/flume_data_oct2026_metric.csv>):

| Column | Meaning | Units |
| --- | --- | --- |
| `depth_m` | Observed stage/head | m |
| `flow_lps` | Observed reference discharge | L/s |
| `flume_type` | `3D1`, `3D2`, or `FIBER` | Identifier |

The original imperial CSV and adjacent conversion-provenance CSV are retained.
Metric conversion preserves the existing observations and their row order.

The current saved native metric run is in
[`output/bayes/oct2026_metric/`](output/bayes/oct2026_metric/). It includes the
fit, parameter summaries, chain/iteration diagnostics, run metadata, rating
curves, point comparisons, uncertainty-component plots, and prediction tables.
The [output guide](output/README.md) distinguishes these from historical fits.

**Diagnostic status:** the included native metric fit has 623 divergences,
7,201 maximum-treedepth hits, and E-BFMI below 0.3 in three of four chains. These
counts come from its
[chain diagnostic CSV](output/bayes/oct2026_metric/oct2026_metric_flume_model_3_chain_diagnostics.csv).
The archived outputs are provisional; this code release does not establish
convergence or publication-ready posterior inference. Diagnostic interpretation
follows the [Stan guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
No calibration was run as part of preparing this release.

## Repository layout

| Location | Purpose |
| --- | --- |
| `code/bayes_flume.R` | Main calibration and result-generation script |
| `code/supporting_code/` | Unit and reporting functions imported by the main script |
| `code/tools/` | Optional input preparation, validation, and saved-fit figure regeneration |
| `code/_unused_old_r/`, `code/_unused_old_python/` | Historical scripts; not the current analysis |
| `data/` | Calibration CSVs, source workbooks, and conversion provenance |
| `docs/` | Reference documents and reproducibility notes |
| `output/` | Current and historical saved fits, diagnostics, figures, and tables |
| `CITATION.cff` | Citation metadata for the complete multi-script analysis |
| `RELEASING.md` | GitHub–Zenodo publication instructions |

`tmp/` and R session history are local preparation artifacts and are excluded
from Git tracking and release archives. Existing local copies are retained.

## Cite this code

Citation metadata are in [CITATION.cff](CITATION.cff). The proposed first release
is **v1.0.0**; its Zenodo DOI is pending. Follow
[RELEASING.md](RELEASING.md) to connect the repository and publish the release.
Once archived, cite the **version-specific Zenodo DOI** for the code used in the
paper. Until then, identify this GitHub repository and the exact commit.
Do not substitute an example DOI for an assigned identifier.

Zenodo can import CFF metadata from a connected GitHub release, and GitHub
provides a citation suggestion from the same file. See
[Zenodo's citation guidance](https://help.zenodo.org/docs/github/describe-software/citation-file/).

## License and references

The existing [GPL version 2 license](LICENSE) is retained. Third-party reference
documents retain their original authorship and notices.

- Openchannelflow. *Large 60-Degree-V Trapezoidal Flume Discharge Table*, pp. 1-2;
  [included reference](docs/large_60-degree-v_trapezoidal_flume_discharge_table.pdf).
- Badger Meter. *ModMAG M2000 User Manual*, MAG-UM-01272-EN-43, August 2026,
  pp. 66 and 75; [manufacturer product/documentation page](https://www.badgermeter.com/products/meters/electromagnetic/m2000-electromagnetic-flow-meter/).
- McElreath, R. (2020). *Statistical Rethinking*, 2nd ed.
  [Book DOI](https://doi.org/10.1201/9780429029608).
- [rethinking source](https://github.com/rmcelreath/rethinking) and
  [CmdStanR installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).
