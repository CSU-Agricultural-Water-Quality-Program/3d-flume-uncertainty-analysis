# Which code to run

**`bayes_flume.R` is the main analysis.** Running it fits the October calibration
and exports the saved fit, diagnostics, figures, prediction tables, and manuscript
parameter table. Stage is in meters and discharge in L/s throughout. It contains
one Bayesian model; there is no fitting switch or environment variable to set.

In RStudio, open `bayes_flume.R` and click **Source**, or run this from the
repository root:

```r
source("code/bayes_flume.R")
```

From PowerShell at the repository root:

```powershell
Rscript code/bayes_flume.R
```

Running selected lines also works when the working directory is the repository
or one of its subfolders. Results go to `output/bayes/oct2026_metric/`.
Edit the labeled `dat <- list(...)` section in `bayes_flume.R` to change prior
medians, prior widths, or the point-gauge accuracy bound. Comments identify the
point gauge, ModMAG M2000 flow meter, lab infiltration allowance, rating-curve
reference, and process scatter, including units, sources, and modeling choices.
Rerun the data/settings section, then the model; supporting functions do not
need to be sourced again. The main script imports the reporting functions
automatically when run from the top.

| Location | Purpose | Run directly? |
| --- | --- | --- |
| `bayes_flume.R` | Fit the model and generate results | Yes |
| `supporting_code/` | Unit and reporting functions used by the main script | No |
| `tools/prepare-flume-metric-data.R` | Prepare the metric copy of the original input CSV | Only when preparing input |
| `tools/validate-flume-metric.R` | Check units, compile the model, and check reporting without sampling | Optional |
| `tools/regenerate-flume-figures.R` | Recreate figures from a matching saved fit | Optional |
| `_unused_old_r/base_flume.R` | Previous main script, including the fitting switch and alternative models | No; historical reference |
| `_unused_old_python/` | Obsolete early Python analyses | No; historical reference |

The result filenames retain `flume_model_3` for compatibility with the existing
reporting helpers. This is the single model fitted in the main script.

The archived `base_flume.R` is unchanged. The new main script preserves its
active likelihood, priors, sampling settings, and reporting functions; models
1, 2, and 4 are omitted. The native metric model has not yet been sampled.
