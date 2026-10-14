# Analysis outputs

## Native metric model

The active R workflow now calibrates and reports in **m and L/s**. It reads
`data/New Data Oct 2026/flume_data_oct2026_metric.csv`, with `depth_m`, `flow_lps`,
and `flume_type`. The original CSV remains unchanged, and the adjacent provenance
CSV records both file hashes and exact unit factors. The manufacturer reference
is `Q (L/s) = 941.2 * H (m)^2.58`. These coefficient values are fixed prior
medians, with the retained log-scale SD of 0.5. Prior settings are data constants,
not estimated hyperparameters.

The saved native metric fit is in `bayes/oct2026_metric/`, with its run metadata
and chain diagnostics. Running the main script fits again into this folder.
The included run has 623 divergences, 7,201 maximum-treedepth hits, and E-BFMI
below 0.3 in three chains; its posterior outputs remain provisional. Tables
report `Depth_m`; point and component data
report `stage_m`; stage axes and captions report m. Prediction grids retain
the previous physical extent at 0.01-m intervals from 0 to 0.24 m. Reporting
does not convert native posterior draws again. Additive/process SDs are L/s;
the coefficient `a` has units `(L/s)/m^b`, while `b` and proportional fractions
are dimensionless. Displaying `cv_Q` as percent still multiplies its fraction by 100.

The stage bound is +/-0.0015 m and enters the calibration through latent stage.
The empirical proportional-SD prior median is 0.015 (1.5%). The additive-SD prior
median is 0.56633693184 L/s, retaining the existing interpretation of the lab
allowance. The process prior rate is 0.0353146667214886 per L/s, retaining the
physical scale of the former imperial exponential prior. These error assumptions
do not reproduce the complete manufacturer's velocity-error specification.

`code/tools/validate-flume-metric.R` checks conversion, active-model compilation,
native table/figure values,
and preservation of historical products without sampling. Its temporary
previews use converted historical draws solely as reporting fixtures.
Saved-fit regeneration defaults to the included native metric fit. It can also
inspect historical imperial fits with explicit arguments and a separate output
folder; those draws are converted at the driver boundary, preserving paired
scale/exponent draws. Historical fits cannot be relabeled as native metric fits.

## October 2026 point comparisons

The saved October calibration outputs are in `bayes/oct2026/` with the
`oct2026_` prefix. The existing files below describe the earlier fixed-stage fit
with ft/cfs calibration and cm/L/s reporting. September settings are retained
in historical scripts and run metadata.
The existing October files have also been renamed to match. The figures use
the saved model 3 posterior without refitting or changing the model or tables:

- `bayes/oct2026/oct2026_flume_points_with_xy_errors.{png,jpg}`: calibration
  points and horizontal/vertical bars around a solid black manufacturer curve.
  Separate flume panels share axes because all three flumes have identical
  measured coordinates in the October CSV.
- `bayes/oct2026/oct2026_flume_points_relative_difference.{png,jpg}`: the same
  points and bars expressed as percent relative difference, using
  `100 * (flow - manufacturer) / manufacturer` at each observed stage.
  Zero indicates agreement; negative values indicate lower flow.
- `bayes/oct2026/oct2026_flume_point_plot_data.csv`: unrounded point summaries,
  interval endpoints, observation indices, and manufacturer-relative differences.
- `bayes/oct2026/oct2026_flume_point_plot_notes.md`: definitions and input/fit hashes.
- `bayes/oct2026/oct2026_flume_figure_captions.md`: draft manuscript captions
  for both point figures and both curve figures. Captions are exported separately;
  no captions are embedded in the images. Manufacturer references appear in legends.

White-filled points are posterior medians of observation-specific latent true
flow `Q_true`, labelled "Estimated latent true flow" and drawn underneath the
smaller solid colored measured-flow points. Vertical bars are 95% equal-tailed
parameter credible intervals from all 8,000 paired `a,b` draws at observed
stage. These summarize the same quantity as the original parameter ribbons;
they exclude additional process and measurement scatter. The parameter
posterior still conditions on the fitted measurement likelihood. Points and
vertical bars describe different estimands: the bars are intervals for mean
rating-curve flow, not intervals for `Q_true`, and their centers need not coincide
with the white points. Both latent and parameter-mean summaries are retained
separately in the CSV. Horizontal bars show the
reported October gauge accuracy of +/- 0.01 cm, without assigning a
coverage probability. The model treats stage as fixed: gauge accuracy is
displayed separately and is not propagated into these vertical intervals.
Manufacturer-relative differences also condition on the observed stage and
the fixed manufacturer equation. A zero manufacturer prediction makes a
relative difference undefined; those values are exported as NA and omitted
from the percent plot. Coordinates are not jittered or aggregated.

The October combined and faceted curve plots also use white latent true-flow
points and colored measured-flow points. Their point intervals now use the
same parameter summaries as the ribbons and the point-only figures. All
8,000 draws are used so the interval definitions and endpoints match.

## Predictive uncertainty components

`code/bayes_flume.R` sources `code/supporting_code/flume-uncertainty-components.R` and calls its
export function using the same joint posterior draws as the other figures.
Running the main script automatically produces these permanent outputs:

- `bayes/oct2026/oct2026_flume_uncertainty_components_sd.{png,jpg}`: component
  standard-deviation scales in L/s and as percentages of the manufacturer prediction.
- `bayes/oct2026/oct2026_flume_uncertainty_variance_shares.{png,jpg}`: stacked
  shares of total predictive variance, summing to 100% at each stage.
- `bayes/oct2026/oct2026_flume_uncertainty_components_data.csv`: all component
  variances, standard deviations, normalized scales, and variance shares.
- `bayes/oct2026/oct2026_flume_uncertainty_components_notes.md`: exact calculation,
  input/fit hashes, flowmeter coefficient summaries, and current diagnostic limitations.
- `bayes/oct2026/oct2026_flume_uncertainty_figure_captions.md`: manuscript captions
  for both component figures. No captions are embedded in the images.

These figures describe a new observed discharge at fixed stage, including mean
rating-curve parameter uncertainty, process variation, infiltration uncertainty,
and proportional flowmeter uncertainty. They do not decompose the existing 95%
credible bars for mean flow. Components add as variances; standard deviations
combine by root-sum-of-squares. The proportional flowmeter contribution integrates
over new latent flow and includes `E[cv_Q^2 * sigma_proc^2]` as well as
`E[cv_Q^2 * mu^2]`. The exact definitions are exported in the notes. Each panel
stays within that flume's observed stage range. The coefficient `cv_Q` is estimated;
its 2% prior median does not fix the fitted line at 2%. These are model allocations,
not experimentally isolated physical causes; resolve sampling warnings before
final inference.

The earlier `uncertainty_components_preview/` folder is retained as a historical
preview. Current figures and captions are exported directly to `bayes/oct2026/`.
The old `code/_unused_old_r/preview-flume-uncertainty-components.R` command is retained as a
compatibility driver and now calls the same helper for the permanent outputs.

To regenerate all six point, curve, and component figures from the saved October fit, run from the
repository root:

```powershell
Rscript code/tools/regenerate-flume-figures.R
```

The driver verifies the input hash against the fit's run metadata, validates
latent medians against the `Q_true` draws and parameter centers/endpoints
against the paired `a,b` draws, checks
percent summaries against draw-level relative differences, and checks equality
with curve-plot point intervals. The fit, diagnostics, and tables are protected
against changes. The driver also accepts three positional arguments:
`input_csv fit_rds output_dir`. The normal `bayes_flume.R` workflow calls the same
plot helpers after exporting the curve figures. Component variance sums and
stacked shares are also checked. Subsequent fits use the explicit
`output_prefix` (currently `oct2026_metric`) for every export; switch the input,
output folder, and prefix together to analyze a different batch.

## Existing products

| Folder | Contents | Generating code |
| --- | --- | --- |
| `bayes/sep2026/` | September 2026 model 3 calibration: metric tables, figures, dashboard, Stan diagnostic CSVs, run metadata, and saved fit | Historical Bayes runs |
| `bayes/legacy/` | Earlier Bayes tables, figures, and diagnostics retained for reference | Historical Bayes runs |
| `python_monte_carlo/` | Monte Carlo refitting plot | `code/_unused_old_python/flume-montecarlo.py` |
| `python_measurement_error/` | Error-bar and propagated measurement-band plots | `code/_unused_old_python/flume-analysis.py`, `code/_unused_old_python/flume-analysis-with-band.py` |
| `simulation/` | True parameters used to generate simulated calibration data | `code/_unused_old_python/simulate_flume_data_model3.py` |
| `manuscript/` | Manually assembled manuscript summary workbook | Manual synthesis |

The September manuscript products remain in `bayes/sep2026/`. Model 4 is disabled and
its output files have been removed. Read that folder's `sep2026_export_notes.md`
and chain diagnostic CSV before interpreting a rerun.

The September stage-flow comparison is saved as
`bayes/sep2026/sep2026_flume_posterior_curves_final.jpg`. It uses distinct curve
patterns, colored circles for estimated flow, white-filled circles for measured flow, and
12% opacity for uncertainty ribbons. The companion
`bayes/sep2026/sep2026_flume_posterior_curves_by_flume.png` separates the three
flumes into panels with identical axis limits to reveal overlapping uncertainty.
Both retain the existing 95% intervals and zero-to-extrapolation stage range.

The prediction tables `sep2026_flume_generated_display.csv` and
`sep2026_flume_generated_numeric.csv` evaluate discharge at integer depths
0-24 cm, in 1-cm increments. The display table contains posterior medians and
95% equal-tailed intervals; its numeric companion retains unrounded discharge.
Predictions are recomputed at these depths, rather than relabelling converted
feet as rounded centimetres.

`bayes/sep2026/sep2026_manuscript_table4.csv` follows manuscript Table 4's five
parameter rows and three flume columns, with posterior means and 95% HPDIs
computed from all draws. The scale coefficient uses stage in cm and flow in
L/s; standard deviations are L/s, b is dimensionless, and cv is percent of
reading. The shared cv result is repeated across the three flume columns.
The formatted CSV uses UTF-8 with a BOM for Excel; its numeric companion
`sep2026_manuscript_table4_numeric.csv` retains full precision, units, interval
probability, draw count, and shared versus flume-specific scope. The additive
term is labelled a measurement-error standard deviation, consistent with its
role in the likelihood rather than a directional loss or mean offset.

The Python driver retains `data/flume_data.csv` as its default input and reports
in its original ft/cfs units. It does not use the September dataset unless it is
passed explicitly with `--csv`. Output folders identify analysis type; use run
metadata to identify the Bayes dataset and fit.

The Python scripts are now archived in `code/_unused_old_python/`. They are
obsolete for the final analysis and retain their historical assumptions.
`code/bayes_flume.R` is the main Bayesian workflow. Running it fits the model
and exports results without an environment-variable switch. Supporting functions
are imported automatically from `code/supporting_code/`. Preparation, validation,
and saved-fit regeneration are optional commands in `code/tools/`.
