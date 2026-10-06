# Analysis outputs

| Folder | Contents | Generating code |
| --- | --- | --- |
| `bayes/sep2026/` | Current September 2026 model 3 calibration: metric tables, figures, dashboard, Stan diagnostic CSVs, run metadata, and saved fit | `code/bayesflume.r` |
| `bayes/legacy/` | Earlier Bayes tables, figures, and diagnostics retained for reference | Historical Bayes runs |
| `python_monte_carlo/` | Monte Carlo refitting plot | `code/flume-montecarlo.py` |
| `python_measurement_error/` | Error-bar and propagated measurement-band plots | `code/flume-analysis.py`, `code/flume-analysis-with-band.py` |
| `simulation/` | True parameters used to generate simulated calibration data | `code/simulate_flume_data_model3.py` |
| `manuscript/` | Manually assembled manuscript summary workbook | Manual synthesis |

The current manuscript analysis is in `bayes/sep2026/`. Model 4 is disabled and
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

Run `python code/main.py` from the repository root to regenerate the three Python
plots. Source `code/bayesflume.r` in R to refit model 3 and regenerate the current
Bayes products. The simulation generator writes `data/sim_data.csv`; it is a
separate workflow and is not run by the Python plot driver.
