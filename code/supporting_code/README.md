# Supporting code

These modules are imported by `../bayes_flume.R`. Loading them defines unit
constants or functions; it does not fit a model or generate results. Calibration
priors and uncertainty settings are edited in the main script's `dat` list.

- `flume-metric-units.R`: unit constants, unit assertions, and conversion of historical imperial draws; no calibration priors.
- `flume-reporting.R`: parameter summaries, sampler diagnostics, rating-curve figures, and prediction/manuscript tables.
- `flume-point-plots.R`: calibration point figures, manufacturer comparisons, and their captions.
- `flume-uncertainty-components.R`: predictive variance components and their figures/captions.
- `project-paths.R`: locate the repository without relying on the current working directory.

Run the main entry point rather than these files individually.
