# Point-plot definitions

- 30 calibration observations; 8000 posterior draws per observation.
- White-filled points: posterior medians of observation-specific latent true flow Q_true, drawn underneath solid colored measured-flow points.
- Vertical bars: 95% equal-tailed parameter credible intervals of a*H^b, conditional on observed stage; the same estimand as the parameter ribbons.
- The point and bar estimands differ: points summarize Q_true; vertical bars summarize a*H^b and are not credible intervals for Q_true. The parameter mean need not coincide with the latent point. No additional process or measurement scatter is added to the bars.
- Horizontal bars: recorded gauge accuracy bound +/- 0.0015 m; no probability coverage is assigned.
- The fitted calibration integrates uniform stage uncertainty within the recorded gauge bound. Vertical bars evaluate the resulting coefficient draws at observed stage; they are not joint H_true/Q_true intervals.
- Manufacturer: Q (L/s) = 941.2 * H (m)^2.58.
- Relative difference (%) = 100 * (flow - manufacturer) / manufacturer at observed stage. Plotted white points use Q_true; vertical bars use a*H^b. Both parameter-mean and latent medians are retained separately in the CSV. Negative values indicate lower flow; zero-reference rows are NA in the CSV and omitted from the percent plot.
- Relative differences evaluate flow at observed stage. No additional stage uncertainty or manufacturer-curve uncertainty is added to those evaluations.
- All observations retain their measured stages; no jitter, dodge, aggregation, or refitting is applied. Facets share axis limits.
- The manufacturer curve extends from zero to 10% above maximum observed stage.
- Calibration input: `C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data/New Data Oct 2026/flume_data_oct2026_metric.csv`; MD5 `d94cd13bde4b36f88c84a034cd2126c1`.
- Saved fit: `C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/output/bayes/oct2026_metric/oct2026_metric_flume_model_3_cmdstanfit.rds`; MD5 `69be91a6dd4f74f0164294feb1715305`.
