# Point-plot definitions

- 30 calibration observations; 8000 posterior draws per observation.
- White-filled points: posterior medians of observation-specific latent true flow Q_true, drawn underneath solid colored measured-flow points.
- Vertical bars: 95% equal-tailed parameter credible intervals of a*H^b, conditional on observed stage; the same estimand as the parameter ribbons.
- The point and bar estimands differ: points summarize Q_true; vertical bars summarize a*H^b and are not credible intervals for Q_true. The parameter mean need not coincide with the latent point. No additional process or measurement scatter is added to the bars.
- Horizontal bars: recorded gauge accuracy +/- 0.000328084 ft (+/- 0.01 cm); no probability coverage is assigned.
- The current likelihood fixes stage. Horizontal accuracy bars are not propagated into vertical credible intervals.
- Manufacturer: Q (cfs) = 1.55 * H (ft)^2.58; converted to cm and L/s.
- Relative difference (%) = 100 * (flow - manufacturer) / manufacturer at observed stage. Plotted white points use Q_true; vertical bars use a*H^b. Both parameter-mean and latent medians are retained separately in the CSV. Negative values indicate lower flow; zero-reference rows are NA in the CSV and omitted from the percent plot.
- Gauge accuracy and manufacturer uncertainty are not propagated into relative differences.
- All observations retain their measured stages; no jitter, dodge, aggregation, or refitting is applied. Facets share axis limits.
- The manufacturer curve extends from zero to 10% above maximum observed stage.
- Calibration input: `C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data/New Data Oct 2026/flume_data_oct2026.csv`; MD5 `b6215d73b8dc603ae0bcc7be9beb97c2`.
- Saved fit: `C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/output/bayes/oct2026/oct2026_flume_model_3_cmdstanfit.rds`; MD5 `b1ecdd89ac93623262123ea9b4fe94e3`.
