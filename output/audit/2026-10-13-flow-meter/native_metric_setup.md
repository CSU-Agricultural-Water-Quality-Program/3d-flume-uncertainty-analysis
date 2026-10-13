# Native metric calibration setup

The active R analysis now uses stage in metres and discharge in L/s throughout. The original October CSV is preserved. Its metric copy contains all 30 observations with unchanged row order and flume identifiers; the adjacent provenance CSV records original/metric hashes and exact conversion factors.

The rating priors are independent lognormal distributions with fixed medians 941.2 and 2.58 and fixed log-scale SD 0.5, matching the coefficient values in the supplied Openchannelflow equation Q (L/s) = 941.2 H (m)^2.58. Their hyperparameters are data constants, not fitted parameters. This deliberately specifies priors in the published units rather than preserving the old imperial joint prior exactly.

| Quantity | Fixed prior or uncertainty setting |
| --- | --- |
| Rating coefficient a | LogNormal(log(941.2), 0.5), coefficient units (L/s)/m^b |
| Rating exponent b | LogNormal(log(2.58), 0.5), dimensionless |
| Proportional measurement SD fraction | LogNormal(log(0.002), 0.5); median 0.2% |
| Additive measurement SD | LogNormal(log(0.56633693184), 0.5), L/s |
| Process SD | Exponential(rate 0.0353146667214886), L/s |
| Stage uncertainty | Uniform conditional interval H_obs +/-0.0015 m |

The additive/process priors retain their former physical scales. The additive prior retains the accepted interpretation of the lab allowance as a prior median for an SD, not a strict bound or an estimated mean loss. The empirical error model does not explicitly reproduce the manufacturer's velocity-error specification. The stage model applies the reported bound to computed depth independently across observations.

All prior settings are in `code/flume-metric-units.r`. The commented log-parameterization model remains in the same metric units. Main-model reporting no longer converts fitted draws or predictions; percentage display still converts a dimensionless fraction to percent. Tables use Depth_m; plots and component tables use stage_m; axes and captions use m and L/s. Prediction stages are 0 to 0.24 m at 0.01-m increments.

Validation passed for all 30 input conversions, the active model's C++ compilation, the commented alternate model's Stan semantics, draw-level unit equivalence, all 15 Table 4 summaries, component variance sums, and preservation of historical file hashes. Metre-valued point, curve, and component previews were visually inspected. Those previews used converted historical draws strictly as reporting fixtures, not new native-model results.

No native metric posterior has been sampled. Sampling is disabled by default and requires explicitly setting FLUME_RUN_CALIBRATION=true. A future authorized run writes to `output/bayes/oct2026_metric/`. Existing posterior files and the manuscript were preserved. Validation does not establish convergence; that will require an authorized sampling run.
