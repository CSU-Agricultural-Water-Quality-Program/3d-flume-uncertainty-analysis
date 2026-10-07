# Predictive uncertainty component definitions

These companion figures describe a NEW observed discharge at known, fixed stage.
They do not decompose the existing 95% parameter credible bars for mean rating-curve flow.
Infiltration retains the parameter's terminology; its SD does not estimate a directional loss.

## Calculation

Let theta denote one joint posterior draw and mu = a*H^b, in internal ft/cfs units.
V_curve = Var_posterior(mu).
V_process = E_posterior(sigma_proc^2).
V_infiltration = E_posterior(sigma_Q_add^2).
V_flowmeter = E_posterior[cv_Q^2 * (mu^2 + sigma_proc^2)].
V_total = V_curve + V_process + V_infiltration + V_flowmeter.
The law of total variance integrates over a new latent Q_true; no random simulation is needed.
The cv_Q^2*sigma_proc^2 interaction is allocated to the flowmeter component, where it originates in the likelihood.
Expectations use paired joint draws; products of marginal summaries are not used.
Noise-scale parameter uncertainty is integrated through these expectations, not added as another independent variance term.
Variances are converted to (L/s)^2; each component line is sqrt(V_component) and the combined line is sqrt(V_total).
Relative SD = 100*SD/manufacturer prediction at the same stage. These are not relative biases or 95% intervals.
Stacked shares = 100*V_component/V_total. Shares sum to 100% at every stage.
Neither stage nor manufacturer uncertainty is propagated; panels remain within each specimen's observed stage range.
The allocations describe the fitted model; physical sources are not isolated experimentally.
The flowmeter coefficient is estimated. Its 2% prior median is not a fixed 2% measurement SD.
Posterior mean cv_Q: 1.3037%; root mean square cv_Q: 1.3986%.
The flowmeter line also depends on the fitted mean, its uncertainty, and process scatter, so it need not equal either coefficient summary.

## Provenance and diagnostics

Calibration input: `C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data/New Data Oct 2026/flume_data_oct2026.csv`; MD5 `b6215d73b8dc603ae0bcc7be9beb97c2`.
Saved fit: `C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/output/bayes/oct2026/oct2026_flume_model_3_cmdstanfit.rds`; MD5 `b1ecdd89ac93623262123ea9b4fe94e3`.
8000 joint posterior draws; 504 flume-stage combinations.
Current fit: 744 divergences among 8000 draws; E-BFMI 0.166 to 0.201. Component estimates are provisional until sampling diagnostics are resolved.
Figure captions are exported separately in the companion uncertainty_figure_captions.md file; images contain no captions.
