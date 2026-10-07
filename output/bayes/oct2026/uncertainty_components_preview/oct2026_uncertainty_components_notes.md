# Uncertainty component previews

These are companion figures for a NEW observed discharge at a known, fixed stage.
They are not a decomposition of the existing 95% credible bars for mean rating-curve flow.
Infiltration retains the user's parameter terminology; its SD does not estimate a directional loss.

## Definitions

Let theta denote one joint posterior draw and mu = a*H^b, in internal ft/cfs units.
V_curve = Var_posterior(mu).
V_process = E_posterior(sigma_proc^2).
V_infiltration = E_posterior(sigma_Q_add^2).
V_flowmeter = E_posterior[cv_Q^2 * (mu^2 + sigma_proc^2)].
V_total = V_curve + V_process + V_infiltration + V_flowmeter.
This follows the law of total variance, integrating over a new latent Q_true.
The cv_Q^2*sigma_proc^2 interaction is allocated to the flowmeter component, where it originates in the likelihood.
Expectations use paired joint draws; no products of marginal summaries are used.
Noise-scale parameter uncertainty is integrated through these expectations, not treated as an extra independent variance term.
Variances are converted to (L/s)^2 before computing their square roots.
Each colored SD line is sqrt(V_component); the combined line is sqrt(V_total). SD lines must not be added.
Relative SD = 100*SD/manufacturer prediction at the same observed stage; it is not a relative bias or a 95% interval.
The stacked figure shows 100*V_component/V_total, summing to 100% at every stage.
No stage uncertainty or manufacturer uncertainty is propagated; each panel stays within that specimen's observed stage range.
These are estimated model allocations, not experimentally isolated causal contributions.

## Draft captions

Uncertainty scales for a new discharge measurement at fixed stage, separated into infiltration, proportional flowmeter, process, and rating-curve parameter contributions. The black line is the combined predictive standard deviation. Upper panels report L/s; lower panels normalize by the manufacturer prediction at the same stage. Colored lines are square roots of posterior variance contributions and do not add linearly. These are standard deviations, not 95% credible interval half-widths.

Estimated shares of predictive discharge variance attributable to each modeled component at fixed stage. Shares sum to 100%; a declining share need not imply declining absolute variance. The flowmeter contribution includes its interaction with process variation. The allocations reflect the specified model and fitted posterior.

## Provenance and limitations

Input MD5: b6215d73b8dc603ae0bcc7be9beb97c2
Saved fit MD5: b1ecdd89ac93623262123ea9b4fe94e3; 8000 joint posterior draws.
Current fit has 744 divergent transitions out of 8000 draws; E-BFMI is 0.166 to 0.201. These previews are provisional until sampling diagnostics are resolved.
No fitting was performed; existing calibration figures, fits, tables, and diagnostics were preserved.
