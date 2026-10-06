# Original model 3 fit, October 5, 2026

Recovered from original CmdStan chain CSV files; seed 1815795276.
All 31 parameter summary rows matched the original exported summary.
This is a preserved earlier fit, not a new calibration or an endorsed final fit.
302 divergences, 4888 maximum-treedepth hits; E-BFMI 0.239-0.282.
cmdstanfit.rds stores complete original draws and sampler diagnostics.
Sampling used four chains, 2000 warmup and 2000 sampling iterations per chain,
adapt_delta 0.95, maximum treedepth 10, random initialization range 2,
and CmdStan 2.36.0. Inputs are the September 2026 ft/cfs data.

## Reproduction audit, October 6, 2026

An isolated rerun of the current model 3 call with seed 1815795276 reproduced
every original posterior draw and sampler diagnostic exactly (maximum absolute
draw difference zero). Current manuscript products were not overwritten.
The current data list matches the original input list exactly. The model 3
formula matches both commit 91aac15 and the earlier commit 1785166. The latter
commit's script selected simulated data; the original September fit and this
audit both use the September data. Four chains, 2000 warmup and 2000 sampling
iterations, adapt_delta 0.95, and maximum treedepth 10 match the original run.
The recent model 3 sampling-call edit adds an explicit seed; unit conversions,
plots, diagnostics, and file organization operate outside calibration.

In the original posterior, correlations between log(a) and b are 0.98338,
0.98403, and 0.98144 for flumes 1, 2, and 3. These measured dependencies are
consistent with difficult sampling geometry, but do not establish its sole
cause. Near-one R-hat alone does not resolve the original divergence warnings.
See [Stan diagnostic guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
