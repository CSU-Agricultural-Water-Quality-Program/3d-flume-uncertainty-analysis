# September 2026 metric exports

Stage is reported in cm and discharge in L/s. Calibration inputs, priors, and
likelihoods remain in ft and cfs. The exact reporting factors are 30.48 cm/ft
and 28.316846592 L/s per cfs, derived from the international foot.
Source: [NIST SP 811, Appendix B.8](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8).

Model 3 outputs use the original four calibration chains from the RStudio run
started on October 5, 2026 at 14:05 America/Denver. Before regeneration, all 31
rows and six numerical columns of the original parameter-summary CSV matched
the recovered chains, and every Stan data input matched the September data.
The metric coefficient was transformed within paired posterior draws:
`a_cm = a * 28.316846592 / 30.48^b`. The original posterior was retained.

The stage-discharge figure evaluates all fitted curves and the manufacturer
equation from 0 to 18.507456 cm, 10% above the maximum observed stage. Observed
points retain their measured stages. The prediction tables retain the previous
reporting stages in cm and add zero; flow columns, intervals, and manufacturer
predictions are in L/s. The parameter summaries retain their original 89%
intervals; the rating-curve figure and prediction tables retain 95% intervals.

Dashboard rendering corrects the positional diagnostic-column assumptions in
the installed rethinking 2.42 function: `divergent__` and `energy__` are selected
by name and supplied in the expected positions. The fitted backend is unchanged.
Source: [rethinking dashboard source](https://github.com/rmcelreath/rethinking/blob/master/R/dashboard.r).

The original model 3 fit has 302 divergent transitions among 8,000 post-warmup
draws, 4,888 maximum-treedepth hits, and E-BFMI below 0.3 in all four chains.
Model 4 exports come from a verification rerun of the unchanged model, with
2,985 divergent transitions among 8,000 post-warmup draws. Its summaries are
verification outputs requiring diagnostic review. These sampling issues require
assessment before treating either model's estimates as reliable final inference.
Source: [Stan diagnostic guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
