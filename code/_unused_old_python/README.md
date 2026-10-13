# Unused old Python code

These scripts belong to the initial analysis and are obsolete for the current
Bayesian calibration. They retain historical ft/cfs inputs and uncertainty
assumptions. The final workflow is `../bayes_flume.R`.

| File | Historical purpose |
| --- | --- |
| `main.py` | Driver for the early Python figures |
| `flume-analysis.py` | Power-law fits with measurement error bars |
| `flume-analysis-with-band.py` | Measurement-band comparison figure |
| `flume-montecarlo.py` | Monte Carlo perturbation and refitting |
| `simulate_flume_data_model3.py` | Simulated data for the earlier Model 3 specification |

Files are preserved for provenance. Their paths were adjusted for this folder;
their statistical assumptions were not updated. They are not imported by the
main R workflow. Running them can recreate or overwrite historical products.
