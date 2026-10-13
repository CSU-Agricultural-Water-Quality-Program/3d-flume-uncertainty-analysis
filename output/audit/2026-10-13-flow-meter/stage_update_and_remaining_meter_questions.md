# Stage uncertainty update and remaining ModMAG questions

This follow-up records the user's decisions after the initial specification audit. The analysis now concerns only the October 2026 ten-point calibration and its ModMAG M2000 reference meter. Parameterization and metric conversion of the entire model remain later steps. The lab's reported infiltration allowance of +/-0.02 cfs is accepted as an input magnitude; its mathematical interpretation is explained below rather than rederived.

## Stage change implemented

1.5 mm equals 0.15 cm, 0.0015 m, or 0.00492126 ft. The proposed 0.015-cm conversion was a factor of ten too small. The earlier 0.005-ft gauge value equals 1.524 mm and is approximately consistent with 1.5 mm; it is unrelated to a flow-meter quantity reported in cfs.

The active `code/bayesflume.r` now uses the conditional uncertainty model

H_true,i | H_obs,i ~ Uniform(H_obs,i - delta_H, H_obs,i + delta_H),

with delta_H = 0.15 cm, converted to feet for the existing calibration code. The rating relation uses H_true rather than H_obs. This treats the reported +/- accuracy as a bound and assumes equal probability within it, independent errors across observations, and a locally flat stage prior. It applies the allowance to computed depth, not separately to the water-surface and bottom readings. A shared uncertain datum is not introduced. These assumptions are explicit and revisable if more information is provided.

Latent stages have vector bounds matching their respective observed-depth intervals. Metric parameter exports convert H_true to cm. Point-plot and component notes distinguish uncertainty incorporated during calibration from predictions evaluated at a new known stage. Saved-fit regeneration reads the stage bound from new run metadata while preserving the previous fixed-stage fit's definitions.

The revised script writes to `output/bayes/oct2026_stage_error/` when it is eventually run. The original October model and results were not overwritten. Validation covered R syntax, all 30 positive stage intervals, generated Stan semantics, old-fit metadata compatibility, and revised figure-caption wording. No posterior sampling or convergence validation has been performed for the revision. Meter and infiltration priors have not yet changed.

## Metric flume reference verified

The supplied [Openchannelflow discharge table](<D:/OneDrive - Colostate/Cabot,Perry's files - 3D Printing Flumes/REFERENCES/large_60-degree-v_trapezoidal_flume_discharge_table.pdf>) states, on both pages,

Q (L/s) = 941.2 [H (m)]^2.58.

Page 1 was also visually inspected. In centimetres this is Q (L/s) = 0.0065115131 [H (cm)]^2.58. The coefficient 941.2 differs by approximately 0.0147% from the exact transformation of its separately rounded imperial coefficient 1.55. For a metric-native reference, use the documented 941.2 and specify metres explicitly. This document does not by itself determine the probability distribution or dependence of the priors on fitted scale and exponent.

## Infiltration allowance and mean loss

The accepted +/-0.02-cfs allowance is +/-0.5663369 L/s. An additive SD in the observation likelihood represents uncertainty around its mean, including uncertainty attributable to infiltration. It does not estimate how much discharge was lost. A directional loss would enter the mean through a mass balance between locations. For example, if upstream inflow is 10 L/s and known leakage is 0.2 L/s, expected downstream flow is 9.8 L/s. A symmetric error term centered at 10 L/s only broadens uncertainty around 10 L/s.

Retaining the lab allowance as an infiltration-related uncertainty contribution is consistent with the user's interpretation. What remains a modeling choice is whether +/-0.02 cfs denotes a bound or an SD. If used as a uniform half-width, its SD is 0.3269748 L/s. If a lognormal prior remains centered at 0.5663369 L/s, that states that the SD has that prior median; it does not encode a strict +/-0.5663369-L/s bound. The manual cannot select the prior width for this facility allowance. Calling the term infiltration-related additive uncertainty avoids claiming that the model estimates a mean loss.

## A flow-dependent standard deviation

The meter function supplies the scale of a probability distribution; it does not replace the standard deviation concept. The current code already has a flow-dependent SD through cv_Q multiplied by Q_true. The manufacturer's full expression adds a velocity-related absolute contribution.

For provisional DN100 internal diameter and flow in L/s, the ModMAG specification half-width is E(Q) = 0.002 |Q| + 0.00785398 L/s. If symmetric uniform bounds are assumed, u_meter(Q) = E(Q)/sqrt(3). This can then enter a Normal observation model as an explicit moment approximation. For independent infiltration-related scatter with SD sigma_add, the measurement SD would be sqrt(u_meter(Q)^2 + sigma_add^2). The two components of the manufacturer's stated envelope add before division by sqrt(3). Actual diameter and certificate interpretation remain to be confirmed.

References: [ModMAG manual](<D:/OneDrive - Colostate/Cabot,Perry's files - 3D Printing Flumes/ModMag Manual.pdf>), pp. 39 and 66; [Badger Meter product specifications](https://www.badgermeter.com/products/meters/electromagnetic/m2000-electromagnetic-flow-meter/); [JCGM 100:2008, 4.3.7](https://www.iso.org/sites/JCGM/GUM/JCGM100/C045315e-html/C045315e_FILES/MAIN_C045315e/04_e.html).

## Items to request from Ben

- Calibration certificate matching the installed sensor and transmitter, including identifiers, date, error results, and the stated uncertainty/coverage factor if provided.
- Full Scale Flow, Low Flow Cutoff percentage, and selected flow units actually used during the October calibration. These are programmable meter settings, not quantities that must come from the calibration certificate. A certificate or setup report may record them, but the historical run settings are needed.
- Sensor size, measuring-tube internal diameter or liner/model information sufficient to establish it, and confirmation of the straight-pipe installation used. The 0.20%-plus-1-mm/s specification assumes an appropriate installation.
- Whether discharge was read from the display, digital output, pulse output, or a 4-20 mA logger. This establishes whether an additional output-chain specification applies.
- Identification of the spreadsheet's +/-0.005-cfs flow-meter entry: its intended meaning and whether it is a separate operational contribution or already includes manufacturer uncertainty. The user accepts reported calibration allowances; the goal is to establish this entry's identity and avoid misinterpretation or double counting, not require it to be rederived.

No Promag information is requested for the new fit. No correspondence has been sent.
