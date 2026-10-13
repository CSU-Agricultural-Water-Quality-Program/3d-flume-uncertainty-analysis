# Flow meter specifications and Bayesian calibration review

Reviewed 13 October 2026. The evidence supports an advantage for the smaller reference meter at low discharge. It does not support representing its uncertainty as a constant 0.2% coefficient of variation, or treating its listed minimum flow as the programmed shutoff. The observation model needs revision before selecting final priors. No calibration code, fitted model, workbook, or manuscript was edited or rerun.

## Sources and verification

The attached Badger Meter manual is MAG-UM-01272-EN-43, August 2026. Pages 66 and 76 were extracted and visually inspected. Pages 6, 12, 33, and 39 were read for the measurement principle, installation, and cutoff settings. The copy under REFERENCES has the same SHA-256 hash as the attachment. The Promag technical information in REFERENCES was read at pages 9, 11, 34-36, and 91; pages 9 and 34 were visually inspected. Page numbering matches the PDF page numbers for these locations.

- [ModMAG manual](<D:/OneDrive - Colostate/Cabot,Perry's files - 3D Printing Flumes/ModMag Manual.pdf>).
- [Promag technical information](<D:/OneDrive - Colostate/Cabot,Perry's files - 3D Printing Flumes/REFERENCES/Promag Manual (1).pdf>).
- [October workbook](<C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/data/New Data Oct 2026/Oct 2026 Trapezoidal Flume Calibration Final.xlsx>), Detailed Calibration (10-point), B2:C5 and calibration rows 10:19.
- [Active code](<C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis/code/bayesflume.r>), especially lines 194-207 and 282-312.
- [Manuscript V2](<D:/OneDrive - Colostate/Cabot,Perry's files - 3D Printing Flumes/2_Second Submission/Deleon et al_Additive_Manufacturing_for_Irrigation_MD_V2.docx>). The file was locked. Its open Word document was read through COM without modifying it. The extracted current text excludes tracked deletions and includes tracked insertions. This is the open document state, which can include unsaved revisions. Equations were also inspected through their Office Math text. This was a content audit, not a page-layout audit.

Manufacturer cross-check: [Badger Meter product specifications](https://www.badgermeter.com/products/meters/electromagnetic/m2000-electromagnetic-flow-meter/) explicitly give 0.20% of rate plus 1 mm/s with appropriate straight pipe runs. Metrological conversion guidance: [JCGM 100:2008, section 4.3](https://www.iso.org/sites/JCGM/GUM/JCGM100/C045315e-html/C045315e_FILES/MAIN_C045315e/04_e.html). Modeling and diagnostics: [Stan measurement error guidance](https://mc-stan.org/docs/stan-users-guide/measurement-error.html), [Stan efficiency tuning](https://mc-stan.org/docs/stan-users-guide/efficiency-tuning.html), and [Stan convergence guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).

## Verified specifications and interpretation

| Quantity | ModMAG M2000 | Promag W400 |
| --- | --- | --- |
| Standard accuracy/error specification | ±0.20% of measured value plus ±1 mm/s, p. 66 | Maximum measured error ±0.5% of reading plus ±1 mm/s under reference conditions, p. 34 |
| Separate repeatability specification | ±0.1%, p. 66 | Maximum ±0.1% of reading plus ±0.5 mm/s, p. 36 |
| Ben's cited low-flow value | 0.85 m³/h is the listed lower flow-range endpoint for 4-in/DN100, p. 76 | 10 m³/h is the listed factory cutoff for 12-in/DN300, p. 9 |
| Actual cutoff | Programmable 0-10% of Full Scale Flow; displayed default 0.2%, p. 39 | Factory value does not establish the value used during the experiment |

The Promag also has an optional 0.2% of reading plus 2 mm/s specification and a Flat Spec option (pp. 34-35). The installed option must be established before applying either. Its US-unit factory table is rounded separately; retain the cited SI table value when reproducing Ben's conversion. The ModMAG's separate zero-straight-run accuracy is not interchangeable with its normal accuracy specification. Installation evidence is therefore necessary.

The email's unit conversions are correct to its reported precision:

- 10 m³/h = 2.777778 L/s = 0.0980963 cfs.
- 0.85 m³/h = 0.236111 L/s = 0.00833819 cfs.

Use m³/h rather than the ambiguous abbreviation cmh. A cutoff suppresses the reported reading below a configured threshold. It is not a standard deviation, an additive uncertainty, or necessarily a fundamental physical detection limit. In particular, 0.85 m³/h is not established as the actual ModMAG cutoff.

The manual supports Faraday's induction principle and requires a full pipe (p. 6). Continuity gives velocity v=Q/A. At nominal diameters 100 and 300 mm, velocity is nine times greater in the smaller tube at equal discharge. The manual explicitly notes that a smaller meter installed with reducers may improve low-flow accuracy (p. 12). This supports the reasoning in the email. It does not establish the installed uncertainty or imply that induced voltage rises by exactly the same factor: the electromagnetic measurement relation also contains electrode spacing and field strength.

## Expressing meter error in L/s

Let D be the actual measuring-tube internal diameter in metres, A=πD²/4, Q be discharge in L/s, and r the proportional specification as a fraction. A velocity error of 1 mm/s corresponds to A L/s because 0.001 m/s × 1000 L/m³=1. Thus the specification half-width is

E(Q)=r|Q|+A, in L/s.

The terms add in the manufacturer's error envelope. Treating them as independent standard deviations and combining them in quadrature would be a different assumption.

For provisional nominal DN100 and DN300 calculations:

- ModMAG: E(Q)=0.002|Q|+0.00785398 L/s.
- Standard Promag: E(Q)=0.005|Q|+0.0706858 L/s.

Nominal diameter is an approximation. For a literal 4-in diameter, the ModMAG additive term is 0.00810732 L/s, and for 102 mm it is 0.00817128 L/s. Promag p. 91 lists internal diameters that depend on liner and design. Use the installed sensor dimensions/certificate for final values.

| New calibration flow (cfs) | Flow (L/s) | Nominal ModMAG error half-width (L/s) | Equivalent percentage of reading |
| --- | --- | --- | --- |
| 0.02 | 0.566337 | 0.00898666 | 1.58680% |
| 0.10 | 2.831685 | 0.01351735 | 0.477361% |
| 0.32 | 9.061391 | 0.02597676 | 0.286675% |

These values are calculated specification envelopes, not observed errors or posterior intervals. They demonstrate why 0.2% alone understates the specified error, especially at the lowest flow.

## Operational tolerance and unresolved probability interpretation

The October workbook reports ±0.005 cfs for the ten-point calibration and ±0.01 cfs for the base calibration. These equal ±0.141584 and ±0.283168 L/s. The ten-point operational half-width corresponds to 25% of the lowest 0.02-cfs reading and 1.5625% of the highest 0.32-cfs reading. It exceeds the nominal manufacturer's envelope throughout the new dataset.

The workbook labels these as accuracy but does not identify a distribution, coverage factor, repeatability calculation, or whether they include the meter specification. Therefore the final operational-error SD cannot be selected from the workbook alone. In particular, do not treat ±0.005 cfs as 0.5%, and do not add it to the manufacturer uncertainty if it already represents total installed uncertainty.

For a provisional Type B assessment, if ±E is interpreted as symmetric bounds with equal probability within them, standard uncertainty is u=E/√3 (JCGM 100:2008, 4.3.7). A Gaussian likelihood using that SD is a moment approximation to this bounded uncertainty, not a statement that the manufacturer's error specification is Gaussian. If a certificate instead supplies expanded uncertainty with coverage factor k, use u=U/k. Divide by 1.96 only if the stated interval really has 95% Gaussian coverage.

Under the rectangular assumption:

- Manufacturer-only ModMAG u(Q)=(0.002|Q|+0.00785398)/√3 L/s.
- Operational ±0.005 cfs gives u_oper=0.0817437 L/s.
- If operational tolerance is a separate independent contribution, combine these standard uncertainties in quadrature.
- If it is total installed uncertainty, use that total instead of counting its included components again.

A calibration uncertainty can contain systematic components common to all observations from the same meter. Repeated independent Normal noise does not reproduce a shared calibration bias. Where supported by the certificate, represent a shared gain/zero correction separately from repeatability. Do not infer independent contributions from an undivided total-error envelope.

## Data checks relevant to the prior decision

The active CSV contains 30 observations, ten for each flume, at 0.02-0.32 cfs (0.566337-9.061391 L/s). All exceed the ModMAG's listed minimum of 0.236111 L/s. With nominal DN100, velocities span 0.0721083-1.153732 m/s. The actual programmed cutoff remains unknown.

Four observations in the original CSV fall below the Promag's published factory cutoff: fiberglass at 0.04 and 0.08 cfs, and 3D Printed 2 at 0.03 and 0.08 cfs. This does not prove they are invalid: the cutoff may have been changed, or their flow provenance may differ. Ben should explain how they were obtained before they are reused or compared as equivalent reference measurements. Do not enlarge a Normal-error prior to compensate for an actual suppressed/censored output.

The new CSV rounds stage to 0.001 ft. Comparison against the source workbook finds a maximum stage-rounding difference of 0.0132 cm, exceeding the workbook's stated ±0.01-cm gauge tolerance. Use the original water-surface-minus-bottom readings when preparing metric inputs. All 30 flow values match the corresponding workbook values, and all stage values agree within the existing CSV rounding tolerance. These checks validate the transfer, not the experimental measurement accuracy.

## Required code and prior changes

1. **Replace the current meter interpretation.** `cv_Q ~ LogNormal(log(0.02),0.5)` has median 2% as an SD fraction. It is not the M2000's full specified error. A single change to 0.002 still omits the velocity-dependent absolute term and equates an accuracy envelope with an SD. Under the provisional rectangular interpretation, the proportional SD coefficient is 0.002/√3=0.00115470 and the nominal additive SD contribution is 0.00453450 L/s; these are combined as the single envelope divided by √3, not as independent variance components.

2. **Use instrument-informed uncertainty as data before freely estimating it.** For a manufacturer-only sensitivity case, use the complete u_meter(Q) above in the measurement likelihood. An operational-total sensitivity case can use 0.0817437 L/s under its explicitly stated rectangular interpretation. These cases can be compared while waiting for Ben; neither should be described as the verified installed SD. If a free multiplier is retained, its prior represents uncertainty about the applicability of the specified envelope, and its width needs independent justification. The manual does not supply a lognormal prior on a CV.

3. **Resolve the additive term's meaning.** The current likelihood adds `sigma_Q_add[F]` inside a zero-mean error SD. This is symmetric measurement scatter. It does not remove discharge and cannot estimate a directional infiltration loss or an additive mean offset. Its prior median of 0.02 cfs is 0.566337 L/s, equal to the lowest calibration flow. Retain that prior only if independent facility evidence supports that quantity as an SD. If leakage is a directional loss between meter and flume, write a mass-balance mean relation with a nonnegative loss parameter, with sign determined by meter location. If it is additional nonspecific scatter, name it accordingly.

4. **Reduce unsupported variance decomposition.** Process SD and additional constant measurement SD are difficult to distinguish with one observation per stage. When proportional error tends to zero, marginal variance depends on their squared sum, giving exact confounding; a small proportional term does not reliably resolve it. Fix well-supported instrument uncertainty and retain a single residual component unless independent measurements support a separate loss/error component. The meter manual cannot determine the process or leakage priors.

5. **Resolve and propagate stage uncertainty.** The workbook/code use ±0.01 cm; the current manuscript says ±1.5 mm. The former is 0.1 mm, below the manuscript's stated 1-mm ruler graduation, and needs clarification as accuracy versus resolution. `sigma_H` is stored and plotted but absent from the active likelihood. A latent-stage model can represent gauge error explicitly. Account for both water-surface and bottom readings, and for any shared datum uncertainty; do not assume their uncertainties are independent without evidence. Absorbing unspecified stage error into a constant process SD does not explicitly propagate it through the nonlinear curve.

6. **Improve parameterization before rerunning.** The saved October diagnostics show 744 divergences, 562 maximum-treedepth hits, and E-BFMI below 0.3 in all four chains. Review the latent-discharge parameterization and variance confounding. A noncentered latent parameterization or appropriate marginalization may help, but must preserve the intended likelihood. Metric units alone cannot establish convergence. Check divergences, R-hat, effective sample sizes, energy diagnostics, prior predictions, and posterior predictions before updating scientific conclusions.

7. **Write the corrected run into a new output directory.** Preserve the October fit and its exports for provenance. Record the meter model, diameter, specification components, uncertainty conversion, operational-error interpretation, input hash, priors, and cutoff provenance. Update the plotting/component scripts when likelihood terms change; their present labels and variance calculations assume the old model.

## Metric priors and an equivalent change of units

Using cm and L/s throughout is appropriate. Conversion of the data alone is insufficient because dimensional priors and the nonlinear coefficient must also change.

For exactly equivalent priors to the current independent a_ft and b specification:

log(a_cm)|b ~ Normal(log(1.55)+log(28.316846592)-b log(30.48), 0.5),

b ~ LogNormal(log(2.58),0.5).

The manufacturer coefficient in cm and L/s is 0.00651055 at b=2.58. Using an independent LogNormal(log(0.00651055),0.5) for a_cm does not preserve the original joint prior when b is uncertain. The conditional expression above does preserve it. Equivalently, use a dimensionless positive deviation k_j with median 1 and log-SD 0.5, and set

μ_j(H)=k_j × 1.55 × 28.316846592 × (H/30.48)^{b_j}, with H in cm.

For better interpretability, a reference-stage parameterization μ_j(H)=Q_ref,j(H/10 cm)^{b_j} places the scale prior on flow at 10 cm; the manufacturer value there is 2.475244 L/s. Making log(Q_ref,j) and b_j independent is a deliberate new prior. Assess it with prior predictions rather than describing it as a unit conversion. Sharing prior settings across flumes is not, by itself, population-level partial pooling.

Dimensional prior equivalents:

- `sigma_proc_cfs ~ Exponential(rate=1)` becomes `sigma_proc_lps ~ Exponential(rate=1/28.316846592)`. Its mean is 28.31685 L/s, more than three times the largest calibration flow; that is an equivalence calculation, not a recommendation to retain this broad prior.
- `sigma_Q_add_cfs ~ LogNormal(log(0.02),0.5)` becomes `sigma_Q_add_lps ~ LogNormal(log(0.56633693184),0.5)`, if retained with the same meaning.
- Proportional fractions and b remain dimensionless. A meter error coefficient is not multiplied by a unit conversion.
- Gauge half-widths of 0.01 cm and 1.5 mm equal 0.01 and 0.15 cm respectively. Decide which applies before assigning a probability distribution.

Update input validation, run metadata, manufacturer curves, posterior transformations, prediction grids, figure captions, and table units together. Remove downstream ft/cfs conversions from metric-native draws to prevent double conversion. Preserve physical meaning when comparing old and new priors.

## Manuscript correction register

| Location | Current inconsistency | Needed revision |
| --- | --- | --- |
| Hydraulic Calibration; Table 2 | ModMag 109-1808 and ±1.5% of reading; discussion later compares against ±2% | Confirm sensor/transmitter identifiers; report M2000's full error specification and the separate operational tolerance with its meaning |
| Stage measurement description; Table 2 | ±1.5 mm differs from workbook ±0.01 cm | Resolve with Ben and match the actual model input |
| Eq. 3 and error-component descriptions | Additive SD described as discharge loss | Describe symmetric scatter or implement an explicit mass-balance loss model |
| Eqs. 6-10 | Imperial scale prior; additive prior says 0.2 while active code uses 0.02; CV prior says 2×10^-4 while active code uses 0.02 | Replace with the final implemented metric model and its probability assumptions |
| Process interpretation | Claims point-gauge uncertainty is included | Implement it explicitly or state that stage is treated as fixed; avoid claiming complete propagation |
| Calibration sample description; Table 3 | Six points per printed flume, four for fiberglass | Report ten observations per flume and the actual measured metric range |
| Table 4 caption and subsequent discussion | Manufacturer exponent 2.55 differs from 2.58 in methods/code | Make the reference equation consistent with the verified flume source |
| Meter-error result | Posterior ±0.022% interpreted as evidence of superior meter precision | Replace after a validated refit; estimated residual CV is not a direct calibration test of the meter, particularly under strong priors and confounded error components |
| Convergence paragraph and supplement | Says all analyses converged | Replace with diagnostics for the corrected fit; current October diagnostics do not support this statement |
| Results, Tables 4-5, Figure 6, nomenclature | Old-fit summaries, imperial units, and additive-offset interpretation | Regenerate after model validation; distinguish mean-curve credible intervals from predictive intervals for latent or observed discharge |
| Fiberglass comparison | Agreement used as accuracy justification | Present as corroborating consistency with a reference curve, not a traceable independent meter calibration; account for the manufacturer's curve being used to inform priors |

These are content corrections. No tracked changes were accepted or modified.

## Details to request from Ben

1. Sensor and transmitter nameplate photographs or exact identifiers, installed internal diameter/liner, and the calibration certificate, including uncertainty definition and coverage factor if supplied.
2. ModMAG Full Scale Flow, Low Flow Cutoff percentage, flow units, and settings actually used during the ten-point runs. These establish the programmed zero threshold: cutoff fraction × Full Scale Flow.
3. The basis of ±0.005 cfs: instrument resolution, fluctuations, repeated-measurement SD, maximum deviation, or total installed accuracy; whether it already includes the manufacturer's specification; and whether listed flow values are actual averaged readings or nominal targets.
4. The applicable point-gauge accuracy, its resolution, treatment of the bottom datum, and whether ±0.01 cm is intended or a transcription error.
5. Promag cutoff/settings and provenance of old nonzero readings below its factory threshold, plus confirmation of installed accuracy option.
6. Independent evidence for the stated infiltration quantity and whether it applies to the recalibration setup. Confirm its interpretation as a loss, a bound, or a standard deviation and the meter's location relative to the potential loss.

No correspondence was sent. The user stated that the certificate/settings are not currently available and will request specifics from Ben.

## Audit outputs

`meter_error_by_calibration_point.csv` contains all 30 points, metric conversions, source-stage comparisons, nominal-diameter specification envelopes, and explicitly assumed rectangular standard uncertainties. `audit_calculations.json` records counts, conversions, and unresolved assumptions. `source_manifest.json` records source hashes, including the read-only snapshot of the open Word document. These outputs quantify the audit; they are not replacement model results.
