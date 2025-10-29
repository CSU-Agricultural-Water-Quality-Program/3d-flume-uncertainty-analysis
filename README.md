# 3D-Printed Trapezoidal Flume Calibration and Sensitivity Analysis  
**Author:** AJ Brown and Emmanuel Deleon [Agricultural Water Quality Program, CSU]  
**Collaborators:** Daniel White (Hydraulics Lab), Perry Cabot, Erik Wardle 
**Project Year:** 2025-2026  

---

## Overview  

This repository supports analysis for the CSU Agricultural Water Quality Program’s **novel 3D-printed trapezoidal flume**, developed to provide a low-cost, field-deployable alternative to commercial fiberglass units.  

Building on hydraulic testing performed at the **CSU Hydraulics Laboratory**, this project uses Python to:  
1. Quantify **measurement uncertainty** from flow-rate and depth sensors;  
2. Perform **sensitivity analysis** to evaluate the effects of each measurement error on flume calibration;  
3. Statistically compare measured discharge against the **manufacturer’s rating equation** for a 60° V-notch trapezoidal flume.  

The ultimate goal is to validate whether 3D-printed PETG flumes can reliably reproduce manufacturer discharge relationships while accounting for real-world error sources such as **flume infiltration** and **gauge precision**.

---

## Background  

### Research Context  
Additive manufacturing allows rapid, low-cost fabrication of large-scale hydraulic structures.  
The CSU team printed full-size trapezoidal and Parshall flumes using **PETG filament** on an **Elegoo OrangeStorm Giga (800×800×1000 mm)** printer, modeled in *Autodesk Fusion 360* following Teledyne ISCO’s master dimensions.  

Hydraulic calibration at CSU’s Hydraulics Lab compared two 3D-printed flumes with a commercial fiberglass unit using a magnetic flow meter (±0.02 %) and point gauge (±0.005 ft). All units exhibited nearly identical power-law stage-discharge relationships:  

\[
Q = a h^b
\]
with \(R^2 \approx 1.0\) across all flume types.  

### Manufacturer’s Reference Equation  
For the large 60° trapezoidal flume:  
\[
Q = 1.55 \, h^{2.58} \, \text{(cfs)}
\]
where \(h\) = head in feet.  

---

## Research Objectives  

| Objective | Description |
|------------|-------------|
| **1. Quantify Measurement Uncertainty** | Integrate flow-meter, point-gauge, and infiltration error sources into propagation models. |
| **2. Sensitivity Analysis** | Determine which uncertainties (depth vs flow rate vs infiltration) most affect calibration accuracy. |
| **3. Equation Comparison** | Test whether observed discharge differs significantly from the manufacturer equation using confidence intervals and hypothesis testing. |
| **4. Confidence Bands** | Build 95 % prediction intervals for the calibration curve \(Q = a h^b\). |
| **5. Communication** | Produce publication-ready figures and scripts for replication in the 3D flume manuscript (target submission: *Transactions of the ASABE*, Jan 2026). |

---

## Repository Structure  

```
flume_sensitivity_analysis/
│
├── data/
│   ├── flume_data.csv
│   ├── large_60-degree-v_trapezoidal_flume_discharge_table.pdf
│   └── Trapezoidal Flume Testing.pdf
│
├── src/
│   ├── sensitivity_analysis.py
│   ├── compare_to_manufacturer.py
│   └── plot_calibration_curves.py
│
├── notebooks/
│   └── 01_uncertainty_analysis.ipynb
│
├── results/
│   ├── calibration_fit_summary.csv
│   ├── uncertainty_contributions.png
│   └── stage_discharge_comparison.png
│
└── README.md
```

---

## Analytical Approach  

1. **Import and Clean Data**
   ```python
   import pandas as pd
   df = pd.read_csv("data/flume_data.csv")
   ```
2. **Fit Power-Law Rating Curves**
   \[
   Q = a h^b
   \]
   via nonlinear least squares (`scipy.optimize.curve_fit`).
3. **Propagate Measurement Uncertainty**
   - Flow meter: ±0.02 %  
   - Point gauge: ±0.005 ft  
   - Infiltration: ±0.02 cfs  
   - Additive error model to estimate total uncertainty in \(Q\).
4. **Sensitivity Analysis**
   - Use known measurement uncertainties to create a stacked confidence interval.
5. **Statistical Comparison**
   - Compare empirical \(a, b\) to manufacturer constants (1.55, 2.58).  
6. **Visualization**
   - Overlay experimental data, fitted curve, and manufacturer equation.  
   - Shaded uncertainty bands showing propagated measurement error.


---

## Expected Outputs  

| Output | Description |
|--------|-------------|
| `uncertainty_contributions.png` | Tornado plot showing which variables dominate uncertainty. |
| `stage_discharge_comparison.png` | Overlay of 3D-printed vs fiberglass vs manufacturer equation. |
| `calibration_fit_summary.csv` | Table of fitted coefficients \(a, b\), standard errors, and \(R^2\). |

---

## License
This project is licensed under the GNU General Public License v2.0. See the [LICENSE](LICENSE) file for details.
