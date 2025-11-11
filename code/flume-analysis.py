#!/usr/bin/env python3
"""
flume-analysis.py

Generates a stage–discharge plot with horizontal (depth) and vertical (flow) error bars
for flume calibration data, fits power-law curves for each flume type, and overlays
the manufacturer equation for a 60° V trapezoidal flume.

Usage:
    python flume-analysis.py "../data/flume_data.csv" --out "../output/flume_stage_discharge_with_errors.png"

CSV must have: depth_ft, flow_cfs, flume_type
Example flume_type values: 3D1, 3D2, FIBER, EQ (EQ is optional)
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------- Constants ---------------------------------
DEPTH_ERR_FT = 0.005          # ±0.005 ft, point gauge
FLOWMETER_REL_ERR = 0.0002    # ±0.02 % of reading
INFILTRATION_ERR_CFS = 0.02   # ±0.02 cfs

LABELS = {
    "3D1": "3D Printed 1",
    "3D2": "3D Printed 2",
    "FIBER": "Fiberglass",
    "EQ": "Manufacturer Eq."
}

LEGEND_ORDER = ["FIBER", "3D1", "3D2", "EQ"]  # optional, keeps legend stable

# ------------------------- Helpers -----------------------------------
def manufacturer_Q(h_ft: np.ndarray) -> np.ndarray:
    """Manufacturer equation: Q = 1.55 * H^2.58 (cfs)."""
    return 1.55 * (h_ft ** 2.58)

def fit_power(x, y):
    """Fit y = a * x^b via log–log least squares; return (a, b, R2)."""
    x = np.asarray(x); y = np.asarray(y)
    m = (x > 0) & (y > 0)
    x, y = x[m], y[m]
    lx, ly = np.log(x), np.log(y)
    b, loga = np.polyfit(lx, ly, 1)
    a = np.exp(loga)
    yhat = a * x**b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return a, b, r2

# ------------------------- Main --------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Flume calibration uncertainty analysis plotter")
    ap.add_argument("csv", help="Path to flume_data.csv")
    ap.add_argument("--out", default="../output/flume_stage_discharge_with_errors.png", help="Output figure path")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    # Error bars
    df["depth_err_ft"] = DEPTH_ERR_FT
    df["flow_err_cfs"] = FLOWMETER_REL_ERR * df["flow_cfs"].abs() + INFILTRATION_ERR_CFS

    # Plot scaffolding
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fit_params = {}

    # Markers and colors
    markers = {"3D1": "o", "3D2": "s", "FIBER": "D"}
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    flumes = [f for f in df["flume_type"].unique() if f != "EQ"]
    color_map = {f: color_cycle[i % len(color_cycle)] for i, f in enumerate(flumes)}

    # Per-flume points, error bars, and fits
    for key, grp in df.groupby("flume_type"):
        if key == "EQ":
            continue

        color = color_map.get(key, "gray")
        marker = markers.get(key, "o")

        ax.errorbar(
            grp["depth_ft"], grp["flow_cfs"],
            xerr=grp["depth_err_ft"], yerr=grp["flow_err_cfs"],
            fmt=marker, markersize=5, mfc=color, mec="black", mew=0.6,
            ecolor=color, elinewidth=1, capsize=3, alpha=0.9, label=LABELS.get(key, key)
        )

        a, b, r2 = fit_power(grp["depth_ft"], grp["flow_cfs"])
        fit_params[key] = (a, b, r2)

        hgrid = np.linspace(max(1e-6, grp["depth_ft"].min() * 0.9), grp["depth_ft"].max() * 1.05, 200)
        ax.plot(hgrid, a * hgrid ** b, color=color, linewidth=1.8, alpha=0.9)

    # Manufacturer curve
    hmin = max(0.0, df["depth_ft"].min() * 0.8)
    hmax = df["depth_ft"].max() * 1.07
    hgrid_full = np.linspace(hmin, hmax, 300)
    ax.plot(hgrid_full, manufacturer_Q(hgrid_full), linestyle="--", label="Manufacturer Eq. (Q = 1.55 H$^{2.58}$)")

    # Labels, grid
    ax.set_xlabel("Stage (ft)")
    ax.set_ylabel("Flow rate (cfs)")
    ax.grid(True, alpha=0.3)

    # Stable legend order
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(LABELS[k]) for k in LEGEND_ORDER if LABELS.get(k) in labels]
    ax.legend([handles[i] for i in order], [labels[i] for i in order])

    # Fit annotations
    y_frac = 0.92
    for k in [kk for kk in LEGEND_ORDER if kk in fit_params]:
        a, b, r2 = fit_params[k]
        txt = f"{LABELS.get(k,k)}: Q = {a:.2f} H^{b:.2f}, R² = {r2:.2f}"
        ax.text(0.02, y_frac, txt, transform=ax.transAxes, fontsize=8)
        y_frac -= 0.06

    fig.tight_layout()
    fig.savefig(args.out, bbox_inches="tight")
    print(f"[OK] Figure saved to {args.out}")

if __name__ == "__main__":
    main()
