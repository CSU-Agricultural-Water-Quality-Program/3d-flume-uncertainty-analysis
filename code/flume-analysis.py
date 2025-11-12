#!/usr/bin/env python3
"""
flume-analysis.py

Stage–discharge plot with per-point error bars and per-flume power-law fits.
Formatting is unified with the other project plots:
  - Manufacturer equation is black, dotted, and appears in the legend.
  - Point markers are shaped by flume type.
  - Colors are consistent across scripts (blue=3D1, orange=3D2, green=Fiberglass).
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEPTH_ERR_FT = 0.005
FLOWMETER_REL_ERR = 0.0002
INFILTRATION_ERR_CFS = 0.02

LABELS = {"3D1": "3D Printed 1", "3D2": "3D Printed 2", "FIBER": "Fiberglass", "EQ": "Manufacturer Eq."}
LEGEND_ORDER = ["FIBER", "3D1", "3D2", "EQ"]

COLOR_MAP = {"3D1": plt.rcParams['axes.prop_cycle'].by_key()['color'][0],
             "3D2": plt.rcParams['axes.prop_cycle'].by_key()['color'][1],
             "FIBER": plt.rcParams['axes.prop_cycle'].by_key()['color'][2]}

MARKERS = {"3D1": "o", "3D2": "s", "FIBER": "D"}

def manufacturer_Q(h_ft: np.ndarray) -> np.ndarray:
    return 1.55 * (h_ft ** 2.58)

def fit_power(x, y):
    x = np.asarray(x); y = np.asarray(y)
    m = (x > 0) & (y > 0)
    x, y = x[m], y[m]
    lx, ly = np.log(x), np.log(y)
    b, loga = np.polyfit(lx, ly, 1)
    a = np.exp(loga)
    yhat = a * x**b
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return a, b, r2

def main():
    ap = argparse.ArgumentParser(description="Flume plot with error bars and power-law fits")
    ap.add_argument("csv", help="Path to flume_data.csv")
    ap.add_argument("--out", default="../output/flume_stage_discharge_with_errors.png", help="Output figure path")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    df["depth_err_ft"] = DEPTH_ERR_FT
    df["flow_err_cfs"] = FLOWMETER_REL_ERR * df["flow_cfs"].abs() + INFILTRATION_ERR_CFS

    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fit_params = {}

    for key, grp in df.groupby("flume_type"):
        if key == "EQ":
            continue
        color = COLOR_MAP.get(key, "gray")
        marker = MARKERS.get(key, "o")

        ax.errorbar(
            grp["depth_ft"], grp["flow_cfs"],
            xerr=grp["depth_err_ft"], yerr=grp["flow_err_cfs"],
            fmt=marker, markersize=5, mfc=color, mec="black", mew=0.6,
            ecolor=color, elinewidth=1, capsize=3, alpha=0.9, label=LABELS.get(key, key)
        )

        a, b, r2 = fit_power(grp["depth_ft"], grp["flow_cfs"])
        fit_params[key] = (a, b, r2)
        hgrid = np.linspace(max(1e-6, grp["depth_ft"].min()*0.9), grp["depth_ft"].max()*1.05, 200)
        ax.plot(hgrid, a*hgrid**b, color=color, linewidth=1.8, alpha=0.95)

    hmin = max(0.0, df["depth_ft"].min()*0.8)
    hmax = df["depth_ft"].max()*1.07
    hgrid_full = np.linspace(hmin, hmax, 300)
    ax.plot(hgrid_full, manufacturer_Q(hgrid_full), linestyle=":", color="black",
            alpha=0.9, label="Manufacturer Eq. (Q = 1.55 H$^{2.58}$)")

    ax.set_xlabel("Stage (ft)")
    ax.set_ylabel("Flow rate (cfs)")
    ax.grid(True, alpha=0.3)

    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(LABELS[k]) for k in LEGEND_ORDER if LABELS.get(k) in labels]
    ax.legend([handles[i] for i in order], [labels[i] for i in order])

    y_frac = 0.92
    for k in [kk for kk in LEGEND_ORDER if kk in fit_params]:
        a, b, r2 = fit_params[k]
        ax.text(0.02, y_frac, f"{LABELS.get(k,k)}: Q = {a:.2f} H^{b:.2f}, R² = {r2:.2f}",
                transform=ax.transAxes, fontsize=8)
        y_frac -= 0.06

    fig.tight_layout()
    fig.savefig(args.out, bbox_inches="tight")
    print(f"[OK] Figure saved to {args.out}")

if __name__ == "__main__":
    main()
