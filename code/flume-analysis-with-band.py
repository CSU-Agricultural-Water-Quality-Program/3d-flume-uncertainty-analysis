#!/usr/bin/env python3
"""
flume-analysis-with-band.py

Stage–discharge plot with per-point error bars, per-flume fits, and shaded bands
based on additive measurement uncertainty propagated through the curve.
"""

import argparse
from pathlib import Path
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

def total_vertical_uncertainty(a, b, h):
    Q = a * h**b
    dQdh = np.abs(a * b * h**(b - 1))
    dQ_flow = FLOWMETER_REL_ERR * np.abs(Q) + INFILTRATION_ERR_CFS
    dQ_depth = dQdh * DEPTH_ERR_FT
    return dQ_flow + dQ_depth

def main():
    parser = argparse.ArgumentParser(description="Flume plot with shaded measurement bands")
    parser.add_argument("csv", help="Path to flume_data.csv")
    parser.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "output" /
                                            "python_measurement_error" / "flume_stage_discharge_with_band.png"),
                        help="Output figure path")
    parser.add_argument("--band", action="store_true", help="Render shaded measurement-uncertainty bands")
    parser.add_argument("--manuf-band", action="store_true", help="Also render a band for the manufacturer curve")
    parser.add_argument("--alpha", type=float, default=0.18, help="Alpha for shaded bands")
    parser.add_argument("--no-errorbars", dest="no_errorbars", action="store_true",
                        help="Suppress per-point error bars while still plotting points")
    args = parser.parse_args()

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

        if not args.no_errorbars:
            ax.errorbar(grp["depth_ft"], grp["flow_cfs"],
                        xerr=grp["depth_err_ft"], yerr=grp["flow_err_cfs"],
                        fmt=marker, markersize=5, mfc=color, mec="black", mew=0.6,
                        ecolor=color, elinewidth=1, capsize=3, alpha=0.9, label=LABELS.get(key, key))
        else:
            ax.plot(grp["depth_ft"], grp["flow_cfs"], linestyle="none",
                    marker=marker, markersize=5, mfc=color, mec="black", mew=0.6,
                    alpha=0.9, label=LABELS.get(key, key))

        a, b, r2 = fit_power(grp["depth_ft"], grp["flow_cfs"])
        fit_params[key] = (a, b, r2)

        hgrid = np.linspace(max(1e-6, grp["depth_ft"].min()*0.9), grp["depth_ft"].max()*1.05, 400)
        Qhat = a * hgrid**b
        ax.plot(hgrid, Qhat, color=color, linewidth=1.8, alpha=0.95)

        if args.band:
            dQ = total_vertical_uncertainty(a, b, hgrid)
            ax.fill_between(hgrid, Qhat - dQ, Qhat + dQ, color=color, alpha=args.alpha, linewidth=0.0)

    hmin = max(0.0, df["depth_ft"].min()*0.8)
    hmax = df["depth_ft"].max()*1.07
    hgrid_full = np.linspace(hmin, hmax, 600)
    Qm = manufacturer_Q(hgrid_full)
    ax.plot(hgrid_full, Qm, linestyle=":", color="black", alpha=0.9,
            label="Manufacturer Eq. (Q = 1.55 H$^{2.58}$)")

    if args.manuf_band and args.band:
        dQm = total_vertical_uncertainty(1.55, 2.58, hgrid_full)
        ax.fill_between(hgrid_full, Qm - dQm, Qm + dQm, color="black", alpha=args.alpha*0.6, linewidth=0.0)

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
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    print(f"[OK] Figure saved to {args.out}")

if __name__ == "__main__":
    main()
