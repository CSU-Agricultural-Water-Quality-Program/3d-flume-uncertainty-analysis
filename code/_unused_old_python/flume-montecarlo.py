#!/usr/bin/env python3
"""
flume-montecarlo.py

Monte Carlo refitting of Q = a * H^b using known measurement uncertainties.
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

def fit_power(x, y):
    m = (x > 0) & (y > 0)
    lx, ly = np.log(x[m]), np.log(y[m])
    b, loga = np.polyfit(lx, ly, 1)
    return np.exp(loga), b

def manufacturer_Q(h_ft):
    return 1.55 * (h_ft ** 2.58)

def montecarlo_refit(df, n_iter=1000, random_state=None):
    rng = np.random.default_rng(random_state)
    results = {}
    for ftype, grp in df.groupby("flume_type"):
        if ftype == "EQ":
            continue
        H_obs = grp["depth_ft"].to_numpy()
        Q_obs = grp["flow_cfs"].to_numpy()

        a_draws = np.full(n_iter, np.nan)
        b_draws = np.full(n_iter, np.nan)

        for i in range(n_iter):
            H_sim = rng.normal(H_obs, DEPTH_ERR_FT)
            Q_sigma = FLOWMETER_REL_ERR * Q_obs + INFILTRATION_ERR_CFS
            Q_sim = rng.normal(Q_obs, Q_sigma)
            mask = (H_sim > 0) & (Q_sim > 0)
            if np.sum(mask) < 3:
                continue
            a_i, b_i = fit_power(H_sim[mask], Q_sim[mask])
            a_draws[i], b_draws[i] = a_i, b_i

        h_grid = np.linspace(H_obs.min()*0.9, H_obs.max()*1.05, 220)
        Qmat = np.array([a*h_grid**b for a, b in zip(a_draws, b_draws)])
        Q_mean = np.nanmean(Qmat, axis=0)
        Q_lo, Q_hi = np.nanpercentile(Qmat, [5, 95], axis=0)

        results[ftype] = dict(a=a_draws, b=b_draws, grid=h_grid,
                              mean=Q_mean, lo=Q_lo, hi=Q_hi)
    return results

def plot_results(df, results, out_path, alpha_band=0.2):
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    for ftype, res in results.items():
        color = COLOR_MAP.get(ftype, "gray")
        marker = MARKERS.get(ftype, "o")
        sub = df[df["flume_type"] == ftype]
        ax.scatter(sub["depth_ft"], sub["flow_cfs"], label=LABELS.get(ftype, ftype),
                   alpha=0.8, color=color, edgecolor="k", s=30, marker=marker)
        ax.fill_between(res["grid"], res["lo"], res["hi"], color=color, alpha=alpha_band)
        ax.plot(res["grid"], res["mean"], color=color, linewidth=1.8)

    h_full = np.linspace(df["depth_ft"].min()*0.8, df["depth_ft"].max()*1.1, 300)
    ax.plot(h_full, manufacturer_Q(h_full), linestyle=":", color="black",
            alpha=0.9, label="Manufacturer Eq. (Q = 1.55 H$^{2.58}$)")

    ax.set_xlabel("Stage (ft)")
    ax.set_ylabel("Flow rate (cfs)")
    ax.grid(True, alpha=0.3)

    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(LABELS[k]) for k in LEGEND_ORDER if LABELS.get(k) in labels]
    # Include manufacturer last if present
    if "Manufacturer Eq. (Q = 1.55 H$^{2.58}$)" in labels:
        order += [labels.index("Manufacturer Eq. (Q = 1.55 H$^{2.58}$)")]
    ax.legend([handles[i] for i in order], [labels[i] for i in order])

    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    print(f"[OK] Saved Monte Carlo ensemble plot to {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Monte Carlo refitting of flume rating curves")
    parser.add_argument("csv", help="Path to flume_data.csv")
    parser.add_argument("--out", default=str(Path(__file__).resolve().parents[2] / "output" /
                                            "python_monte_carlo" / "flume_montecarlo_band.png"),
                        help="Output PNG path")
    parser.add_argument("--n", type=int, default=1000, help="Monte Carlo iterations")
    parser.add_argument("--alpha", type=float, default=0.2, help="Transparency for shaded bands")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    results = montecarlo_refit(df, n_iter=args.n)
    plot_results(df, results, args.out, alpha_band=args.alpha)

if __name__ == "__main__":
    main()
