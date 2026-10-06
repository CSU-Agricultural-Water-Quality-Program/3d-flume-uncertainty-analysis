#!/usr/bin/env python3
"""
simulate_flume_data_model3.py

Simulate flume calibration data consistent with "Model 3" in bayesflume.r.

Model 3 (from bayesflume.r) assumptions:
  mu = a[F] * H_obs^b[F]
  Q_true ~ Normal(mu, sigma_proc[F])
  Q_obs  ~ Normal(Q_true, sqrt((cv_Q * Q_true)^2 + sigma_Q_add[F]^2))

This script:
  - reads your real data (data/flume_data.csv) to mimic:
      * which flume types exist (excluding EQ)
      * how many points per flume
      * depth ranges per flume
  - generates simulated data from the Model 3 generative process
  - writes:
      data/sim_data.csv
      output/simulation/sim_truth_model3.csv  (true parameters used)

Run from repo root:
  python code/simulate_flume_data_model3.py
"""

from __future__ import annotations

import math
from pathlib import Path
import numpy as np
import pandas as pd


# ----------------------------
# User-tunable simulation knobs
# ----------------------------

SEED = 20260211

# If True, depths are sampled within each flume's observed depth range in your real data.
# If False, depths will be generated on an evenly spaced grid within each flume's range.
DEPTHS_RANDOM_WITHIN_RANGE = True

# Small jitter added to depths (ft) to avoid identical values in case you use grids.
# Set to 0.0 if you want perfectly clean H values.
DEPTH_JITTER_FT = 0.002

# Enforce physically meaningful positivity:
MIN_POSITIVE_Q = 1e-6
MIN_POSITIVE_H = 1e-6

# If a simulated draw produces negative Q_true or Q_obs, we resample up to this many times,
# then fall back to clipping at MIN_POSITIVE_Q.
RESAMPLE_TRIES = 100


def _repo_root() -> Path:
    """Return the repository root assuming this file lives in <root>/code/."""
    return Path(__file__).resolve().parents[1]


def _positive_resample_normal(rng: np.random.Generator, mean: float, sd: float) -> float:
    """
    Draw from Normal(mean, sd) but enforce positivity (by resampling).
    If it fails repeatedly (rare unless sd is huge), clip to MIN_POSITIVE_Q.
    """
    if sd <= 0:
        return max(mean, MIN_POSITIVE_Q)

    for _ in range(RESAMPLE_TRIES):
        x = rng.normal(mean, sd)
        if x > MIN_POSITIVE_Q:
            return x

    # Fallback: clip if repeated negative draws occur
    return max(mean, MIN_POSITIVE_Q)


def main() -> None:
    rng = np.random.default_rng(SEED)

    root = _repo_root()
    data_dir = root / "data"
    out_dir = root / "output" / "simulation"
    data_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    real_path = data_dir / "flume_data.csv"
    if not real_path.exists():
        raise FileNotFoundError(
            f"Could not find {real_path}. Run from repo root, and confirm data/flume_data.csv exists."
        )

    # ----------------------------
    # 1) Read real data for structure
    # ----------------------------
    df_real = pd.read_csv(real_path)

    # Match your R preprocessing: drop EQ rows for the model fit
    df_real = df_real[df_real["flume_type"].astype(str).str.upper() != "EQ"].copy()

    # Standardize flume_type labels (preserve original casing in output if you want)
    df_real["flume_type"] = df_real["flume_type"].astype(str).str.strip()

    # Count points and depth ranges per flume type
    grp = df_real.groupby("flume_type", observed=True)
    flumes = list(grp.groups.keys())

    counts = grp.size().to_dict()
    depth_min = grp["depth_ft"].min().to_dict()
    depth_max = grp["depth_ft"].max().to_dict()

    # ----------------------------
    # OPTIONAL: Inflate N to stress-test structural identifiability
    # ----------------------------

    # Strategy B: Use a fixed N for every flume (overrides multiplier if set)
    N_PER_FLUME = 30        # set to an int like 30, 50, 100
    if N_PER_FLUME is not None:
        counts = {fl: int(N_PER_FLUME) for fl in flumes}

    print("Simulated N per flume:", counts)
    
    # ----------------------------
    # 2) Choose "true" parameters
    # ----------------------------
    # These are intentionally in the neighborhood of your priors in bayesflume.r:
    #   a ~ LogNormal(log(1.55), 0.5)
    #   b ~ LogNormal(log(2.58), 0.5)
    #
    # If you want an even stricter test, set these to values you *expect* from prior fits.

    # Global flow-meter percent-of-reading term (cv_Q)
    # (In your R code: cv_Q ~ LogNormal(log(2e-4), 0.5))
    cv_Q_true = 2.2e-4

    # Per-flume parameters (a, b, sigma_proc, sigma_Q_add)
    # We'll sample them from reasonable distributions so each simulated dataset is distinct,
    # but still plausible and identifiable.
    truth_rows = []
    flume_truth = {}

    for fl in flumes:
        # Draw a and b from lognormal-ish neighborhood (you can hardcode if preferred)
        a_true = float(np.exp(rng.normal(np.log(1.55), 0.25)))
        b_true = float(np.exp(rng.normal(np.log(2.58), 0.15)))

        # Process scatter around rating curve (cfs)
        # In your model: sigma_proc ~ Exponential(1)
        # We'll keep it modest so the curve is still learnable with n~4-6 per flume.
        sigma_proc_true = float(rng.exponential(0.015) + 0.005)  # ~0.005–0.05 typical

        # Additive discharge uncertainty (cfs), per flume
        # In your model: sigma_Q_add ~ LogNormal(log(0.02), 0.5)
        sigma_Q_add_true = float(np.exp(rng.normal(np.log(0.02), 0.35)))

        flume_truth[fl] = dict(
            a=a_true,
            b=b_true,
            sigma_proc=sigma_proc_true,
            sigma_Q_add=sigma_Q_add_true,
        )
        truth_rows.append(
            dict(
                flume_type=fl,
                a_true=a_true,
                b_true=b_true,
                sigma_proc_true=sigma_proc_true,
                sigma_Q_add_true=sigma_Q_add_true,
                cv_Q_true=cv_Q_true,  # repeated for convenience
            )
        )

    # ----------------------------
    # 3) Simulate data following Model 3
    # ----------------------------
    sim_rows = []

    for fl in flumes:
        n = int(counts[fl])
        h_lo = float(depth_min[fl])
        h_hi = float(depth_max[fl])

        if DEPTHS_RANDOM_WITHIN_RANGE:
            H = rng.uniform(h_lo, h_hi, size=n)
        else:
            # Evenly spaced depths within each flume's observed range
            H = np.linspace(h_lo, h_hi, n)

        # Optional jitter to avoid identical depths
        if DEPTH_JITTER_FT > 0:
            H = H + rng.normal(0.0, DEPTH_JITTER_FT, size=n)

        # Enforce positive head
        H = np.maximum(H, MIN_POSITIVE_H)

        a = flume_truth[fl]["a"]
        b = flume_truth[fl]["b"]
        sigma_proc = flume_truth[fl]["sigma_proc"]
        sigma_add = flume_truth[fl]["sigma_Q_add"]

        # Mean curve
        mu = a * (H ** b)

        for i in range(n):
            # Q_true ~ Normal(mu, sigma_proc)
            q_true = _positive_resample_normal(rng, mean=float(mu[i]), sd=float(sigma_proc))

            # Observation SD: sqrt((cv_Q * Q_true)^2 + sigma_Q_add^2)
            sd_obs = math.sqrt((cv_Q_true * q_true) ** 2 + (sigma_add ** 2))

            # Q_obs ~ Normal(Q_true, sd_obs)
            q_obs = _positive_resample_normal(rng, mean=float(q_true), sd=float(sd_obs))

            sim_rows.append(
                dict(
                    depth_ft=float(H[i]),
                    flow_cfs=float(q_obs),
                    flume_type=fl,
                )
            )

    df_sim = pd.DataFrame(sim_rows).sort_values(["flume_type", "depth_ft"]).reset_index(drop=True)

    # ----------------------------
    # 4) Write outputs
    # ----------------------------
    sim_path = data_dir / "sim_data.csv"
    truth_path = out_dir / "sim_truth_model3.csv"

    df_sim.to_csv(sim_path, index=False)
    pd.DataFrame(truth_rows).to_csv(truth_path, index=False)

    # Helpful console output so you don't have to open files to see what happened
    print(f"Wrote simulated data:   {sim_path}")
    print(f"Wrote true parameters:  {truth_path}")
    print("\nTrue parameters used (Model 3 targets):")
    print(pd.DataFrame(truth_rows).to_string(index=False))


if __name__ == "__main__":
    main()
