#!/usr/bin/env python3
"""
main.py

Run all project plots with unified styling:
  1) flume-analysis.py
  2) flume-analysis-with-band.py
  3) flume-montecarlo.py

By default, running main.py produces all three outputs.
"""

from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CODE = HERE
DATA = REPO / "data"
OUTPUT = REPO / "output"

DEFAULT_CSV = DATA / "flume_data.csv"
OUT_BASIC   = OUTPUT / "python_measurement_error" / "flume_stage_discharge_with_errors.png"
OUT_BAND    = OUTPUT / "python_measurement_error" / "flume_stage_discharge_with_band.png"
OUT_MC      = OUTPUT / "python_monte_carlo" / "flume_montecarlo_band.png"

SCRIPT_BASIC = CODE / "flume-analysis.py"
SCRIPT_BAND  = CODE / "flume-analysis-with-band.py"
SCRIPT_MC    = CODE / "flume-montecarlo.py"

def run_script(script: Path, args: list[str]) -> None:
    if not script.exists():
        raise FileNotFoundError(f"Missing helper script: {script}")
    cmd = [sys.executable, str(script)] + args
    print(f"[RUN] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def str_path(p: Path | str) -> str:
    return str(p)

def main():
    ap = argparse.ArgumentParser(description="Archived initial Python workflow; unused by the current Bayesian analysis")
    ap.add_argument("--csv", default=str_path(DEFAULT_CSV), help="Input CSV path")
    ap.add_argument("--alpha", type=float, default=0.18, help="Alpha for shaded bands")
    ap.add_argument("--n", type=int, default=1000, help="Monte Carlo iterations")
    ap.add_argument("--manuf-band", action="store_true", help="Shaded band for manufacturer curve in band plot")
    ap.add_argument("--no-errorbars", dest="no_errorbars", action="store_true",
                    help="Suppress per-point error bars in band plot")
    ap.add_argument("--only", choices=["basic","band","mc","all"], default="all",
                    help="Limit to a specific product, default runs all")
    args = ap.parse_args()

    OUTPUT.mkdir(parents=True, exist_ok=True)

    to_run = []
    if args.only in ("basic", "all"):
        to_run.append((SCRIPT_BASIC, [str_path(args.csv), "--out", str_path(OUT_BASIC)]))
    if args.only in ("band", "all"):
        band_args = [str_path(args.csv), "--out", str_path(OUT_BAND), "--band", "--alpha", str(args.alpha)]
        if args.no_errorbars:
            band_args.append("--no-errorbars")
        if args.manuf_band:
            band_args.append("--manuf-band")
        to_run.append((SCRIPT_BAND, band_args))
    if args.only in ("mc", "all"):
        to_run.append((SCRIPT_MC, [str_path(args.csv), "--out", str_path(OUT_MC), "--n", str(args.n)]))

    try:
        for script, argv in to_run:
            run_script(script, argv)
        print("[OK] All requested plots generated.")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Subprocess failed with exit code {e.returncode}: {e}")
        sys.exit(e.returncode)
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        sys.exit(2)

if __name__ == "__main__":
    main()
