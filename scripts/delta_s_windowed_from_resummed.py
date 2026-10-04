#!/usr/bin/env python3
"""Windowed Delta_s(a) from the l_max=8 harmonic-resummed real-space
correlator (reconstruct_real_space_from_Fl.py's C_resum(theta)), across
push7's full ladder -- per user request 2026-08-25 ("get the delta_s(a)
from F0 F1 from the corrected integration method from the real space
correlators").

This tests whether the l_max=8 truncated resummation's known shape
distortion (journal.md's "BREAKTHROUGH" entry: factor ~100-250x
disagreement vs. the direct measurement, worst at small theta, not a
constant ratio) is confined enough to short distances that windowing it
out (same theta_min/theta_max cut as fit_real_space_2pt.py /
windowed_f0_f1_delta.py) recovers a Delta close to 1/8 -- i.e. whether
"reconstruct then window" gives the same answer as "measure directly then
window".

Method, reusing windowed_f0_f1_delta.py's model side unchanged:
  1. For each ladder point, compute F_l(a) (l=0..8) via
     cft_symmetry_test.analyze (same jackblocks pooling as
     run_cft_ladder.py).
  2. Resum to C_resum(theta) on a fine grid (same formula as
     reconstruct_real_space_from_Fl.py).
  3. Compute windowed F_0,F_1 from C_resum via trapezoidal integration
     over z=cos(theta) in [cos(theta_max), cos(theta_min)] (same window
     convention as windowed_f0_f1_delta.py).
  4. Solve for Delta matching the amplitude-independent ratio
     2*F_1/(F_1+F_0) against the model (scipy.integrate.quad, no
     recursion).
  5. Plot this "windowed-resummed" Delta_s(a) vs 1/n_refine next to the
     plain unwindowed Delta_s(a) (already computed by run_cft_ladder.py).

Usage:
  scripts/delta_s_windowed_from_resummed.py campaign_runs/push7_2026-08-25_dual_naive_eqarea \
    --ladder 2,3,4,6,8,12,16,24,32,48,64,96,128 --mode naive
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import eval_legendre
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
from cft_symmetry_test import analyze
from windowed_f0_f1_delta import model_Fl

MODE_SUFFIX = {"naive": "", "equal_area": "_eqarea", "equal_rp": "_eqrp"}


def pooled_paths(run_dir, k, mode):
    suffix = MODE_SUFFIX.get(mode, "")
    pat = os.path.join(run_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_jackblocks_*.dat")
    return sorted(glob.glob(pat))


def windowed_delta_from_curve(theta_grid, C, theta_min, theta_max):
    mask = (theta_grid >= theta_min) & (theta_grid <= theta_max)
    z = np.cos(theta_grid[mask])
    order = np.argsort(z)
    z_sorted = z[order]
    C_sorted = C[mask][order]

    F0_data = np.trapz(C_sorted, z_sorted)
    F1_data = np.trapz(C_sorted * z_sorted, z_sorted)
    if F1_data + F0_data == 0:
        return None
    ratio_data = 2.0 * F1_data / (F1_data + F0_data)

    z_min, z_max = np.cos(theta_max), np.cos(theta_min)

    def resid(Delta):
        F0m = model_Fl(Delta, z_min, z_max, 0)
        F1m = model_Fl(Delta, z_min, z_max, 1)
        return 2.0 * F1m / (F1m + F0m) - ratio_data

    lo, hi = 1e-3, 0.99
    if resid(lo) * resid(hi) > 0:
        return None
    return brentq(resid, lo, hi, xtol=1e-6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--ladder", default="2,3,4,6,8,12,16,24,32,48,64,96,128")
    ap.add_argument("--mode", default="naive")
    ap.add_argument("--theta-min", type=float, default=0.3)
    ap.add_argument("--theta-max", type=float, default=2.8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ladder = [int(x) for x in args.ladder.split(",")]
    out_path = args.out or os.path.join(args.run_dir, "plots", "delta_s_windowed_vs_unwindowed.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    theta_grid = np.linspace(0.02, np.pi - 0.02, 2000)
    z_grid = np.cos(theta_grid)

    invk, ds_unwindowed, ds_windowed = [], [], []
    print(f"{'n_refine':>9} {'Delta_s_unwindowed':>20} {'Delta_s_windowed_resum':>24}")
    for k in ladder:
        paths = pooled_paths(args.run_dir, k, args.mode)
        if not paths:
            print(f"n_refine={k}: no jackblocks found, skipping", file=sys.stderr)
            continue
        res = analyze(paths)
        F_l = res["F_full"]
        l_max = res["l_max"]

        C_resum = np.zeros_like(theta_grid)
        for l in range(l_max + 1):
            C_resum += (2 * l + 1) / 2.0 * F_l[l] * eval_legendre(l, z_grid)

        d_win = windowed_delta_from_curve(theta_grid, C_resum, args.theta_min, args.theta_max)
        d_unwin = res["delta_s"]

        invk.append(1.0 / k)
        ds_unwindowed.append(d_unwin)
        ds_windowed.append(d_win if d_win is not None else np.nan)
        print(f"{k:>9} {d_unwin:>20.5f} {'--' if d_win is None else f'{d_win:.5f}':>24}")

    invk = np.array(invk)
    ds_unwindowed = np.array(ds_unwindowed)
    ds_windowed = np.array(ds_windowed)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(invk, ds_unwindowed, "o-", color="#4C72B0", label="Delta_s(a) (unwindowed, full sphere)")
    ax.plot(invk, ds_windowed, "s-", color="#DD8452", label=f"Delta_s(a) (windowed resum, theta in [{args.theta_min},{args.theta_max}])")
    ax.axhline(0.125, color="gray", ls="--", lw=1, label="exact 1/8")
    ax.set_xlabel("1/n_refine")
    ax.set_ylabel("Delta_s(a)")
    ax.set_title(f"Delta_s(a) vs 1/n_refine, push7 {args.mode} mesh:\nunwindowed vs. windowed-resummed F0/F1")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
