#!/usr/bin/env python3
"""Reconstruct the real-space two-point function C(theta) from push7's
harmonic-MC F_l(a) data (l_max=8), by resumming the truncated Legendre
series:

    C_resum(theta) = sum_{l=0}^{l_max} (2l+1)/2 * F_l * P_l(cos(theta))

per user request 2026-08-25 ("reconstruct the real space correlators for
the push7 data") -- this is the same resummation used earlier in the
session (journal.md's "BREAKTHROUGH" entry) to show the l_max=8 harmonic
reconstruction disagrees in SHAPE, not just normalization, with the
directly-measured real-space correlator. Re-derived here across the full
push7 ladder (not just one point) and plotted vs. the direct real-space
measurements where available (n_refine=16, 32, journal.md's validated
spot checks -- same physical setup: naive mesh, exact_sinh coupling).

Usage:
  scripts/reconstruct_real_space_from_Fl.py campaign_runs/push7_2026-08-25_dual_naive_eqarea \
    --ladder 8,16,32,64,128 --mode naive
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

sys.path.insert(0, os.path.dirname(__file__))
from cft_symmetry_test import analyze

MODE_SUFFIX = {"naive": "", "equal_area": "_eqarea", "equal_rp": "_eqrp"}

# validated direct real-space measurements (naive mesh, exact_sinh),
# journal.md 2026-08-25 "BREAKTHROUGH" entry -- for overlay comparison
REAL_SPACE_FILES = {
    16: "campaign_runs/real_space_smoke_2026-08-25/q5k16_real_space_2pt_00000001.dat",
    32: "campaign_runs/real_space_2pt_2026-08-25/q5k32_real_space_2pt_0000002A.dat",
}


def pooled_paths(run_dir, k, mode):
    suffix = MODE_SUFFIX.get(mode, "")
    pat = os.path.join(run_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_jackblocks_*.dat")
    return sorted(glob.glob(pat))


def load_real_space(path):
    theta, mean = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            theta.append(float(p[4]))
            mean.append(float(p[5]))
    return np.array(theta), np.array(mean)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--ladder", default="8,16,32,64,128")
    ap.add_argument("--mode", default="naive")
    ap.add_argument("--out", default=None)
    ap.add_argument("--theta-ref", type=float, default=1.5,
                     help="reference angle to normalize all curves to 1.0 at, for shape-only "
                          "comparison (the F_l convention's normalization differs from the "
                          "real-space measurement's by a known ~1/(4pi)^2-ish l-independent "
                          "factor -- irrelevant to the shape question this plot is for)")
    args = ap.parse_args()

    ladder = [int(x) for x in args.ladder.split(",")]
    out_path = args.out or os.path.join(args.run_dir, "plots", "reconstructed_real_space_Fl.png")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    theta_grid = np.linspace(0.05, np.pi - 0.05, 400)
    z_grid = np.cos(theta_grid)
    theta_ref = args.theta_ref

    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0.1, 0.9, len(ladder)))

    print(f"{'n_refine':>9} {'F_0':>10} {'F_1':>10} {'F_1/F_0':>10}")
    for k, color in zip(ladder, colors):
        paths = pooled_paths(args.run_dir, k, args.mode)
        if not paths:
            print(f"n_refine={k}: no jackblocks found, skipping", file=sys.stderr)
            continue
        res = analyze(paths)
        F_l = res["F_full"]
        l_max = res["l_max"]
        print(f"{k:>9} {F_l[0]:>10.4f} {F_l[1]:>10.4f} {F_l[1] / F_l[0]:>10.4f}")

        C_resum = np.zeros_like(theta_grid)
        for l in range(l_max + 1):
            C_resum += (2 * l + 1) / 2.0 * F_l[l] * eval_legendre(l, z_grid)

        # normalize to 1.0 at theta_ref -- the F_l convention's amplitude
        # differs from the direct real-space measurement's by a known
        # l-independent ~1/(4pi)^2-ish factor (journal.md's code audit),
        # irrelevant to the shape question this plot is for.
        resum_ref = np.interp(theta_ref, theta_grid, C_resum)
        C_resum_norm = C_resum / resum_ref
        ax.plot(theta_grid, C_resum_norm, color=color, lw=2, label=f"n_refine={k} (l_max={l_max} resum)")

        if k in REAL_SPACE_FILES and os.path.exists(REAL_SPACE_FILES[k]):
            rs_theta, rs_mean = load_real_space(REAL_SPACE_FILES[k])
            rs_ref = np.interp(theta_ref, rs_theta[::-1], rs_mean[::-1])
            ax.plot(rs_theta, rs_mean / rs_ref, "o", color=color, ms=5, mfc="none", mew=1.5,
                     label=f"n_refine={k} (direct real-space)")

    ax.axvline(theta_ref, color="black", ls=":", lw=1, alpha=0.5)
    ax.set_xlabel("theta")
    ax.set_ylabel(f"C(theta) / C(theta_ref={theta_ref}) -- normalized for shape comparison")
    ax.set_title(f"Reconstructed real-space correlator SHAPE from F_l(a), push7 {args.mode} mesh\n"
                 f"(solid = l_max=8 harmonic resum, open circles = direct real-space measurement;\n"
                 f"amplitudes differ by a known ~1/(4pi)^2 convention factor, normalized out here)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
