#!/usr/bin/env python3
"""Two-panel kappa2_r plot: left = naive (all l), right = equal_area (all l).

Just the kappa2_r half of symmetry_check_lean_kappa2.py's output (no R_l),
naive and equal_area side by side instead of in two separate files.

Usage:
  scripts/plot_kappa2r_naive_vs_eqarea_alll.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out.png]
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symmetry_check_lean_kappa2 import point_kappa2, COLORS  # noqa: E402


def collect(data_dir, l_max):
    files = glob.glob(os.path.join(data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)(?:_eqarea|_eqrp)?_lean", f).group(1))
        by_n.setdefault(n, []).append(f)
    n_refines = sorted(by_n)
    kappa2r = np.full((len(n_refines), l_max + 1), np.nan)
    kappa2r_err = np.full((len(n_refines), l_max + 1), np.nan)
    for i, n in enumerate(n_refines):
        _, _, results = point_kappa2(sorted(by_n[n]), l_max)
        for r in results:
            l = r["l"]
            if l == 0:
                continue
            kappa2r[i, l] = r["kappa2_r"]
            kappa2r_err[i, l] = r["kappa2_r_err"]
    return np.array(n_refines, dtype=float), kappa2r, kappa2r_err


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("naive_dir")
    ap.add_argument("eqarea_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out", default="kappa2r_naive_vs_eqarea_alll.png")
    args = ap.parse_args()

    n_naive, k_naive, k_naive_err = collect(args.naive_dir, args.l_max)
    n_eqarea, k_eqarea, k_eqarea_err = collect(args.eqarea_dir, args.l_max)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for l in range(1, args.l_max + 1):
        color = COLORS[(l - 1) % len(COLORS)]
        axes[0].errorbar(1.0 / n_naive, k_naive[:, l], yerr=k_naive_err[:, l],
                          fmt="o-", ms=4, color=color, label=f"l={l}")
        axes[1].errorbar(1.0 / n_eqarea, k_eqarea[:, l], yerr=k_eqarea_err[:, l],
                          fmt="o-", ms=4, color=color, label=f"l={l}")

    axes[0].set_title("naive")
    axes[1].set_title("equal_area")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$1/n_{\rm refine}$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=8, ncol=2)
    axes[0].set_ylabel(r"$\kappa_{2,r}$")
    fig.suptitle(r"Real $\kappa_{2,r}$ (off-diagonal $M_l$), all $l$")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
