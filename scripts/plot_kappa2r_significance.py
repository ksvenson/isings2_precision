#!/usr/bin/env python3
"""Two-panel significance (z = kappa2_r/err) plot: left = naive, right =
equal_area, all l. Answers "is the l>=3 signal distinguishable from the
l=1,2 MC-noise floor" directly, instead of inferring it from central
values on a log axis where large error bars are easy to miss.

Usage:
  scripts/plot_kappa2r_significance.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out.png]
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
    z = np.full((len(n_refines), l_max + 1), np.nan)
    for i, n in enumerate(n_refines):
        _, _, results = point_kappa2(sorted(by_n[n]), l_max)
        for r in results:
            l = r["l"]
            if l == 0:
                continue
            z[i, l] = r["kappa2_r"] / r["kappa2_r_err"] if r["kappa2_r_err"] > 0 else np.nan
    return np.array(n_refines, dtype=float), z


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("naive_dir")
    ap.add_argument("eqarea_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out", default="kappa2r_significance.png")
    ap.add_argument("--ymax", type=float, default=None,
                     help="clip y-axis (the small-n_refine points can be 10-100x larger, "
                          "swamping the noise-floor region otherwise)")
    args = ap.parse_args()

    n_naive, z_naive = collect(args.naive_dir, args.l_max)
    n_eqarea, z_eqarea = collect(args.eqarea_dir, args.l_max)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for l in range(1, args.l_max + 1):
        color = COLORS[(l - 1) % len(COLORS)]
        axes[0].plot(1.0 / n_naive, z_naive[:, l], "o-", ms=4, color=color, label=f"l={l}")
        axes[1].plot(1.0 / n_eqarea, z_eqarea[:, l], "o-", ms=4, color=color, label=f"l={l}")

    axes[0].set_title("naive")
    axes[1].set_title("equal_area")
    for ax in axes:
        ax.set_xscale("log")
        ax.axhspan(0, 2, color="gray", alpha=0.15, zorder=0)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xlabel(r"$1/n_{\rm refine}$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=8, ncol=2)
    if args.ymax is not None:
        axes[0].set_ylim(top=args.ymax)
    axes[0].set_ylabel(r"$z = \kappa_{2,r} / \sigma(\kappa_{2,r})$")
    fig.suptitle(r"$\kappa_{2,r}$ significance above zero (shaded band = $z<2$, the l=1,2 noise floor)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
