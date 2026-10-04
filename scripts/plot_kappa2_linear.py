#!/usr/bin/env python3
"""Linear-linear (not log-log) naive-vs-equal_area kappa2_r comparison vs
1/L -- companion to plot_mesh_compare_full.py's log-log kappa2_r panel,
for reading off the plateau's absolute shape/spacing rather than its
scaling behavior. Same underlying data/reduction
(symmetry_check_lean_kappa2.point_kappa2), just linear axes.

Usage:
  scripts/plot_kappa2_linear.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out_dir]
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
from symmetry_check_lean_kappa2 import point_kappa2  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
COLOR_NAIVE = COLORS[0]
COLOR_EQAREA = COLORS[1]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)


def collect(data_dir):
    files = glob.glob(os.path.join(data_dir, "q5k*_lean_*.dat"))
    by_n = {}
    for f in files:
        m = re.search(r"q5k(\d+)(?:_eqarea|_eqrp)?_lean", f)
        by_n.setdefault(int(m.group(1)), []).append(f)
    return by_n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("naive_dir")
    ap.add_argument("eqarea_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--min_n_refine", type=int, default=32,
                     help="only plot ladder points with n_refine strictly greater than this (default 32)")
    ap.add_argument("-o", "--out_dir", default=None)
    args = ap.parse_args()

    naive_by_n = collect(args.naive_dir)
    eqarea_by_n = collect(args.eqarea_dir)
    ns = sorted(n for n in (set(naive_by_n) & set(eqarea_by_n)) if n > args.min_n_refine)
    inv_L = 1.0 / np.array(ns, dtype=float)

    out_dir = args.out_dir or os.path.join(os.path.dirname(args.naive_dir.rstrip("/")), "plots")
    os.makedirs(out_dir, exist_ok=True)

    ls = list(range(3, args.l_max + 1))
    naive_rows, eqarea_rows = {}, {}
    for n in ns:
        _, _, rn = point_kappa2(sorted(naive_by_n[n]), args.l_max)
        _, _, re_ = point_kappa2(sorted(eqarea_by_n[n]), args.l_max)
        naive_rows[n] = {r["l"]: r for r in rn}
        eqarea_rows[n] = {r["l"]: r for r in re_}

    ncols = 3
    nrows = int(np.ceil(len(ls) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for i, l in enumerate(ls):
        ax = axes[i]
        naive_vals = np.array([naive_rows[n][l]["kappa2_r"] for n in ns], dtype=float)
        naive_errs = np.array([naive_rows[n][l]["kappa2_r_err"] for n in ns], dtype=float)
        eq_vals = np.array([eqarea_rows[n][l]["kappa2_r"] for n in ns], dtype=float)
        eq_errs = np.array([eqarea_rows[n][l]["kappa2_r_err"] for n in ns], dtype=float)
        ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=COLOR_NAIVE,
                    label="naive", ms=5, linewidth=1.6, capsize=3)
        ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=COLOR_EQAREA,
                    label="equal_area", ms=5, linewidth=1.6, capsize=3)
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
        ax.set_title(f"l = {l}", fontsize=11)
        style_axes(ax)
        ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
        if i % ncols == 0:
            ax.set_ylabel(r"$\kappa_{2,r}$")
        if i == 0:
            ax.legend(fontsize=9, frameon=False)
    for j in range(len(ls), len(axes)):
        axes[j].axis("off")
    fig.suptitle(
        r"$\kappa_{2,r}$ (linear-linear axes): naive vs equal_area "
        f"({min(ns)} $\\leq$ n_refine $\\leq$ {max(ns)})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(out_dir, "mesh_compare_kappa2_r_vs_invL_linear.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
