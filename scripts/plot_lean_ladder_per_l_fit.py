#!/usr/bin/env python3
"""Continuum-extrapolation fit + plot of Delta_l_pair(a), one panel per
harmonic level l=1..l_max, from a lean_harmonic_stats ladder. Same
Delta(a) = Delta_inf + c*n_refine^-p model as fit_lean_ladder.py (reuses
its fit_power()), applied independently to each l's own Delta_l_pair(a)
and jackknife error (from analyze_lean_harmonic_stats.analyze()) instead
of just Delta_s(a) -- each l is an independent scaling-dimension estimator
via the same one-step recursion inversion (delta_l_pair_formula), so each
gets its own Delta_inf/chi2 rather than assuming they share one.

Usage:
  scripts/plot_lean_ladder_per_l_fit.py <dir_with_*_lean_*.dat files> [--l_max 8] [-o out.png]
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
from analyze_lean_harmonic_stats import analyze  # noqa: E402
from fit_lean_ladder import fit_power  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.data_dir, "plots", "delta_l_pair_continuum_fit.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n.setdefault(n, []).append(f)

    n_refines = []
    dp_all, dp_err_all = [], []  # list of arrays, index l-1
    for n in sorted(by_n):
        r = analyze(sorted(by_n[n]), args.l_max)
        n_refines.append(n)
        dp_all.append(r["delta_l_pair"])
        dp_err_all.append(r["delta_l_pair_err"])
    n_refines = np.array(n_refines, dtype=float)
    dp_all = np.array(dp_all)       # (n_points, l_max)
    dp_err_all = np.array(dp_err_all)

    ncols = 3
    nrows = int(np.ceil(args.l_max / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.3 * nrows), squeeze=False)
    inv_L = 1.0 / n_refines
    xx = np.linspace(0, inv_L.max() * 1.05, 200)

    print(f"{'l':>3} {'p':>2} {'Delta_inf':>12} {'err':>10} {'sigma':>8} {'chi2/dof':>10}")
    for l in range(1, args.l_max + 1):
        ax = axes[(l - 1) // ncols][(l - 1) % ncols]
        y = dp_all[:, l - 1]
        yerr = dp_err_all[:, l - 1]
        color = COLORS[(l - 1) % len(COLORS)]
        ax.errorbar(inv_L, y, yerr=yerr, fmt="o", ms=4, capsize=2, color=color, label="data")

        d1, e1, c1, chi1 = fit_power(n_refines, y, yerr, 1.0)
        d2, e2, c2, chi2v = fit_power(n_refines, y, yerr, 2.0)
        sig1 = (d1 - 0.125) / e1 if e1 == e1 else float("nan")
        sig2 = (d2 - 0.125) / e2 if e2 == e2 else float("nan")
        print(f"{l:3d} {1:2d} {d1:12.5f} {e1:10.5f} {sig1:8.2f} {chi1:10.3f}")
        print(f"{l:3d} {2:2d} {d2:12.5f} {e2:10.5f} {sig2:8.2f} {chi2v:10.3f}")

        ax.plot(xx, d1 + c1 * xx, "-", color="#D55E00", linewidth=1.3,
                 label=rf"p=1: {d1:.4f}({e1*1e4:.0f}), $\chi^2$/dof={chi1:.2f}")
        ax.plot(xx, d2 + c2 * xx ** 2, "--", color="#009E73", linewidth=1.3,
                 label=rf"p=2: {d2:.4f}({e2*1e4:.0f}), $\chi^2$/dof={chi2v:.2f}")
        ax.axhline(0.125, color="gray", linestyle=":", linewidth=1.0)
        ax.set_title(f"l={l}", fontsize=10)
        ax.set_xlabel(r"$1/n_{\rm refine}$", fontsize=8)
        ax.set_ylabel(r"$\Delta_l^{\rm pair}(a)$", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=6, loc="best")

    for idx in range(args.l_max, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle(r"Per-$l$ $\Delta_l^{\rm pair}(a)$ continuum extrapolation")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
