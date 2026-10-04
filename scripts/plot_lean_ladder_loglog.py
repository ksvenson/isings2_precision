#!/usr/bin/env python3
"""Log-linear and log-log diagnostic plots of Delta_l_pair(a) vs
n_refine, one small-multiples figure per style, for all available
harmonic levels l=1..l_max of a lean_harmonic_stats ladder.

Log-linear (semilogx): Delta_l_pair(a) vs n_refine, x on log scale, y
linear -- shows the raw convergence trend and how close each level sits
to the exact 1/8 line as resolution increases, without assuming any
extrapolation model.

Log-log: |Delta_l_pair(a) - 1/8| vs n_refine, both axes log -- a power
law Delta(a) = Delta_inf + c*n^-p shows up as an approximately straight
line here (with Delta_inf=1/8 exactly) with slope -p, so this is a
model-free visual check of the effective local power at each l, directly
comparable to fit_lean_ladder_family.py's fitted p values without
depending on any single fit.

Usage:
  scripts/plot_lean_ladder_loglog.py <dir_with_*_lean_*.dat files> [--l_max 8] [-o out_prefix]
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

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
EXACT = 0.125


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out_prefix", default=None)
    args = ap.parse_args()
    out_prefix = args.out_prefix or os.path.join(args.data_dir, "plots", "delta_l_pair")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n.setdefault(n, []).append(f)

    n_refines = []
    dp_all, dp_err_all = [], []
    for n in sorted(by_n):
        r = analyze(sorted(by_n[n]), args.l_max)
        n_refines.append(n)
        dp_all.append(r["delta_l_pair"])
        dp_err_all.append(r["delta_l_pair_err"])
    n_refines = np.array(n_refines, dtype=float)
    dp_all = np.array(dp_all)
    dp_err_all = np.array(dp_err_all)

    ncols = 3
    nrows = int(np.ceil(args.l_max / ncols))

    # 1) log-linear: Delta_l_pair(a) vs n_refine, x log / y linear
    fig1, axes1 = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.3 * nrows), squeeze=False)
    for l in range(1, args.l_max + 1):
        ax = axes1[(l - 1) // ncols][(l - 1) % ncols]
        y = dp_all[:, l - 1]
        yerr = dp_err_all[:, l - 1]
        color = COLORS[(l - 1) % len(COLORS)]
        ax.errorbar(n_refines, y, yerr=yerr, fmt="o-", ms=4, capsize=2, color=color)
        ax.axhline(EXACT, color="gray", linestyle=":", linewidth=1.0, label="exact 1/8")
        ax.set_xscale("log")
        ax.set_title(f"l={l}", fontsize=10)
        ax.set_xlabel(r"$n_{\rm refine}$ (log)", fontsize=8)
        ax.set_ylabel(r"$\Delta_l^{\rm pair}(a)$", fontsize=8)
        style_axes(ax)
        ax.legend(frameon=False, fontsize=6)
    for idx in range(args.l_max, nrows * ncols):
        axes1[idx // ncols][idx % ncols].axis("off")
    fig1.suptitle(r"$\Delta_l^{\rm pair}(a)$ vs $n_{\rm refine}$ (log-linear)")
    fig1.tight_layout()
    p1 = out_prefix + "_loglin.png"
    fig1.savefig(p1, dpi=150)
    plt.close(fig1)
    print(f"wrote {p1}", file=sys.stderr)

    # 2) log-log: |Delta_l_pair(a) - 1/8| vs n_refine, both axes log
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.3 * nrows), squeeze=False)
    xx = np.geomspace(n_refines.min(), n_refines.max(), 100)
    for l in range(1, args.l_max + 1):
        ax = axes2[(l - 1) // ncols][(l - 1) % ncols]
        y = dp_all[:, l - 1]
        yerr = dp_err_all[:, l - 1]
        dev = y - EXACT
        # error on |dev| is just yerr (dev's own error), sign doesn't affect |.|'s local error to first order
        pos = dev > 0
        color = COLORS[(l - 1) % len(COLORS)]
        ax.errorbar(n_refines[pos], np.abs(dev[pos]), yerr=yerr[pos], fmt="o", ms=5, capsize=2,
                    color=color, label=r"$\Delta(a)-1/8 > 0$")
        ax.errorbar(n_refines[~pos], np.abs(dev[~pos]), yerr=yerr[~pos], fmt="s", ms=5, capsize=2,
                    color=color, mfc="white", label=r"$\Delta(a)-1/8 < 0$")
        # reference slope guides through the finest-resolution point
        anchor_n, anchor_y = n_refines[-1], np.abs(dev[-1])
        if anchor_y > 0:
            for p, ls in ((1, "--"), (2, ":")):
                ax.plot(xx, anchor_y * (xx / anchor_n) ** (-p), ls, color="gray", linewidth=1.0, alpha=0.7,
                        label=f"slope -{p} ref" if l == 1 else None)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"l={l}", fontsize=10)
        ax.set_xlabel(r"$n_{\rm refine}$ (log)", fontsize=8)
        ax.set_ylabel(r"$|\Delta_l^{\rm pair}(a)-1/8|$ (log)", fontsize=8)
        style_axes(ax)
        ax.legend(frameon=False, fontsize=6)
    for idx in range(args.l_max, nrows * ncols):
        axes2[idx // ncols][idx % ncols].axis("off")
    fig2.suptitle(r"$|\Delta_l^{\rm pair}(a)-1/8|$ vs $n_{\rm refine}$ (log-log; slope $\approx -p$)")
    fig2.tight_layout()
    p2 = out_prefix + "_loglog.png"
    fig2.savefig(p2, dpi=150)
    plt.close(fig2)
    print(f"wrote {p2}", file=sys.stderr)


if __name__ == "__main__":
    main()
