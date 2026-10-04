#!/usr/bin/env python3
"""Plot bin/lean_harmonic_stats observables vs 1/n_refine across a ladder
of independent runs. Reuses analyze_lean_harmonic_stats.analyze() for the
per-point Delta_s/Delta_l_pair/<M^2>/Binder-U4 reduction (with jackknife
errors) -- no reimplementation of the F_l/(2l+1) physics here.

Usage:
  scripts/plot_lean_ladder.py <dir_with_*_lean_*.dat files> --l_max 8 -o <out_dir>
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

# palette: colorblind-safe categorical set (Okabe-Ito), consistent across plots
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out_dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or os.path.join(args.data_dir, "plots")
    os.makedirs(out_dir, exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")

    # group shard files by n_refine (multiple independent-seed shards per
    # point are pooled by analyze(), same as this campaign's other
    # multi-shard ladder analyses -- see cluster/sge/lean_harmonic_stats_task.sh)
    by_n_refine = {}
    for f in files:
        n_refine = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n_refine.setdefault(n_refine, []).append(f)

    results = []
    for n_refine in sorted(by_n_refine):
        shard_files = sorted(by_n_refine[n_refine])
        r = analyze(shard_files, args.l_max)
        r["n_refine"] = n_refine
        results.append(r)
        print(f"n_refine={n_refine:3d}  n_shards={len(shard_files)}  "
              f"Delta_s={r['delta_s']:.5f}+/-{r['delta_s_err']:.5f}  "
              f"<M^2>={r['m2']:.4f}  U4={r['u4']:.4f}", file=sys.stderr)

    inv_L = np.array([1.0 / r["n_refine"] for r in results])

    # 1) Delta_s(a) vs 1/L
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ds = np.array([r["delta_s"] for r in results])
    ds_err = np.array([r["delta_s_err"] for r in results])
    ax.errorbar(inv_L, ds, yerr=ds_err, fmt="o-", color=COLORS[0], ms=5, capsize=3, label=r"$\Delta_s(a)$ (lean sim, Trace($M_l$) method)")
    ax.axhline(0.125, color=COLORS[1], linestyle="--", linewidth=1.2, label=r"exact $1/8$")
    ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
    ax.set_ylabel(r"$\Delta_s(a)$")
    ax.set_title(r"$\Delta_\sigma$ continuum extrapolation")
    style_axes(ax)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "delta_s_vs_invL.png"), dpi=150)
    plt.close(fig)

    # 2) Delta_l_pair(a) vs 1/L, one line per l
    fig, ax = plt.subplots(figsize=(6.5, 5))
    l_max_common = min(r["l_max"] for r in results)
    for l in range(1, l_max_common + 1):
        y = np.array([r["delta_l_pair"][l - 1] for r in results])
        yerr = np.array([r["delta_l_pair_err"][l - 1] for r in results])
        color = COLORS[(l - 1) % len(COLORS)]
        ax.errorbar(inv_L, y, yerr=yerr, fmt="o-", ms=4, capsize=2, color=color, label=f"l={l}")
    ax.axhline(0.125, color="gray", linestyle="--", linewidth=1.0, label=r"exact $1/8$")
    ax.set_xlabel(r"$1/L$")
    ax.set_ylabel(r"$\Delta_l^{\rm pair}(a)$")
    ax.set_title(r"Per-$l$ scaling-dimension estimator (single-step recursion inversion)")
    style_axes(ax)
    ax.legend(frameon=False, fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "delta_l_pair_vs_invL.png"), dpi=150)
    plt.close(fig)

    # 3) F_l(a) vs 1/L, one panel per l (small multiples)
    n_panels = l_max_common + 1
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    for l in range(n_panels):
        ax = axes[l // ncols][l % ncols]
        y = np.array([r["F_full"][l] for r in results])
        yerr = np.array([r["F_err"][l] for r in results])
        ax.errorbar(inv_L, y, yerr=yerr, fmt="o-", ms=4, capsize=2, color=COLORS[l % len(COLORS)])
        ax.set_title(f"l={l}", fontsize=10)
        ax.set_xlabel(r"$1/L$", fontsize=8)
        style_axes(ax)
    for l in range(n_panels, nrows * ncols):
        axes[l // ncols][l % ncols].axis("off")
    fig.suptitle(r"$F_l(a)$ (fixed, $(2l{+}1)$-divided) vs $1/L$")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "F_l_vs_invL.png"), dpi=150)
    plt.close(fig)

    # 4) <M^2> and Binder U4 vs 1/L
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    m2 = np.array([r["m2"] for r in results])
    m2_err = np.array([r["m2_err"] for r in results])
    u4 = np.array([r["u4"] for r in results])
    u4_err = np.array([r["u4_err"] for r in results])
    axes[0].errorbar(inv_L, m2, yerr=m2_err, fmt="o-", color=COLORS[0], ms=5, capsize=3)
    axes[0].set_xlabel(r"$1/L$")
    axes[0].set_ylabel(r"$\langle M^2 \rangle$")
    axes[0].set_title(r"Mean-square magnetization")
    style_axes(axes[0])
    axes[1].errorbar(inv_L, u4, yerr=u4_err, fmt="o-", color=COLORS[2], ms=5, capsize=3)
    axes[1].set_xlabel(r"$1/L$")
    axes[1].set_ylabel(r"$U_4 = 1 - \langle M^4\rangle/(3\langle M^2\rangle^2)$")
    axes[1].set_title(r"Binder cumulant")
    style_axes(axes[1])
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "magnetization_vs_invL.png"), dpi=150)
    plt.close(fig)

    print(f"wrote plots to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
