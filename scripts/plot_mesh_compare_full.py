#!/usr/bin/env python3
"""Full naive-vs-equal_area comparison, all measurement metrics, vs 1/L.

Produces:
  - one small-multiples figure per symmetry-breaking metric (R_l, kappa2_r,
    kappa3_r, kappa4_r), panels = l=3..l_max, x-axis = 1/L (=1/n_refine)
  - one figure for the scalar CFT-data/bulk metrics (Delta_s(a), Binder U4),
    x-axis = 1/L
  - one ratio-summary figure (naive/eqarea) across all four symmetry metrics

All from the SAME underlying *_lean_*.dat files: the per-l symmetry metrics
come from the real off-diagonal M_l (symmetry_check_lean_kappa2.point_kappa2
-> symmetry_test.analyze_blocks), the scalar metrics come from the diagonal
|S_lm|^2 + magnetization-moment reduction (analyze_lean_harmonic_stats.analyze).
Both reuse the campaign's existing, already-verified reduction code -- no
reimplementation. Style (palette, axis conventions, style_axes) matches
scripts/plot_lean_ladder.py so every plot in this campaign reads as one
system.

Usage:
  scripts/plot_mesh_compare_full.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out_dir]
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
from analyze_lean_harmonic_stats import analyze as analyze_scalar  # noqa: E402

# palette: colorblind-safe categorical set (Okabe-Ito), consistent across
# every plot script in this campaign (see plot_lean_ladder.py)
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
COLOR_NAIVE = COLORS[0]
COLOR_EQAREA = COLORS[1]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)


SYMMETRY_METRICS = [
    ("R_l", "R_l_err", r"$R_l$"),
    ("kappa2_r", "kappa2_r_err", r"$\kappa_{2,r}$"),
    ("kappa3_r", "kappa3_r_err", r"$\kappa_{3,r}$"),
    ("kappa4_r", "kappa4_r_err", r"$\kappa_{4,r}$"),
]


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
    ap.add_argument("-o", "--out_dir", default=None)
    args = ap.parse_args()

    naive_by_n = collect(args.naive_dir)
    eqarea_by_n = collect(args.eqarea_dir)
    ns = sorted(set(naive_by_n) & set(eqarea_by_n))
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
        print(f"n_refine={n:3d}  " + "  ".join(
            f"l={l}:naive_k2={naive_rows[n][l]['kappa2_r']:.2e},eq_k2={eqarea_rows[n][l]['kappa2_r']:.2e}"
            for l in ls), file=sys.stderr)

    # ---- one small-multiples figure per symmetry-breaking metric ----
    for metric, metric_err, ylabel in SYMMETRY_METRICS:
        ncols = 3
        nrows = int(np.ceil(len(ls) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), sharex=True)
        axes = np.atleast_1d(axes).ravel()
        for i, l in enumerate(ls):
            ax = axes[i]
            naive_vals = np.array([naive_rows[n][l][metric] for n in ns], dtype=float)
            naive_errs = np.array([naive_rows[n][l][metric_err] for n in ns], dtype=float)
            eq_vals = np.array([eqarea_rows[n][l][metric] for n in ns], dtype=float)
            eq_errs = np.array([eqarea_rows[n][l][metric_err] for n in ns], dtype=float)
            ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=COLOR_NAIVE,
                        label="naive", ms=5, linewidth=1.6, capsize=3)
            ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=COLOR_EQAREA,
                        label="equal_area", ms=5, linewidth=1.6, capsize=3)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(f"l = {l}", fontsize=11)
            style_axes(ax)
            ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
            if i % ncols == 0:
                ax.set_ylabel(ylabel)
            if i == 0:
                ax.legend(fontsize=9, frameon=False)
        for j in range(len(ls), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"{ylabel}: naive vs equal_area vs $1/L$ (n_refine $\\leq$ {max(ns)})", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        out = os.path.join(out_dir, f"mesh_compare_{metric}_vs_invL.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"wrote {out}")

    # ---- scalar metrics: Delta_s(a), Binder U4 ----
    naive_scalar, eqarea_scalar = [], []
    for n in ns:
        naive_scalar.append(analyze_scalar(sorted(naive_by_n[n]), l_max=args.l_max))
        eqarea_scalar.append(analyze_scalar(sorted(eqarea_by_n[n]), l_max=args.l_max))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ds_n = np.array([r["delta_s"] for r in naive_scalar])
    ds_n_err = np.array([r["delta_s_err"] for r in naive_scalar])
    ds_e = np.array([r["delta_s"] for r in eqarea_scalar])
    ds_e_err = np.array([r["delta_s_err"] for r in eqarea_scalar])
    ax.errorbar(inv_L, ds_n, yerr=ds_n_err, fmt="o-", color=COLOR_NAIVE, label="naive",
                ms=5, linewidth=1.6, capsize=3)
    ax.errorbar(inv_L, ds_e, yerr=ds_e_err, fmt="s-", color=COLOR_EQAREA, label="equal_area",
                ms=5, linewidth=1.6, capsize=3)
    ax.axhline(0.125, color=COLORS[2], linestyle="--", linewidth=1.2, label=r"exact $1/8$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
    ax.set_ylabel(r"$\Delta_s(a)$")
    ax.set_title(r"$\Delta_s(a) = 2F_1/(F_1+F_0)$")
    style_axes(ax)
    ax.legend(fontsize=9, frameon=False)

    ax = axes[1]
    u4_n = np.array([r["u4"] for r in naive_scalar])
    u4_n_err = np.array([r["u4_err"] for r in naive_scalar])
    u4_e = np.array([r["u4"] for r in eqarea_scalar])
    u4_e_err = np.array([r["u4_err"] for r in eqarea_scalar])
    ax.errorbar(inv_L, u4_n, yerr=u4_n_err, fmt="o-", color=COLOR_NAIVE, label="naive",
                ms=5, linewidth=1.6, capsize=3)
    ax.errorbar(inv_L, u4_e, yerr=u4_e_err, fmt="s-", color=COLOR_EQAREA, label="equal_area",
                ms=5, linewidth=1.6, capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
    ax.set_ylabel(r"Binder $U_4$")
    ax.set_title(r"Binder cumulant $U_4 = 1 - \langle M^4\rangle/(3\langle M^2\rangle^2)$")
    style_axes(ax)
    ax.legend(fontsize=9, frameon=False)

    fig.suptitle(f"Scalar bulk/CFT-data metrics: naive vs equal_area (n_refine $\\leq$ {max(ns)})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = os.path.join(out_dir, "mesh_compare_scalar_vs_invL.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")

    # ---- ratio summary (naive/eqarea) across all four symmetry metrics, l=3..8 ----
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
    axes = axes.ravel()
    l_colors = {l: COLORS[i % len(COLORS)] for i, l in enumerate(ls)}
    for mi, (metric, metric_err, ylabel) in enumerate(SYMMETRY_METRICS):
        ax = axes[mi]
        for l in ls:
            naive_vals = np.array([naive_rows[n][l][metric] for n in ns], dtype=float)
            naive_errs = np.array([naive_rows[n][l][metric_err] for n in ns], dtype=float)
            eq_vals = np.array([eqarea_rows[n][l][metric] for n in ns], dtype=float)
            eq_errs = np.array([eqarea_rows[n][l][metric_err] for n in ns], dtype=float)
            ratio = naive_vals / eq_vals
            # kappa3_r can be negative (odd cumulant) -- error magnitude must
            # stay positive regardless of ratio's sign
            ratio_err = np.abs(ratio) * np.sqrt((naive_errs / naive_vals) ** 2 + (eq_errs / eq_vals) ** 2)
            ax.errorbar(inv_L, ratio, yerr=ratio_err, fmt="o-", color=l_colors[l], label=f"l={l}",
                        ms=4, linewidth=1.3, capsize=2)
        ax.axhline(1.0, color="gray", linewidth=1, linestyle="--", alpha=0.7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
        ax.set_ylabel(f"{ylabel} ratio (naive/eqarea)")
        style_axes(ax)
        if mi == 0:
            ax.legend(fontsize=8, frameon=False, ncol=2)
    fig.suptitle("Mesh-mode ratio across all four symmetry-breaking metrics", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(out_dir, "mesh_compare_ratio_all_metrics_vs_invL.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
