#!/usr/bin/env python3
"""naive-vs-equal_area comparison of ALL magnetization observables, vs 1/L.

Companion to plot_mesh_compare_full.py (which only plots Delta_s(a) and
Binder U4 from the scalar/magnetization reduction) -- this script plots the
rest of what analyze_lean_harmonic_stats.analyze() computes from the same
per-jackknife-block M/M2/M3/M4 sums: <M>, <M^2>, <M^3>, <M^4>, susceptibility
chi = (<M^2>-<M>^2)/n_sites, and Binder U4, each naive vs equal_area vs 1/L.

<M> and <M^3> are odd moments, expected ~0 by the model's exact Z2 symmetry
(cold/hot start symmetric MC does not spontaneously break it in these run
lengths) -- plotted anyway as a direct check that this null expectation
actually holds down to the level of the jackknife error bars, on both
meshes, rather than assumed.

Usage:
  scripts/plot_mesh_compare_magnetization.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out_dir]
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
from analyze_lean_harmonic_stats import analyze as analyze_scalar  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
COLOR_NAIVE = COLORS[0]
COLOR_EQAREA = COLORS[1]

# (result dict key, error key, ylabel, title, yscale, zero-reference-line)
METRICS = [
    ("m1", "m1_err", r"$\langle M\rangle$", r"$\langle M\rangle$ (odd moment, expect 0)", "linear", True),
    ("m2", "m2_err", r"$\langle M^2\rangle$", r"$\langle M^2\rangle$", "log", False),
    ("m3", "m3_err", r"$\langle M^3\rangle$", r"$\langle M^3\rangle$ (odd moment, expect 0)", "linear", True),
    ("m4", "m4_err", r"$\langle M^4\rangle$", r"$\langle M^4\rangle$", "log", False),
    ("chi", "chi_err", r"$\chi = (\langle M^2\rangle-\langle M\rangle^2)/n_{\rm sites}$",
     "Susceptibility per site", "log", False),
    ("u4", "u4_err", r"Binder $U_4$", r"$U_4 = 1-\langle M^4\rangle/(3\langle M^2\rangle^2)$", "linear", False),
]


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
    ap.add_argument("-o", "--out_dir", default=None)
    args = ap.parse_args()

    naive_by_n = collect(args.naive_dir)
    eqarea_by_n = collect(args.eqarea_dir)
    ns = sorted(set(naive_by_n) & set(eqarea_by_n))
    inv_L = 1.0 / np.array(ns, dtype=float)

    out_dir = args.out_dir or os.path.join(os.path.dirname(args.naive_dir.rstrip("/")), "plots")
    os.makedirs(out_dir, exist_ok=True)

    naive_res, eqarea_res = {}, {}
    for n in ns:
        naive_res[n] = analyze_scalar(sorted(naive_by_n[n]), l_max=args.l_max)
        eqarea_res[n] = analyze_scalar(sorted(eqarea_by_n[n]), l_max=args.l_max)
        print(
            f"n_refine={n:3d}  "
            f"naive: <M>={naive_res[n]['m1']:.3e} <M^2>={naive_res[n]['m2']:.3e} "
            f"chi={naive_res[n]['chi']:.3e} U4={naive_res[n]['u4']:.4f}  |  "
            f"eqarea: <M>={eqarea_res[n]['m1']:.3e} <M^2>={eqarea_res[n]['m2']:.3e} "
            f"chi={eqarea_res[n]['chi']:.3e} U4={eqarea_res[n]['u4']:.4f}",
            file=sys.stderr,
        )

    # ---- one small-multiples figure, all 6 magnetization observables ----
    ncols = 3
    nrows = int(np.ceil(len(METRICS) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.6 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for i, (key, err_key, ylabel, title, yscale, zero_ref) in enumerate(METRICS):
        ax = axes[i]
        naive_vals = np.array([naive_res[n][key] for n in ns], dtype=float)
        naive_errs = np.array([naive_res[n][err_key] for n in ns], dtype=float)
        eq_vals = np.array([eqarea_res[n][key] for n in ns], dtype=float)
        eq_errs = np.array([eqarea_res[n][err_key] for n in ns], dtype=float)
        if yscale == "log":
            ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=COLOR_NAIVE,
                        label="naive", ms=5, linewidth=1.6, capsize=3)
            ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=COLOR_EQAREA,
                        label="equal_area", ms=5, linewidth=1.6, capsize=3)
            ax.set_yscale("log")
        else:
            ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=COLOR_NAIVE,
                        label="naive", ms=5, linewidth=1.6, capsize=3)
            ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=COLOR_EQAREA,
                        label="equal_area", ms=5, linewidth=1.6, capsize=3)
        if zero_ref:
            ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
        ax.set_xscale("log")
        ax.set_title(title, fontsize=10.5)
        style_axes(ax)
        ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
        ax.set_ylabel(ylabel)
        if i == 0:
            ax.legend(fontsize=8, frameon=False)
    for j in range(len(METRICS), len(axes)):
        axes[j].axis("off")
    fig.suptitle(
        f"Magnetization observables: naive vs equal_area (n_refine $\\leq$ {max(ns)})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(out_dir, "mesh_compare_magnetization_vs_invL.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")

    # ---- ratio summary (naive/eqarea) for the strictly-positive moments ----
    ratio_metrics = [m for m in METRICS if m[0] in ("m2", "m4", "chi")]
    fig, axes = plt.subplots(1, len(ratio_metrics), figsize=(4.6 * len(ratio_metrics), 4.2), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for i, (key, err_key, ylabel, title, _yscale, _zero_ref) in enumerate(ratio_metrics):
        ax = axes[i]
        naive_vals = np.array([naive_res[n][key] for n in ns], dtype=float)
        naive_errs = np.array([naive_res[n][err_key] for n in ns], dtype=float)
        eq_vals = np.array([eqarea_res[n][key] for n in ns], dtype=float)
        eq_errs = np.array([eqarea_res[n][err_key] for n in ns], dtype=float)
        ratio = naive_vals / eq_vals
        ratio_err = np.abs(ratio) * np.sqrt((naive_errs / naive_vals) ** 2 + (eq_errs / eq_vals) ** 2)
        ax.errorbar(inv_L, ratio, yerr=ratio_err, fmt="o-", color=COLORS[2],
                    ms=5, linewidth=1.6, capsize=3)
        ax.axhline(1.0, color="gray", linewidth=1, linestyle="--", alpha=0.7)
        ax.set_xscale("log")
        ax.set_xlabel(r"$1/L$")
        ax.set_ylabel(f"{ylabel} ratio (naive/eqarea)")
        ax.set_title(title, fontsize=10.5)
        style_axes(ax)
    fig.suptitle("Mesh-mode ratio, magnetization observables", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = os.path.join(out_dir, "mesh_compare_magnetization_ratio_vs_invL.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
