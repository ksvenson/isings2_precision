#!/usr/bin/env python3
"""naive-vs-equal_area comparison of the Owen_section_D recursion-relation
diagnostics -- delta_l(a) (dimension-independent symmetry-breaking measure)
and Delta_l_pair(a) (per-l scaling-dimension estimate from inverting the
single recursion step F_{l-1}->F_l, of which Delta_s(a) at l=1 is the
special case) -- computed from the SAME lean_harmonic_stats data used for
the R_l/kappa2_r comparison (plot_mesh_compare_full.py), not the raw
save_configs jackblocks cft_symmetry_test.py's own analyze() expects.

Reuses cft_symmetry_test.py's compute_F_l/delta_l_formula/
delta_l_pair_formula/jackknife_err verbatim; only the block-loading step
is swapped for analyze_lean_harmonic_stats.combine +
to_analyze_blocks_format (lean files), mirroring the same swap
symmetry_check_lean_kappa2.py already does for R_l/kappa2_r.

Usage:
  scripts/plot_delta_pairs_mesh_compare.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out_dir]
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
from analyze_lean_harmonic_stats import combine, to_analyze_blocks_format  # noqa: E402
from cft_symmetry_test import (  # noqa: E402
    compute_F_l, delta_l_formula, delta_l_pair_formula, jackknife_err,
)

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


def point_recursion(paths, l_max):
    """Mirrors symmetry_check_lean_kappa2.point_kappa2, but for the
    recursion-relation quantities (delta_l, delta_l_pair) instead of
    R_l/kappa_r."""
    _, lm_list, blocks, offdiag_list = combine(paths)
    ab_blocks = to_analyze_blocks_format(lm_list, offdiag_list, blocks)
    F_full, F_jk = compute_F_l(l_max, ab_blocks)

    delta_l = {}
    delta_l_err = {}
    for l in range(2, l_max + 1):
        dl_full = delta_l_formula(F_full[0], F_full[1], F_full[l - 1], F_full[l], l)
        dl_jk = delta_l_formula(F_jk[:, 0], F_jk[:, 1], F_jk[:, l - 1], F_jk[:, l], l)
        delta_l[l] = dl_full
        delta_l_err[l] = jackknife_err(dl_full, dl_jk)

    delta_l_pair = {}
    delta_l_pair_err = {}
    for l in range(1, l_max + 1):
        dp_full = delta_l_pair_formula(F_full[l - 1], F_full[l], l)
        dp_jk = delta_l_pair_formula(F_jk[:, l - 1], F_jk[:, l], l)
        delta_l_pair[l] = dp_full
        delta_l_pair_err[l] = jackknife_err(dp_full, dp_jk)

    return dict(delta_l=delta_l, delta_l_err=delta_l_err,
                delta_l_pair=delta_l_pair, delta_l_pair_err=delta_l_pair_err)


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
        naive_res[n] = point_recursion(sorted(naive_by_n[n]), args.l_max)
        eqarea_res[n] = point_recursion(sorted(eqarea_by_n[n]), args.l_max)
        print(
            f"n_refine={n:3d}  "
            + "  ".join(
                f"l={l}:naive={naive_res[n]['delta_l_pair'][l]:.4f},"
                f"eq={eqarea_res[n]['delta_l_pair'][l]:.4f}"
                for l in range(1, args.l_max + 1)
            ),
            file=sys.stderr,
        )

    # ---- Delta_l_pair(a): one panel per l=1..l_max, naive vs eqarea vs 1/L ----
    ls_pair = list(range(1, args.l_max + 1))
    ncols = 3
    nrows = int(np.ceil(len(ls_pair) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for i, l in enumerate(ls_pair):
        ax = axes[i]
        naive_vals = np.array([naive_res[n]["delta_l_pair"][l] for n in ns])
        naive_errs = np.array([naive_res[n]["delta_l_pair_err"][l] for n in ns])
        eq_vals = np.array([eqarea_res[n]["delta_l_pair"][l] for n in ns])
        eq_errs = np.array([eqarea_res[n]["delta_l_pair_err"][l] for n in ns])
        ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=COLOR_NAIVE,
                    label="naive", ms=5, linewidth=1.6, capsize=3)
        ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=COLOR_EQAREA,
                    label="equal_area", ms=5, linewidth=1.6, capsize=3)
        ax.axhline(0.125, color=COLORS[2], linestyle="--", linewidth=1.2, label=r"exact $1/8$")
        ax.set_xscale("log")
        ax.set_title(f"l-1={l-1} -> l={l}", fontsize=11)
        style_axes(ax)
        ax.set_xlabel(r"$1/L$")
        if i % ncols == 0:
            ax.set_ylabel(r"$\Delta_{l,\rm pair}(a)$")
        if i == 0:
            ax.legend(fontsize=8, frameon=False)
    for j in range(len(ls_pair), len(axes)):
        axes[j].axis("off")
    fig.suptitle(
        r"$\Delta_{l,\rm pair}(a)$ (single-step recursion inversion, all $l$ pairs): "
        f"naive vs equal_area (n_refine $\\leq$ {max(ns)})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(out_dir, "mesh_compare_delta_l_pair_vs_invL.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")

    # ---- delta_l(a): dimension-independent test, one panel per l=2..l_max ----
    ls_dl = list(range(2, args.l_max + 1))
    nrows = int(np.ceil(len(ls_dl) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for i, l in enumerate(ls_dl):
        ax = axes[i]
        naive_vals = np.array([naive_res[n]["delta_l"][l] for n in ns])
        naive_errs = np.array([naive_res[n]["delta_l_err"][l] for n in ns])
        eq_vals = np.array([eqarea_res[n]["delta_l"][l] for n in ns])
        eq_errs = np.array([eqarea_res[n]["delta_l_err"][l] for n in ns])
        ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=COLOR_NAIVE,
                    label="naive", ms=5, linewidth=1.6, capsize=3)
        ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=COLOR_EQAREA,
                    label="equal_area", ms=5, linewidth=1.6, capsize=3)
        ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
        ax.set_xscale("log")
        ax.set_title(f"l = {l}", fontsize=11)
        style_axes(ax)
        ax.set_xlabel(r"$1/L$")
        if i % ncols == 0:
            ax.set_ylabel(r"$\delta_l(a)$")
        if i == 0:
            ax.legend(fontsize=8, frameon=False)
    for j in range(len(ls_dl), len(axes)):
        axes[j].axis("off")
    fig.suptitle(
        r"$\delta_l(a)$ (dimension-independent conformal-symmetry test): "
        f"naive vs equal_area (n_refine $\\leq$ {max(ns)})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(out_dir, "mesh_compare_delta_l_vs_invL.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
