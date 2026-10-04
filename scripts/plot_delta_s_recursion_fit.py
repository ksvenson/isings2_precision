#!/usr/bin/env python3
"""Delta_{l,pair}(a) vs 1/L, one panel per harmonic level l=1..l_max, with
continuum extrapolations from the recursion relation overlaid -- naive vs
equal_area, built from the same lean mesh-compare data (and the same
compute_F_l/delta_l_pair_formula recursion machinery) as
plot_delta_pairs_mesh_compare.py's Delta_{l,pair}(a) panels.
Delta_s(a) is exactly the l=1 panel of this family (built from F_0, F_1
only) -- this script extends that single case to every l pair with fitted
extrapolations to a=0 added on each panel, rather than a new measurement.

Per-panel (per naive/eqarea, per l) two continuum-extrapolation models are
fit and both overlaid, reusing fit_lean_ladder_family's implementations:
  p2      Delta_inf + c/n^2        (fixed quadratic power)
  pfloat  Delta_inf + c/n^p        (p floated, 3 free params)
Per this campaign's established finding (journal.md's per-l family-fit
entry), l=1 alone prefers p~1 while l>=2 generally prefer p~2-2.7 -- so a
fixed p=2 is not assumed correct for every level, and the floated-power
fit is shown alongside it as a check on that assumption rather than in
place of it.

Usage:
  scripts/plot_delta_s_recursion_fit.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out_dir]
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
from plot_delta_pairs_mesh_compare import collect, point_recursion  # noqa: E402
from fit_lean_ladder_family import fit_single_power, fit_floated_power  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
COLOR_NAIVE = COLORS[0]
COLOR_EQAREA = COLORS[1]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)


def fit_pair(n_refines, deltas, errs):
    """p=2 fixed and floated-power fits, both against n_refine (not 1/L)."""
    p2 = fit_single_power(n_refines, deltas, errs, 2.0)
    pfloat, p_val = fit_floated_power(n_refines, deltas, errs)
    pfloat["p_value"] = p_val
    return p2, pfloat


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
    n_arr = np.array(ns, dtype=float)
    inv_L = 1.0 / n_arr

    out_dir = args.out_dir or os.path.join(os.path.dirname(args.naive_dir.rstrip("/")), "plots")
    os.makedirs(out_dir, exist_ok=True)

    naive_res, eqarea_res = {}, {}
    for n in ns:
        naive_res[n] = point_recursion(sorted(naive_by_n[n]), args.l_max)
        eqarea_res[n] = point_recursion(sorted(eqarea_by_n[n]), args.l_max)

    ls = list(range(1, args.l_max + 1))
    ncols = 3
    nrows = int(np.ceil(len(ls) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.8 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()

    n_fit = np.linspace(n_arr.min(), n_arr.max(), 400)
    x_fit = 1.0 / n_fit

    for i, l in enumerate(ls):
        ax = axes[i]
        naive_vals = np.array([naive_res[n]["delta_l_pair"][l] for n in ns])
        naive_errs = np.array([naive_res[n]["delta_l_pair_err"][l] for n in ns])
        eq_vals = np.array([eqarea_res[n]["delta_l_pair"][l] for n in ns])
        eq_errs = np.array([eqarea_res[n]["delta_l_pair_err"][l] for n in ns])

        naive_p2, naive_pfloat = fit_pair(n_arr, naive_vals, naive_errs)
        eq_p2, eq_pfloat = fit_pair(n_arr, eq_vals, eq_errs)

        print(
            f"l={l}: naive p2 Delta_inf={naive_p2['d_inf']:.5f}+/-{naive_p2['d_inf_err']:.5f} "
            f"chi2/dof={naive_p2['chi2']/naive_p2['dof']:.2f}  "
            f"pfloat p={naive_pfloat['p_value']:.2f} Delta_inf={naive_pfloat['d_inf']:.5f}"
            f"+/-{naive_pfloat['d_inf_err']:.5f}  |  "
            f"eqarea p2 Delta_inf={eq_p2['d_inf']:.5f}+/-{eq_p2['d_inf_err']:.5f}  "
            f"pfloat p={eq_pfloat['p_value']:.2f} Delta_inf={eq_pfloat['d_inf']:.5f}"
            f"+/-{eq_pfloat['d_inf_err']:.5f}",
            file=sys.stderr,
        )

        ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o", color=COLOR_NAIVE,
                    label="naive", ms=5, capsize=3)
        ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s", color=COLOR_EQAREA,
                    label="equal_area", ms=5, capsize=3)

        ax.plot(x_fit, naive_p2["model"](n_fit), "--", color=COLOR_NAIVE, linewidth=1.3,
                label=f"naive p=2: {naive_p2['d_inf']:.4f}$\\pm${naive_p2['d_inf_err']:.4f}")
        ax.plot(x_fit, eq_p2["model"](n_fit), "--", color=COLOR_EQAREA, linewidth=1.3,
                label=f"eqarea p=2: {eq_p2['d_inf']:.4f}$\\pm${eq_p2['d_inf_err']:.4f}")
        ax.plot(x_fit, naive_pfloat["model"](n_fit), ":", color=COLOR_NAIVE, linewidth=1.6,
                label=f"naive pfloat (p={naive_pfloat['p_value']:.2f}): {naive_pfloat['d_inf']:.4f}$\\pm${naive_pfloat['d_inf_err']:.4f}")
        ax.plot(x_fit, eq_pfloat["model"](n_fit), ":", color=COLOR_EQAREA, linewidth=1.6,
                label=f"eqarea pfloat (p={eq_pfloat['p_value']:.2f}): {eq_pfloat['d_inf']:.4f}$\\pm${eq_pfloat['d_inf_err']:.4f}")

        for fit, marker_color in ((naive_p2, COLOR_NAIVE), (eq_p2, COLOR_EQAREA),
                                   (naive_pfloat, COLOR_NAIVE), (eq_pfloat, COLOR_EQAREA)):
            ax.errorbar([0.0], [fit["d_inf"]], yerr=[fit["d_inf_err"]], fmt="D",
                        color=marker_color, ms=6, mfc="white", mew=1.3, zorder=5)
        ax.axhline(0.125, color=COLORS[2], linestyle=":", linewidth=1.0, alpha=0.7)

        style_axes(ax)
        ax.set_title(f"l-1={l-1} -> l={l}", fontsize=11)
        ax.set_xlabel(r"$1/L$")
        if i % ncols == 0:
            ax.set_ylabel(r"$\Delta_{l,\rm pair}(a)$")
        ax.legend(fontsize=6, frameon=False, loc="best")

    for j in range(len(ls), len(axes)):
        axes[j].axis("off")

    fig.suptitle(
        r"$\Delta_{l,\rm pair}(a)$ continuum extrapolation (p=2 and floated-power fits), "
        f"all $l$ pairs (n_refine {min(ns)}-{max(ns)})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(out_dir, "delta_l_pair_recursion_extrapolation.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
