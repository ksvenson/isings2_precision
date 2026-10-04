#!/usr/bin/env python3
"""Same output as plot_mesh_compare_full.py (all four symmetry metrics
R_l/kappa2_r/kappa3_r/kappa4_r, the scalar Delta_s/U4 figure, and the
ratio-summary figure) but restricted to a --max_n cutoff, so a full set
of plots can be produced from whatever ladder points are already cheap
to compute without waiting on the slowest (largest n_refine) point.
Reuses plot_mesh_compare_full.py's own collect()/point_kappa2()/
analyze_scalar() plotting code verbatim (imported, not reimplemented) --
this file only adds the --max_n filter on top.

Usage:
  scripts/plot_mesh_compare_partial.py <naive_dir> <eqarea_dir> --max_n 112 [--l_max 8] [-o out_dir]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_mesh_compare_full as full  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("naive_dir")
    ap.add_argument("eqarea_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--max_n", type=int, required=True)
    ap.add_argument("-o", "--out_dir", default=None)
    args = ap.parse_args()

    naive_by_n = full.collect(args.naive_dir)
    eqarea_by_n = full.collect(args.eqarea_dir)
    ns = sorted(n for n in (set(naive_by_n) & set(eqarea_by_n)) if n <= args.max_n)
    inv_L = 1.0 / np.array(ns, dtype=float)

    out_dir = args.out_dir or os.path.join(os.path.dirname(args.naive_dir.rstrip("/")), "plots_partial")
    os.makedirs(out_dir, exist_ok=True)

    ls = list(range(3, args.l_max + 1))
    naive_rows, eqarea_rows = {}, {}
    for n in ns:
        _, _, rn = full.point_kappa2(sorted(naive_by_n[n]), args.l_max)
        _, _, re_ = full.point_kappa2(sorted(eqarea_by_n[n]), args.l_max)
        naive_rows[n] = {r["l"]: r for r in rn}
        eqarea_rows[n] = {r["l"]: r for r in re_}
        print(f"n_refine={n:3d}  " + "  ".join(
            f"l={l}:naive_k2={naive_rows[n][l]['kappa2_r']:.2e},eq_k2={eqarea_rows[n][l]['kappa2_r']:.2e}"
            for l in ls), file=sys.stderr)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for metric, metric_err, ylabel in full.SYMMETRY_METRICS:
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
            ax.errorbar(inv_L, naive_vals, yerr=naive_errs, fmt="o-", color=full.COLOR_NAIVE,
                        label="naive", ms=5, linewidth=1.6, capsize=3)
            ax.errorbar(inv_L, eq_vals, yerr=eq_errs, fmt="s-", color=full.COLOR_EQAREA,
                        label="equal_area", ms=5, linewidth=1.6, capsize=3)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(f"l = {l}", fontsize=11)
            full.style_axes(ax)
            ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
            if i % ncols == 0:
                ax.set_ylabel(ylabel)
            if i == 0:
                ax.legend(fontsize=9, frameon=False)
        for j in range(len(ls), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"PARTIAL ({ylabel}): naive vs equal_area vs $1/L$ "
                     f"(n_refine $\\leq$ {max(ns)}, full ladder goes to 512/... check journal)", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        out = os.path.join(out_dir, f"mesh_compare_{metric}_vs_invL_partial.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"wrote {out}")

    naive_scalar, eqarea_scalar = [], []
    for n in ns:
        naive_scalar.append(full.analyze_scalar(sorted(naive_by_n[n]), l_max=args.l_max))
        eqarea_scalar.append(full.analyze_scalar(sorted(eqarea_by_n[n]), l_max=args.l_max))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ds_n = np.array([r["delta_s"] for r in naive_scalar])
    ds_n_err = np.array([r["delta_s_err"] for r in naive_scalar])
    ds_e = np.array([r["delta_s"] for r in eqarea_scalar])
    ds_e_err = np.array([r["delta_s_err"] for r in eqarea_scalar])
    ax.errorbar(inv_L, ds_n, yerr=ds_n_err, fmt="o-", color=full.COLOR_NAIVE, label="naive",
                ms=5, linewidth=1.6, capsize=3)
    ax.errorbar(inv_L, ds_e, yerr=ds_e_err, fmt="s-", color=full.COLOR_EQAREA, label="equal_area",
                ms=5, linewidth=1.6, capsize=3)
    ax.axhline(0.125, color=full.COLORS[2], linestyle="--", linewidth=1.2, label=r"exact $1/8$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
    ax.set_ylabel(r"$\Delta_s(a)$")
    ax.set_title(r"$\Delta_s(a) = 2F_1/(F_1+F_0)$")
    full.style_axes(ax)
    ax.legend(fontsize=9, frameon=False)

    ax = axes[1]
    u4_n = np.array([r["u4"] for r in naive_scalar])
    u4_n_err = np.array([r["u4_err"] for r in naive_scalar])
    u4_e = np.array([r["u4"] for r in eqarea_scalar])
    u4_e_err = np.array([r["u4_err"] for r in eqarea_scalar])
    ax.errorbar(inv_L, u4_n, yerr=u4_n_err, fmt="o-", color=full.COLOR_NAIVE, label="naive",
                ms=5, linewidth=1.6, capsize=3)
    ax.errorbar(inv_L, u4_e, yerr=u4_e_err, fmt="s-", color=full.COLOR_EQAREA, label="equal_area",
                ms=5, linewidth=1.6, capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
    ax.set_ylabel(r"Binder $U_4$")
    ax.set_title(r"Binder cumulant $U_4 = 1 - \langle M^4\rangle/(3\langle M^2\rangle^2)$")
    full.style_axes(ax)
    ax.legend(fontsize=9, frameon=False)

    fig.suptitle(f"PARTIAL scalar bulk/CFT-data metrics: naive vs equal_area (n_refine $\\leq$ {max(ns)})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = os.path.join(out_dir, "mesh_compare_scalar_vs_invL_partial.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
    axes = axes.ravel()
    l_colors = {l: full.COLORS[i % len(full.COLORS)] for i, l in enumerate(ls)}
    for mi, (metric, metric_err, ylabel) in enumerate(full.SYMMETRY_METRICS):
        ax = axes[mi]
        for l in ls:
            naive_vals = np.array([naive_rows[n][l][metric] for n in ns], dtype=float)
            naive_errs = np.array([naive_rows[n][l][metric_err] for n in ns], dtype=float)
            eq_vals = np.array([eqarea_rows[n][l][metric] for n in ns], dtype=float)
            eq_errs = np.array([eqarea_rows[n][l][metric_err] for n in ns], dtype=float)
            ratio = naive_vals / eq_vals
            ratio_err = np.abs(ratio) * np.sqrt((naive_errs / naive_vals) ** 2 + (eq_errs / eq_vals) ** 2)
            ax.errorbar(inv_L, ratio, yerr=ratio_err, fmt="o-", color=l_colors[l], label=f"l={l}",
                        ms=4, linewidth=1.3, capsize=2)
        ax.axhline(1.0, color="gray", linewidth=1, linestyle="--", alpha=0.7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$1/L$ (= $1/n_{\rm refine}$)")
        ax.set_ylabel(f"{ylabel} ratio (naive/eqarea)")
        full.style_axes(ax)
        if mi == 0:
            ax.legend(fontsize=8, frameon=False, ncol=2)
    fig.suptitle(f"PARTIAL mesh-mode ratio across all four symmetry-breaking metrics (n_refine $\\leq$ {max(ns)})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(out_dir, "mesh_compare_ratio_all_metrics_vs_invL_partial.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
