#!/usr/bin/env python3
"""REAL R_l/kappa2_r spherical-symmetry check from bin/lean_harmonic_stats
data -- requires the 2026-08-27+ format (off-diagonal M_l[m,m'] terms,
see src/lean_harmonic_stats.cc's header and CLAUDE.md's "how many
operators" table). Converts each shard's blocks into the {n, sum:
{(l,m,mp): complex}} structure scripts/symmetry_test.py's analyze_blocks()
expects (analyze_lean_harmonic_stats.to_analyze_blocks_format) and reuses
that already-verified R_l/kappa2_r/kappa3_r/kappa4_r + delete-one-shard
jackknife code directly -- same pipeline save_configs_ladder_2026-08-26's
offline M_l analysis (build_ylm_matrix_from_configs.py) uses, just fed
from online-accumulated lean data instead of raw saved configs.

Supersedes symmetry_check_lean.py's diagonal-only chi2 proxy for any data
written after this driver's off-diagonal upgrade -- that script remains
useful only for pre-upgrade files (e.g. the lean_ladder_overnight_
2026-08-26 ladder), which don't have off-diagonal data to compute the
real kappa2_r from.

Usage:
  scripts/symmetry_check_lean_kappa2.py <dir_with_*_lean_*.dat files> [--l_max 8] [-o out_prefix]
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
from symmetry_test import analyze_blocks  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def point_kappa2(paths, l_max, min_blocks_factor=5):
    header, lm_list, blocks, offdiag_list = combine(paths)
    ab_blocks = to_analyze_blocks_format(lm_list, offdiag_list, blocks)
    n_blocks, total_n, results = analyze_blocks(l_max, ab_blocks, min_blocks_factor)
    return n_blocks, total_n, results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out_prefix", default=None)
    args = ap.parse_args()
    out_prefix = args.out_prefix or os.path.join(args.data_dir, "plots", "kappa2r_from_lean")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)(?:_eqarea|_eqrp)?_lean", f).group(1))
        by_n.setdefault(n, []).append(f)

    n_refines = sorted(by_n)
    kappa2r = np.full((len(n_refines), args.l_max + 1), np.nan)
    kappa2r_err = np.full((len(n_refines), args.l_max + 1), np.nan)
    R_l = np.full((len(n_refines), args.l_max + 1), np.nan)

    print(f"{'n_refine':>9}  " + "  ".join(f"l={l}(kappa2_r)" for l in range(1, args.l_max + 1)))
    for i, n in enumerate(n_refines):
        n_blocks, total_n, results = point_kappa2(sorted(by_n[n]), args.l_max)
        row = []
        for r in results:
            l = r["l"]
            if l == 0:
                continue
            kappa2r[i, l] = r["kappa2_r"]
            kappa2r_err[i, l] = r["kappa2_r_err"]
            R_l[i, l] = r["R_l"]
            row.append(f"{r['kappa2_r']:10.3e}")
        print(f"{n:>9}  " + "  ".join(row))

    inv_L = 1.0 / np.array(n_refines, dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for l in range(1, args.l_max + 1):
        color = COLORS[(l - 1) % len(COLORS)]
        axes[0].errorbar(inv_L, kappa2r[:, l], yerr=kappa2r_err[:, l], fmt="o-", ms=4, color=color, label=f"l={l}")
        axes[1].plot(inv_L, R_l[:, l], "o-", ms=4, color=color, label=f"l={l}")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$1/n_{\rm refine}$")
    axes[0].set_ylabel(r"$\kappa_{2,r}$ (basis-independent)")
    axes[0].set_title(r"Real $\kappa_{2,r}$ from lean off-diagonal data")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$1/n_{\rm refine}$")
    axes[1].set_ylabel(r"$R_l$ (basis-dependent, off-diag power)")
    axes[1].set_title(r"$R_l$")
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    out = out_prefix + ".png"
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
