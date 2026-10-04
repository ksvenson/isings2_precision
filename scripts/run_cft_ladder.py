#!/usr/bin/env python3
"""Run cft_symmetry_test.py's Owen_section_D diagnostics (Delta_s(a),
delta_l(a)) across a mesh-refinement ladder, for one or two --mesh_mode
values, and plot/tabulate vs 1/n_refine. Mirrors compare_push7_modes.py's
structure/pooling convention for kappa2_r.

Usage:
  scripts/run_cft_ladder.py <run_dir> --ladder 2,3,4,6,8,12,16,24,32,48,64,96,128 \
      --modes naive,equal_area --ls 2,3,4,5,6,7,8
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from cft_symmetry_test import analyze

MODE_SUFFIX = {"naive": "", "equal_area": "_eqarea", "equal_rp": "_eqrp"}


def pooled_paths(run_dir, k, mode):
    suffix = MODE_SUFFIX.get(mode, "")
    pat = os.path.join(run_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_jackblocks_*.dat")
    return sorted(glob.glob(pat))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--ladder", default="2,3,4,6,8,12,16,24,32,48,64,96,128")
    ap.add_argument("--modes", default="naive,equal_area")
    ap.add_argument("--ls", default="2,3,4,5,6,7,8", help="l values to plot delta_l for")
    ap.add_argument("--out-prefix", default=None)
    ap.add_argument("--expect-shards", type=int, default=None,
                     help="warn if a ladder point doesn't have this many shard files")
    args = ap.parse_args()

    ladder = [int(x) for x in args.ladder.split(",")]
    modes = args.modes.split(",")
    ls = [int(x) for x in args.ls.split(",")]
    out_prefix = args.out_prefix or os.path.join(args.run_dir, "plots", "cft_owen_secD")

    data = {mode: {"k": [], "invk": [], "delta_s": [], "delta_s_err": [],
                    "delta_l": {l: [] for l in ls}, "delta_l_err": {l: [] for l in ls}}
            for mode in modes}

    for k in ladder:
        for mode in modes:
            paths = pooled_paths(args.run_dir, k, mode)
            if args.expect_shards is not None and len(paths) != args.expect_shards:
                print(f"WARNING: k={k} mode={mode} found {len(paths)}/{args.expect_shards} shard files",
                      file=sys.stderr)
            if not paths:
                continue
            r = analyze(paths)
            d = data[mode]
            d["k"].append(k)
            d["invk"].append(1.0 / k)
            d["delta_s"].append(r["delta_s"])
            d["delta_s_err"].append(r["delta_s_err"])
            for l in ls:
                if l in r["delta_l"]:
                    d["delta_l"][l].append(r["delta_l"][l])
                    d["delta_l_err"][l].append(r["delta_l_err"][l])
                else:
                    d["delta_l"][l].append(np.nan)
                    d["delta_l_err"][l].append(np.nan)

    # table
    print(f"{'n_refine':>9} {'mode':>12} {'Delta_s(a)':>16} {'err':>12} "
          + " ".join(f"delta_{l:<10}" for l in ls))
    for mode in modes:
        d = data[mode]
        for i, k in enumerate(d["k"]):
            dl_str = " ".join(f"{d['delta_l'][l][i]:+.3e}" for l in ls)
            print(f"{k:9d} {mode:>12} {d['delta_s'][i]:16.6f} {d['delta_s_err'][i]:12.6f} {dl_str}")

    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    colors = {"naive": "tab:red", "equal_area": "tab:blue"}

    # Delta_s(a) vs 1/n_refine
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for mode in modes:
        d = data[mode]
        if not d["k"]:
            continue
        ax.errorbar(d["invk"], d["delta_s"], yerr=d["delta_s_err"], marker="o", ms=4,
                    label=mode, color=colors.get(mode), capsize=2, lw=1)
    ax.axhline(0.125, color="gray", lw=1, ls="--", label="exact Delta_sigma=1/8")
    ax.set_xlabel("1/n_refine")
    ax.set_ylabel("Delta_s(a) = 2 F_1/(F_1+F_0)")
    ax.set_title("Owen_section_D: lattice scaling-dimension estimator vs refinement")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_prefix + "_delta_s.png", dpi=150)
    print(f"wrote {out_prefix}_delta_s.png")

    # delta_l(a) vs 1/n_refine, one panel per l
    n_l = len(ls)
    ncols = 4
    nrows = (n_l + ncols - 1) // ncols
    fig2, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)
    for idx, l in enumerate(ls):
        ax = axes[idx // ncols][idx % ncols]
        for mode in modes:
            d = data[mode]
            if not d["k"]:
                continue
            ax.errorbar(d["invk"], d["delta_l"][l], yerr=d["delta_l_err"][l], marker="o", ms=3,
                        label=mode, color=colors.get(mode), capsize=2, lw=1)
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.set_title(f"l={l}")
        ax.set_xlabel("1/n_refine")
        ax.set_ylabel("delta_l(a)")
        if idx == 0:
            ax.legend(fontsize=8)
    for idx in range(n_l, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    fig2.suptitle("Owen_section_D: dimension-independent symmetry-breaking delta_l(a) vs 1/n_refine")
    fig2.tight_layout()
    fig2.savefig(out_prefix + "_delta_l.png", dpi=150)
    print(f"wrote {out_prefix}_delta_l.png")


if __name__ == "__main__":
    main()
