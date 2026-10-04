#!/usr/bin/env python3
"""Direct naive-vs-equal_area comparison for push7's dual mesh_mode ladder.

Pools all 32 shards per (n_refine, mesh_mode) combo via
symmetry_test.analyze_multi, extracts kappa2_r(l) (basis-independent
SO(3)-breaking cumulant) for l=1..8, and reports/plots naive vs
equal_area side by side plus their ratio.

Usage:
  scripts/compare_push7_modes.py <push7_dir> [--ls 1,2,3,4,5,6,7,8]
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
from symmetry_test import analyze_multi

LADDER = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]


def pooled_paths(push7_dir, k, mode):
    suffix = "" if mode == "naive" else "_eqarea"
    pat = os.path.join(push7_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_jackblocks_*.dat")
    paths = sorted(glob.glob(pat))
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("push7_dir")
    ap.add_argument("--ls", default="1,2,3,4,5,6,7,8")
    ap.add_argument("--out-prefix", default=None)
    args = ap.parse_args()

    ls = [int(x) for x in args.ls.split(",")]
    out_prefix = args.out_prefix or os.path.join(args.push7_dir, "plots", "naive_vs_eqarea")

    data = {mode: {l: {"invn": [], "k": [], "k2": [], "k2e": []} for l in ls}
            for mode in ("naive", "equal_area")}

    for k in LADDER:
        for mode in ("naive", "equal_area"):
            paths = pooled_paths(args.push7_dir, k, mode)
            if len(paths) != 32:
                print(f"WARNING: k={k} mode={mode} found {len(paths)}/32 shard files", file=sys.stderr)
            if not paths:
                continue
            l_max, _, _, _, results = analyze_multi(paths)
            by_l = {r["l"]: r for r in results}
            for l in ls:
                if l > l_max:
                    continue
                r = by_l[l]
                d = data[mode][l]
                d["invn"].append(1.0 / k)
                d["k"].append(k)
                d["k2"].append(r["kappa2_r"])
                d["k2e"].append(r["kappa2_r_err"])

    # print table
    print(f"{'l':>3} {'n_refine':>9} {'naive kappa2_r':>18} {'eqarea kappa2_r':>18} {'ratio naive/eq':>15} {'sigma of diff':>14}")
    for l in ls:
        dn, de = data["naive"][l], data["equal_area"][l]
        for i, k in enumerate(dn["k"]):
            if k not in de["k"]:
                continue
            j = de["k"].index(k)
            kn, kne = dn["k2"][i], dn["k2e"][i]
            ke, kee = de["k2"][j], de["k2e"][j]
            diff = kn - ke
            sig = diff / np.sqrt(kne**2 + kee**2) if (kne > 0 or kee > 0) else float("nan")
            ratio = kn / ke if ke != 0 else float("inf")
            print(f"{l:3d} {k:9d} {kn:10.4e}+-{kne:.1e} {ke:10.4e}+-{kee:.1e} {ratio:15.2f} {sig:14.2f}")

    # plot: one panel per l (up to 8), kappa2_r vs 1/n_refine, naive vs eqarea overlaid
    n_l = len(ls)
    ncols = 4
    nrows = (n_l + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)
    for idx, l in enumerate(ls):
        ax = axes[idx // ncols][idx % ncols]
        for mode, color in (("naive", "tab:red"), ("equal_area", "tab:blue")):
            d = data[mode][l]
            if not d["invn"]:
                continue
            ax.errorbar(d["invn"], d["k2"], yerr=d["k2e"], marker="o", ms=3,
                        label=mode, color=color, capsize=2, lw=1)
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.set_title(f"l={l}")
        ax.set_xlabel("1/n_refine")
        ax.set_ylabel("kappa2_r")
        if idx == 0:
            ax.legend(fontsize=8)
    for idx in range(n_l, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    fig.suptitle("push7: naive vs equal_area kappa2_r(l) vs 1/n_refine (l_max=8, matched stats)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    fig.savefig(out_prefix + "_overlay.png", dpi=150)
    print(f"wrote {out_prefix}_overlay.png")

    # log-log version
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)
    for idx, l in enumerate(ls):
        ax = axes2[idx // ncols][idx % ncols]
        for mode, color in (("naive", "tab:red"), ("equal_area", "tab:blue")):
            d = data[mode][l]
            if not d["invn"]:
                continue
            invn = np.array(d["invn"])
            k2 = np.abs(np.array(d["k2"]))
            k2e = np.array(d["k2e"])
            ax.errorbar(invn, k2, yerr=k2e, marker="o", ms=3,
                        label=mode, color=color, capsize=2, lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"l={l}")
        ax.set_xlabel("1/n_refine")
        ax.set_ylabel("|kappa2_r|")
        if idx == 0:
            ax.legend(fontsize=8)
    for idx in range(n_l, nrows * ncols):
        axes2[idx // ncols][idx % ncols].axis("off")
    fig2.suptitle("push7: naive vs equal_area |kappa2_r(l)| vs 1/n_refine (log-log)")
    fig2.tight_layout()
    fig2.savefig(out_prefix + "_overlay_loglog.png", dpi=150)
    print(f"wrote {out_prefix}_overlay_loglog.png")


if __name__ == "__main__":
    main()
