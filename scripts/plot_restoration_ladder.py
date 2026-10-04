#!/usr/bin/env python3
"""Plot kappa2_r (and R_l) vs 1/n_refine for a chosen set of harmonic
levels, reading the *_ylm_2pt_full_jackblocks_*.dat file(s) already on disk
for each ladder point. Multiple shard files for one ladder point are pooled
via symmetry_test.analyze_multi (comma-separated paths).

Usage:
  scripts/plot_restoration_ladder.py <out_png> [--ls 1,2,3,8,12] \\
      <n_refine>:<jackblocks_path>[,<jackblocks_path>...] [...]

Example (single-file points, 2026-08-22):
  .venv_plot/bin/python scripts/plot_restoration_ladder.py \\
    campaign_runs/push_2026-08-22/plots/symmetry_signals_vs_invL.png \\
    2:campaign_runs/push_2026-08-22/q5k2/q5k2_ylm_2pt_full_jackblocks_00001092.dat \\
    4:campaign_runs/push_2026-08-22/q5k4/q5k4_ylm_2pt_full_jackblocks_00001092.dat

Example (pooled multi-shard points, 2026-08-24):
  .venv_plot/bin/python scripts/plot_restoration_ladder.py \\
    campaign_runs/analysis_lmax12_2026-08-24/plots/lmax8_ladder.png --ls 1,2,3,8 \\
    4:"$(echo campaign_runs/production_2026-08-23/q5k4/shard_*/q5k4/*_jackblocks_*.dat | tr ' ' ',')" \\
    8:"$(echo campaign_runs/production_2026-08-23/q5k8/shard_*/q5k8/*_jackblocks_*.dat | tr ' ' ',')"
"""
import sys
import os
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from symmetry_test import analyze, analyze_multi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_png")
    ap.add_argument("--ls", default="1,2,3", help="comma-separated harmonic levels to plot (default 1,2,3)")
    ap.add_argument("--log-log", action="store_true", help="log-scale both axes")
    ap.add_argument("points", nargs="+", help="<n_refine>:<path>[,<path>...]")
    args = ap.parse_args()

    ls = [int(x) for x in args.ls.split(",")]

    points = []
    for arg in args.points:
        n_refine_str, path_str = arg.split(":", 1)
        paths = path_str.split(",")
        points.append((int(n_refine_str), paths))
    points.sort()

    data = {l: {"invn": [], "k2": [], "k2e": [], "R": [], "Re": []} for l in ls}

    l_max_seen = None
    for n_refine, paths in points:
        if len(paths) == 1:
            l_max, _, _, _, results = analyze(paths[0])
        else:
            l_max, _, _, _, results = analyze_multi(paths)
        l_max_seen = l_max
        by_l = {r["l"]: r for r in results}
        for l in ls:
            r = by_l[l]
            data[l]["invn"].append(1.0 / n_refine)
            data[l]["k2"].append(r["kappa2_r"])
            data[l]["k2e"].append(r["kappa2_r_err"])
            data[l]["R"].append(r["R_l"])
            data[l]["Re"].append(r["R_l_err"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    cmap = plt.get_cmap("tab10")
    colors = {l: cmap(i % 10) for i, l in enumerate(ls)}
    for l in ls:
        d = data[l]
        axes[0].errorbar(d["invn"], d["k2"], yerr=d["k2e"], marker="o",
                          label=f"l={l}", color=colors[l], capsize=3)
        axes[1].errorbar(d["invn"], d["R"], yerr=d["Re"], marker="o",
                          label=f"l={l}", color=colors[l], capsize=3)
    axes[0].set_title("kappa2_r (basis-independent)")
    axes[1].set_title("R_l (basis-dependent)")
    for ax in axes:
        ax.set_xlabel("1 / n_refine")
        if args.log_log:
            ax.set_xscale("log")
            ax.set_yscale("log")
        else:
            ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.legend(fontsize=8)
    fig.suptitle(f"Symmetry-breaking signals vs mesh refinement (q=5, exact_sinh, l_max={l_max_seen})")
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out_png), exist_ok=True)
    fig.savefig(args.out_png, dpi=150)
    print(f"wrote {args.out_png}")


if __name__ == "__main__":
    main()
