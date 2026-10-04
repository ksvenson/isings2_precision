#!/usr/bin/env python3
"""push7 kappa2_r(naive)/kappa2_r(equal_area) ratio vs n_refine, one line per l.

Reuses the same pooling logic as compare_push7_modes.py (32-shard pool per
(n_refine, mesh_mode) via symmetry_test.analyze_multi) but plots the ratio
directly instead of the two curves separately.

Usage:
  scripts/plot_push7_ratio.py <push7_dir> [--ls 1,2,3,4,5,6,7,8]
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
# below this, both modes are still ~bare icosahedron / aliasing-contaminated
# at l_max=8 -- excluded from the ratio plot (see journal.md 2026-08-25)
SAFE_MIN_NREFINE = 4


def pooled_paths(push7_dir, k, mode):
    suffix = "" if mode == "naive" else "_eqarea"
    pat = os.path.join(push7_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_jackblocks_*.dat")
    return sorted(glob.glob(pat))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("push7_dir")
    ap.add_argument("--ls", default="3,4,5,6,7,8")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ls = [int(x) for x in args.ls.split(",")]
    out = args.out or os.path.join(args.push7_dir, "plots", "naive_over_eqarea_ratio.png")

    data = {mode: {l: {} for l in ls} for mode in ("naive", "equal_area")}

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
                data[mode][l][k] = (r["kappa2_r"], r["kappa2_r_err"])

    # l=1,2 are excluded by default (see --ls): kappa2_r is consistent with
    # zero for both mesh modes there (icosahedral group has no invariant
    # below l=6), so their ratio is noise divided by noise, not signal --
    # keeping them in a "cleaner" plot just adds decade-spanning error bars
    # that swamp the real l=3..8 trend.

    # ordered sequential ramp (light -> dark), since l is an ordered index,
    # not an unordered category -- avoids a rainbow/viridis categorical read.
    # Floor kept at #6baed6 (not lighter): anything paler fails contrast on
    # a white surface.
    RAMP = ["#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b", "#041c38"]
    plt.rcParams.update({
        "font.size": 11,
        "axes.edgecolor": "#888888",
        "axes.labelcolor": "#333333",
        "xtick.color": "#555555",
        "ytick.color": "#555555",
    })

    fig, ax = plt.subplots(figsize=(8.2, 5.5))
    end_points = []  # (l, x_last, y_last, color) for label collision resolution
    for i, l in enumerate(ls):
        ks, ratios, ratio_errs = [], [], []
        for k in LADDER:
            if k < SAFE_MIN_NREFINE:
                continue
            if k not in data["naive"][l] or k not in data["equal_area"][l]:
                continue
            kn, kne = data["naive"][l][k]
            ke, kee = data["equal_area"][l][k]
            if ke == 0:
                continue
            ratio = kn / ke
            # propagate error on a ratio of two independent noisy quantities
            rel_err = np.sqrt((kne / kn) ** 2 + (kee / ke) ** 2) if kn != 0 else np.nan
            ks.append(k)
            ratios.append(ratio)
            ratio_errs.append(abs(ratio) * rel_err)
        if not ks:
            continue
        color = RAMP[i % len(RAMP)]
        ax.errorbar(ks, ratios, yerr=ratio_errs, marker="o", ms=4.5, lw=1.8,
                     capsize=0, elinewidth=0.8, alpha=0.95, color=color,
                     ecolor=color, zorder=3)
        end_points.append((l, ks[-1], ratios[-1], color))

    # direct end-of-line labels, in place of a legend box -- nudge apart in
    # log-y when several lines land within a decade of each other at the
    # right edge (l=6,7,8 all converge near ratio~1)
    end_points.sort(key=lambda p: p[2])
    placed_log_y = []
    min_gap = 0.09  # in log10(y) units
    x_right = max(x for _, x, _, _ in end_points)
    for l, x, y, color in end_points:
        ly = np.log10(y)
        if placed_log_y and ly - placed_log_y[-1] < min_gap:
            ly = placed_log_y[-1] + min_gap
        placed_log_y.append(ly)
        label_y = 10 ** ly
        # thin leader line from the actual last point to the de-collided
        # label row when they diverge, so the label still reads unambiguously
        if abs(label_y - y) / y > 1e-3:
            ax.plot([x, x_right * 1.06], [y, label_y], color=color, lw=0.6,
                    alpha=0.5, zorder=2, clip_on=False)
        ax.annotate(f"l={l}", xy=(x_right * 1.06, label_y), xytext=(4, 0),
                    textcoords="offset points", xycoords="data",
                    va="center", fontsize=9.5, color=color, fontweight="bold",
                    annotation_clip=False)

    ax.axhline(1.0, color="#999999", lw=1, ls="--", zorder=1)
    ax.text(SAFE_MIN_NREFINE, 1.08, "no difference", fontsize=8.5,
            color="#777777", va="bottom")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("n_refine (mesh resolution)")
    ax.set_ylabel(r"$\kappa_2^r$(naive) / $\kappa_2^r$(equal_area)")
    fig.suptitle("push7: SO(3)-symmetry-breaking ratio, naive vs. equal_area mesh",
                 fontsize=13, fontweight="bold", y=0.99)
    ax.set_title("matched stats, l_max=8  ·  larger = equal_area restores symmetry faster",
                 fontsize=9.5, color="#666666", fontweight="normal", pad=10)
    ax.set_xlim(SAFE_MIN_NREFINE * 0.85, 128 * 2.4)
    ax.grid(True, which="major", axis="both", alpha=0.25, lw=0.6)
    ax.grid(True, which="minor", axis="both", alpha=0.1, lw=0.4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
