#!/usr/bin/env python3
"""One-off: naive vs equal_area kappa2_r comparison table (ratio + sigma)
for the lean_mesh_compare_100shards_2026-08-27 campaign, restricted to the
n_refine=4..16 range where both mesh modes have full 100/100 shards.
"""
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symmetry_check_lean_kappa2 import point_kappa2

BASE = "campaign_runs/lean_mesh_compare_100shards_2026-08-27/analysis_upto16"
L_MAX = 8


def load(data_dir):
    files = glob.glob(os.path.join(data_dir, "q5k*_lean_*.dat"))
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)(?:_eqarea|_eqrp)?_lean", f).group(1))
        by_n.setdefault(n, []).append(f)
    out = {}
    for n, paths in sorted(by_n.items()):
        _, _, results = point_kappa2(paths, L_MAX)
        out[n] = {l: (r["kappa2_r"], r["kappa2_r_err"]) for l, r in enumerate(results)}
    return out


naive = load(os.path.join(BASE, "naive"))
eqarea = load(os.path.join(BASE, "eqarea"))

print(f"{'n':>4} {'l':>3} {'naive':>12} {'eqarea':>12} {'ratio(n/e)':>11} {'sigma(diff)':>12}")
for n in sorted(naive):
    for l in range(3, L_MAX + 1):
        vn, en = naive[n][l]
        ve, ee = eqarea[n][l]
        ratio = vn / ve if ve else float("nan")
        sigma = (vn - ve) / (en ** 2 + ee ** 2) ** 0.5
        print(f"{n:>4} {l:>3} {vn:>12.4e} {ve:>12.4e} {ratio:>11.2f} {sigma:>12.2f}")
