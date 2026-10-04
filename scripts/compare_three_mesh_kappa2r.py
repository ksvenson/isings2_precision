#!/usr/bin/env python3
"""naive vs equal_area vs equal_rp kappa2_r comparison, ladder_2026-08-27's
100-shard run. Reuses symmetry_check_lean_kappa2.point_kappa2 verbatim (no
reimplementation of the cumulant math) against all three mesh-mode
directories under campaign_runs/lean_mesh_compare_100shards_2026-08-27/,
now that the third (equal_rp) production array (job 7358337) has landed
alongside the earlier naive/eqarea pair.

Usage: scripts/compare_three_mesh_kappa2r.py [--l_max 8]
"""
import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symmetry_check_lean_kappa2 import point_kappa2  # noqa: E402


def collect(d):
    by_n = {}
    for f in glob.glob(os.path.join(d, "*.dat")):
        m = re.search(r"q5k(\d+)_", os.path.basename(f))
        if m:
            by_n.setdefault(int(m.group(1)), []).append(f)
    return by_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument(
        "--base",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..",
            "campaign_runs", "lean_mesh_compare_100shards_2026-08-27",
        ),
    )
    args = ap.parse_args()

    naive_by_n = collect(os.path.join(args.base, "naive"))
    eqarea_by_n = collect(os.path.join(args.base, "eqarea"))
    eqrp_by_n = collect(os.path.join(args.base, "eqrp"))

    ns = sorted(set(naive_by_n) & set(eqarea_by_n) & set(eqrp_by_n))
    print(f"n_refine values compared: {ns}")
    header = f"{'n_refine':>9}  " + "  ".join(
        f"{lbl}_l{l}" for l in range(3, args.l_max + 1) for lbl in ("naive", "eqarea", "eqrp")
    )
    print(header)
    for n in ns:
        _, _, rn = point_kappa2(sorted(naive_by_n[n]), args.l_max)
        _, _, re_ = point_kappa2(sorted(eqarea_by_n[n]), args.l_max)
        _, _, rp_ = point_kappa2(sorted(eqrp_by_n[n]), args.l_max)
        row = [f"{n:9d}"]
        for l in range(3, args.l_max + 1):
            kn = rn[l]["kappa2_r"]
            ke = re_[l]["kappa2_r"]
            kp = rp_[l]["kappa2_r"]
            row.append(f"{kn:10.3e} {ke:10.3e} {kp:10.3e}")
        print("  ".join(row))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
