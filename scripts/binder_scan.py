#!/usr/bin/env python3
"""Light Binder-cumulant (U4) scan across a mesh-refinement ladder, to
check scale invariance independent of the harmonic-decomposition method
(cft_symmetry_test.py/symmetry_test.py) -- a genuinely critical point has
U4 = 1.5*(1 - <m^4>/(3*<m^2>^2)) roughly resolution-independent (the
classic Binder-crossing criticality signature); a residual mass scale
shows up as U4 drifting monotonically with n_refine instead of
plateauing.

Reads *_bulk_*.dat files directly (mag^2/mag^4 lines: mean, naive
(non-autocorr-corrected) error, n -- QfeMeasReal::WriteMeasurement's
format, include/statistics.h) and exactly reconstructs sum/sum2/n per
shard (QfeMeasReal::ReadMeasurement's own inverse formula), so shards can
be pooled exactly (sum over shards) rather than error-propagated from
already-averaged per-shard means. Error bars via leave-one-shard-out
jackknife (shard-level, not the finer jack_block_size blocks the ylm
files use -- "light" scan, no new file format needed).

Usage:
  scripts/binder_scan.py <run_dir> --ladder 2,3,4,6,8,12,16,24,32,48,64,96,128 \
      --mode naive
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

MODE_SUFFIX = {"naive": "", "equal_area": "_eqarea", "equal_rp": "_eqrp"}


def read_bulk_moment(path, name):
    """Returns (sum, sum2, n) for a named line ('mag^2' or 'mag^4')."""
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == name:
                mean, err, n = float(parts[1]), float(parts[2]), int(parts[3])
                sum_ = mean * n
                sum2 = (err * err * n + mean * mean) * n
                return sum_, sum2, n
    raise ValueError(f"{name} not found in {path}")


def pooled_paths(run_dir, k, mode):
    suffix = MODE_SUFFIX[mode]
    pat = os.path.join(run_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_bulk_*.dat")
    return sorted(glob.glob(pat))


def u4_from_moments(m2_sum, m2_n, m4_sum, m4_n):
    m2 = m2_sum / m2_n
    m4 = m4_sum / m4_n
    return 1.5 * (1.0 - m4 / (3.0 * m2 * m2))


def analyze_point(paths):
    m2 = []  # (sum, sum2, n) per shard
    m4 = []
    for p in paths:
        m2.append(read_bulk_moment(p, "mag^2"))
        m4.append(read_bulk_moment(p, "mag^4"))

    m2_sum_tot = sum(x[0] for x in m2)
    m2_n_tot = sum(x[2] for x in m2)
    m4_sum_tot = sum(x[0] for x in m4)
    m4_n_tot = sum(x[2] for x in m4)

    u4_full = u4_from_moments(m2_sum_tot, m2_n_tot, m4_sum_tot, m4_n_tot)

    n_shards = len(paths)
    u4_jk = np.empty(n_shards)
    for i in range(n_shards):
        m2_sum_loo = m2_sum_tot - m2[i][0]
        m2_n_loo = m2_n_tot - m2[i][2]
        m4_sum_loo = m4_sum_tot - m4[i][0]
        m4_n_loo = m4_n_tot - m4[i][2]
        u4_jk[i] = u4_from_moments(m2_sum_loo, m2_n_loo, m4_sum_loo, m4_n_loo)

    bar = np.mean(u4_jk)
    u4_err = np.sqrt((n_shards - 1) / n_shards * np.sum((u4_jk - bar) ** 2))
    return u4_full, u4_err, n_shards, m2_n_tot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--ladder", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODE_SUFFIX))
    ap.add_argument("--expect-shards", type=int, default=None)
    args = ap.parse_args()

    ladder = [int(x) for x in args.ladder.split(",")]

    print(f"{'n_refine':>9} {'U4':>12} {'err':>12} {'n_shards':>9} {'pooled_traj':>12}")
    for k in ladder:
        paths = pooled_paths(args.run_dir, k, args.mode)
        if args.expect_shards is not None and len(paths) != args.expect_shards:
            print(f"WARNING: k={k} found {len(paths)}/{args.expect_shards} shard files", file=sys.stderr)
        if not paths:
            print(f"WARNING: k={k} no bulk files found, skipping", file=sys.stderr)
            continue
        u4, err, n_shards, pooled_traj = analyze_point(paths)
        print(f"{k:9d} {u4:12.6f} {err:12.6f} {n_shards:9d} {pooled_traj:12d}")


if __name__ == "__main__":
    main()
