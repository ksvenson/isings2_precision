#!/usr/bin/env python3
"""Compute the MC expectation value and variance of the full two-point
function <s_i s_j>, entirely OFFLINE from bin/save_configs's raw saved
shards -- built 2026-08-26 per user direction: only want the MC
expectation of the two-point function and its variance, computed after
the simulations are done (not accumulated online during MC, and not
built from raw per-configuration Y_lm projections).

Per shard, the mean correlator is one fast matmul (corr_shard =
spins.T @ spins / n_meas -- BLAS-vectorized, much cheaper in practice
than the equivalent O(n_sites^2) work done as n_meas separate rank-1
updates interleaved inside the MC loop, which is what src/full_corr_test.cc
did and is not needed here). Shards are then combined via Welford's
online algorithm (mean and variance updated one shard at a time, matching
this campaign's "accumulated stats, mean/variance of step n updated to
step n+1" convention -- just applied across independent shards instead of
across configs, since each shard's own MC-time average already handles
in-chain autocorrelation and shards are independent MC chains): this
gives the standard error of the pooled mean directly, without needing to
hold every shard's matrix in memory simultaneously (peak memory is ~3
dense matrices: running mean, running M2, current shard -- not
n_shards x n_sites^2).

**Memory is unavoidably O(n_sites^2)** for the final dense result itself
(mean + variance matrices), independent of how/when it's computed -- this
script refuses to run past --max-n-sites (default 20000, ~5.3 GB for the
two matrices) without --force, since going further (e.g. n_refine=128,
n_sites~1.6e5, ~430 GB) is not something to do by accident. Options if
you need larger n_refine: reduce to distinct icosahedral-orbit-class
values only (~120x smaller, not implemented here), or use
build_ylm_matrix_from_configs.py, which never materializes the dense
matrix at all (only ever forms n_sites x n_lm and n_lm x n_lm objects).

Usage:
  scripts/compute_two_point_stats.py \\
    <positions1.dat> <configs1.bin> [<positions2.dat> <configs2.bin> ...] \\
    -o <out_prefix> [--max-n-sites 20000] [--force]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_saved_configs import load_configs, load_positions  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="alternating <positions.dat> <configs.bin> pairs, one pair per shard")
    ap.add_argument("-o", "--out-prefix", required=True, help="writes <prefix>_mean.npy and <prefix>_err.npy")
    ap.add_argument("--max-n-sites", type=int, default=20000,
                     help="refuse to materialize a dense matrix larger than this without --force")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if len(args.files) % 2 != 0:
        raise SystemExit("expected an even number of files (positions.dat configs.bin pairs)")
    shard_pairs = list(zip(args.files[0::2], args.files[1::2]))

    x, y, z, w = load_positions(shard_pairs[0][0])
    n_sites = len(x)
    if n_sites > args.max_n_sites and not args.force:
        bytes_needed = 2 * n_sites * n_sites * 8
        raise SystemExit(
            f"n_sites={n_sites} exceeds --max-n-sites={args.max_n_sites} "
            f"(mean+M2 matrices would need ~{bytes_needed/1e9:.1f} GB) -- pass --force to proceed anyway, "
            f"or use build_ylm_matrix_from_configs.py instead, which never forms the dense matrix"
        )

    mean = np.zeros((n_sites, n_sites))
    M2 = np.zeros((n_sites, n_sites))
    count = 0
    total_configs = 0

    for pos_path, configs_path in shard_pairs:
        xs, ys, zs, ws = load_positions(pos_path)
        if len(xs) != n_sites:
            raise ValueError(f"{pos_path}: n_sites={len(xs)} mismatches first shard's {n_sites}")
        spins = load_configs(configs_path, n_sites)
        n_meas = spins.shape[0]
        total_configs += n_meas

        corr_shard = (spins.T @ spins) / n_meas  # one BLAS matmul, this shard's own MC-time average

        # Welford's algorithm across shards (each shard = one independent sample)
        count += 1
        delta = corr_shard - mean
        mean += delta / count
        delta2 = corr_shard - mean
        M2 += delta * delta2

        print(f"# {pos_path}: n_meas={n_meas}, running shard count={count}", file=sys.stderr)

    if count < 2:
        raise SystemExit(f"need >=2 shards for a variance estimate, got {count}")

    variance = M2 / (count - 1)  # sample variance of the shard-mean estimator
    sem = np.sqrt(variance / count)  # standard error of the pooled mean

    np.save(f"{args.out_prefix}_mean.npy", mean)
    np.save(f"{args.out_prefix}_err.npy", sem)
    np.savez(f"{args.out_prefix}_positions.npz", x=x, y=y, z=z, w=w)

    print(f"# {count} shards, {total_configs} total configs, n_sites={n_sites}")
    print(f"# wrote {args.out_prefix}_mean.npy, {args.out_prefix}_err.npy, {args.out_prefix}_positions.npz")
    print(f"# diag(mean) range: [{np.diag(mean).min():.6f}, {np.diag(mean).max():.6f}] (must be exactly 1.0)")


if __name__ == "__main__":
    main()
