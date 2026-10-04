#!/usr/bin/env python3
"""Print cft_symmetry_test.py's Delta_l_pair(a) table (every adjacent
(l-1,l) recursion-inversion estimate, not just l=0,1) across a mesh
refinement ladder for one --mesh_mode, pooling all shards per ladder
point. Companion to run_cft_ladder.py (which plots delta_s/delta_l but
not the full delta_l_pair table) -- built for push8's equal_rp analysis,
generalizes to any push directory / mode.

Usage:
  scripts/print_delta_l_pair_ladder.py <run_dir> --ladder 4,6,8,12,16,24,32 \
      --mode equal_rp --l_max 8
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from cft_symmetry_test import analyze

MODE_SUFFIX = {"naive": "", "equal_area": "_eqarea", "equal_rp": "_eqrp"}


def pooled_paths(run_dir, k, mode):
    suffix = MODE_SUFFIX[mode]
    pat = os.path.join(run_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}",
                        "*_ylm_2pt_full_jackblocks_*.dat")
    return sorted(glob.glob(pat))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--ladder", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODE_SUFFIX))
    ap.add_argument("--expect-shards", type=int, default=None)
    args = ap.parse_args()

    ladder = [int(x) for x in args.ladder.split(",")]

    per_k = {}
    l_max = None
    for k in ladder:
        paths = pooled_paths(args.run_dir, k, args.mode)
        if args.expect_shards is not None and len(paths) != args.expect_shards:
            print(f"WARNING: k={k} found {len(paths)}/{args.expect_shards} shard files", file=sys.stderr)
        if not paths:
            print(f"WARNING: k={k} no jackblocks files found, skipping", file=sys.stderr)
            continue
        r = analyze(paths)
        per_k[k] = r
        l_max = r["l_max"] if l_max is None else min(l_max, r["l_max"])

    ladder_ok = [k for k in ladder if k in per_k]
    header = "l  | " + " | ".join(f"n_refine={k:<4d}" for k in ladder_ok)
    print(header)
    print("-" * len(header))
    for l in range(1, l_max + 1):
        row = [f"{l:2d} |"]
        for k in ladder_ok:
            r = per_k[k]
            val = r["delta_l_pair"][l]
            err = r["delta_l_pair_err"][l]
            row.append(f" {val:7.4f}({err*1e4:4.0f}e-4)")
        print(" ".join(row))

    print()
    print("# n_blocks per ladder point (jack_block_size pooled across shards):")
    for k in ladder_ok:
        print(f"#   n_refine={k}: n_blocks={per_k[k]['n_blocks']}")


if __name__ == "__main__":
    main()
