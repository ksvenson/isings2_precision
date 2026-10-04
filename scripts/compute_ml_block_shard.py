#!/usr/bin/env python3
"""Compute ONE shard's M_l jackknife block (build_ylm_matrix_from_configs.
shard_to_block) and pickle it to disk, instead of looping over all shards
serially in one process. Split out 2026-08-27 because
build_ylm_matrix_from_configs.py's serial per-point loop hit SGE's 8h
h_rt wall-clock limit at n_refine=64/128 (job 7322113 tasks 5/6, exit 137
after exactly 28801s -- NOT an OOM kill, maxvmem was only 6-18GB, see
journal.md/CLAUDE.md 2026-08-27) before finishing all 96 shards. This
script does the same per-shard work but as an independently
schedulable/parallelizable unit -- pairs with reduce_ml_blocks.py, same
shard-then-reduce pattern as every other production push in this
campaign (production_shard_task.sh, save_configs_shard_task.sh, etc).

Usage:
  scripts/compute_ml_block_shard.py <positions.dat> <configs.bin> --l_max 8 -o <block.pkl>
"""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_ylm_matrix_from_configs import shard_to_block  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("positions")
    ap.add_argument("configs")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--max-configs", type=int, default=None)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    blk, n_sites, n_meas = shard_to_block(args.positions, args.configs, args.l_max, args.max_configs)
    with open(args.out, "wb") as f:
        pickle.dump(dict(block=blk, n_sites=n_sites, n_meas=n_meas, l_max=args.l_max), f)
    print(f"wrote {args.out}: n_sites={n_sites} n_meas={n_meas}", file=sys.stderr)


if __name__ == "__main__":
    main()
