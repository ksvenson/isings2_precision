#!/usr/bin/env python3
"""Combine per-shard M_l block pickles (compute_ml_block_shard.py output)
into the same R_l/kappa2_r/F_l/Delta_l_pair table
build_ylm_matrix_from_configs.py prints -- cheap (analyze_blocks does no
matmul, just sums + a delete-one-shard jackknife), so this runs fine on
the login node even though the per-shard matmuls that produced the
pickles needed to be split across separate SGE tasks.

Usage:
  scripts/reduce_ml_blocks.py <block1.pkl> [<block2.pkl> ...]
"""
import argparse
import pickle
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cft_symmetry_test import compute_F_l, delta_l_pair_formula, delta_s_formula, jackknife_err  # noqa: E402
from symmetry_test import analyze_blocks  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pickles", nargs="+")
    ap.add_argument("--min-blocks-factor", type=int, default=5)
    args = ap.parse_args()

    blocks = []
    l_max = n_sites_seen = None
    for path in args.pickles:
        with open(path, "rb") as f:
            d = pickle.load(f)
        if l_max is None:
            l_max, n_sites_seen = d["l_max"], d["n_sites"]
        elif d["l_max"] != l_max:
            raise ValueError(f"{path}: l_max={d['l_max']} != {l_max}")
        elif d["n_sites"] != n_sites_seen:
            raise ValueError(f"{path}: n_sites={d['n_sites']} != {n_sites_seen}")
        blocks.append(d["block"])

    n_blocks, total_n, results = analyze_blocks(l_max, blocks, args.min_blocks_factor)

    print(f"# {len(args.pickles)} shard(s), n_sites={n_sites_seen}")
    print(f"# l_max={l_max} n_blocks(shards)={n_blocks} total_samples={total_n}")
    print(f"# R_l/chi2 are basis-dependent (mesh's fixed quantization axis); "
          f"kappa{{2,3,4}}_r are basis-independent spectral cumulants of M_l's "
          f"eigenvalues (0 iff M_l is proportional to the identity)")
    print(f"{'l':>3} {'dim':>4} {'R_l':>14} {'R_l_err':>14} {'chi2':>10} {'dof':>4} "
          f"{'chi2/dof':>9} {'c_hat':>14} {'kappa2_r':>14} {'err':>14} "
          f"{'kappa3_r':>14} {'err':>14} {'kappa4_r':>14} {'err':>14}")
    for r in results:
        chi2_str = f"{r['chi2']:10.3f}" if r["chi2"] is not None else f"{'--':>10}"
        chi2dof_str = f"{r['chi2']/r['dof']:9.3f}" if r["chi2"] is not None and r["dof"] > 0 else f"{'--':>9}"
        k3 = f"{r['kappa3_r']:14.6e}" if r["kappa3_r"] is not None else f"{'--':>14}"
        k3e = f"{r['kappa3_r_err']:14.6e}" if r["kappa3_r_err"] is not None else f"{'--':>14}"
        k4 = f"{r['kappa4_r']:14.6e}" if r["kappa4_r"] is not None else f"{'--':>14}"
        k4e = f"{r['kappa4_r_err']:14.6e}" if r["kappa4_r_err"] is not None else f"{'--':>14}"
        c_hat_str = f"{r['c_hat']:14.6e}" if r["c_hat"] is not None else f"{'--':>14}"
        print(f"{r['l']:3d} {r['dim']:4d} {r['R_l']:14.6e} {r['R_l_err']:14.6e} "
              f"{chi2_str} {r['dof']:4d} {chi2dof_str} {c_hat_str} "
              f"{r['kappa2_r']:14.6e} {r['kappa2_r_err']:14.6e} {k3} {k3e} {k4} {k4e}")

    # F_l / Delta_s / Delta_l_pair, same as build_ylm_matrix_from_configs.py
    F_full, F_jk = compute_F_l(l_max, blocks)
    F_err = jackknife_err(F_full, F_jk)
    ds_full = delta_s_formula(F_full[0], F_full[1])
    ds_jk = delta_s_formula(F_jk[:, 0], F_jk[:, 1])
    ds_err = jackknife_err(ds_full, ds_jk)
    print(f"# Delta_s(a) = 2*F_1/(F_1+F_0) [exact value 1/8 = 0.125]: {ds_full:.6f} +/- {ds_err:.6f}")
    print(f"{'l':>3} {'Delta_l_pair':>16} {'err':>16}")
    for l in range(1, l_max + 1):
        dp_full = delta_l_pair_formula(F_full[l - 1], F_full[l], l)
        dp_jk = delta_l_pair_formula(F_jk[:, l - 1], F_jk[:, l], l)
        dp_err = jackknife_err(dp_full, dp_jk)
        print(f"{l:3d} {dp_full:16.6e} {dp_err:16.6e}")


if __name__ == "__main__":
    main()
