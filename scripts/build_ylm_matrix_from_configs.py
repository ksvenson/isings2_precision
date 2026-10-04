#!/usr/bin/env python3
"""Build the M_l[m,m'] = <S_lm S*_lm'> spherical-symmetry matrices OFFLINE
from bin/save_configs's raw saved (positions, bit-packed configs) shards --
the storage-efficient alternative to build_ylm_matrix_from_corr.py's
full_corr_test.cc path (2026-08-26 "use 1" decision -- see journal.md).

Both paths compute the exact same estimator: by linearity of expectation,
<S_lm S*_lm'> = sum_ij w_i w_j Y*_lm(x_i) Y_lm'(x_j) <s_i s_j> is identical
whether you (a) form the smoothed <s_i s_j> correlator first and contract
it with the Y_lm weights afterward (full_corr_test.cc), or (b) project
each saved raw configuration onto the Y_lm basis (S = spins @ A) and only
average the *outer product* S^H @ S over configs at the end -- both are
"average over configs, then done", neither does anything nonlinear before
that average. The only difference is computational: (a) needs
O(n_sites^2) storage for the intermediate correlator; (b) needs only
O(n_meas * n_sites) bits for the raw configs (typically orders of
magnitude smaller at large n_refine -- e.g. n_refine=128, n_meas=10000:
~205 MB here vs. ~215 GB for the dense matrix) and does the Y_lm
contraction via two O(n_meas * n_sites * n_lm) matmuls, never materializing
an n_sites x n_sites matrix at all.

Reuses scripts/analyze_saved_configs.py's load_positions/load_configs/
build_ylm_matrix/compute_Ml_full directly (not duplicated) since that
module already implements exactly the per-shard M_l computation needed
here; this script's only job is to turn one-or-more shards into the
`blocks` structure scripts/symmetry_test.py's analyze_blocks() expects
(one block per shard: blk['n']=n_configs, blk['sum'][(l,m,mp)] =
n_configs * M_l[m,mp] for m<=mp) so the existing, already-verified
R_l/chi2/kappa{2,3,4}_r + delete-one-shard-jackknife code is reused
verbatim, not reimplemented -- same pattern as
build_ylm_matrix_from_corr.py.

**Also computes F_l/Delta_s/Delta_l_pair** (added 2026-08-26, per user
request "do a projection onto ylms and try the delta = 2F0/F1 method"),
reusing cft_symmetry_test.py's compute_F_l/delta_s_formula/
delta_l_pair_formula directly against the SAME `blocks` list already built
for the R_l/kappa2_r table above -- no second matmul pass, since F_l is
just Trace(M_l)/(2*pi)/(2l+1) and everything needed is already in `blocks`.

**UPDATE, 2026-08-26 later same day**: the original ~0.33 plateau this was
expected to reproduce (per the 2026-08-25 "BREAKTHROUGH" pivot) turned out
to be a real bug in `compute_F_l` -- it was missing a division by `(2l+1)`
(see that function's docstring in cft_symmetry_test.py and journal.md's
"RESOLVED: the north-pole-vs-Trace(M_l) disagreement" entry for the full
derivation). With the fix, a full-statistics rerun at n_refine=4 (96
shards, 960,000 configs from save_configs_ladder_2026-08-26) gives
`Delta_s(a) = 0.13955 +/- 0.00024` -- matching the real-space method's
independently-derived n_refine=4 point (~0.139) to 3 decimal places. The
"l-space is fundamentally unreliable due to UV contamination" conclusion
this docstring used to cite is no longer safe to assume; re-derive/check
before trusting either framing on new data.

Usage:
  scripts/build_ylm_matrix_from_configs.py \\
    <positions1.dat> <configs1.bin> [<positions2.dat> <configs2.bin> ...] \\
    --l_max 8 [--min-blocks-factor 5] [--max-configs-per-shard N]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_saved_configs import build_ylm_matrix, compute_Ml_full, load_configs, load_positions  # noqa: E402
from cft_symmetry_test import compute_F_l, delta_l_pair_formula, delta_s_formula, jackknife_err  # noqa: E402
from symmetry_test import analyze_blocks, full_pairs_for_l  # noqa: E402


def shard_to_block(pos_path, configs_path, l_max, max_configs):
    x, y, z, w = load_positions(pos_path)
    n_sites = len(x)
    spins = load_configs(configs_path, n_sites, max_configs)
    n_meas = spins.shape[0]

    lm_list, Y = build_ylm_matrix(x, y, z, l_max)
    M_full, _ = compute_Ml_full(spins, w, lm_list, Y)

    l_index = {}
    for k, (l, m) in enumerate(lm_list):
        l_index.setdefault(l, []).append(k)

    block_sum = {}
    for l in range(l_max + 1):
        idx = l_index[l]
        Ml = M_full[np.ix_(idx, idx)]  # (2l+1)x(2l+1), indexed [m+l, mp+l]
        for mm, mp in full_pairs_for_l(l):
            block_sum[(l, mm, mp)] = n_meas * Ml[mm + l, mp + l]
    return dict(n=n_meas, sum=block_sum), n_sites, n_meas


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="alternating <positions.dat> <configs.bin> pairs, one pair per shard")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--min-blocks-factor", type=int, default=5)
    ap.add_argument("--max-configs-per-shard", type=int, default=None)
    args = ap.parse_args()

    if len(args.files) % 2 != 0:
        raise SystemExit("expected an even number of files (positions.dat configs.bin pairs)")
    shard_pairs = list(zip(args.files[0::2], args.files[1::2]))

    blocks = []
    n_sites_seen = None
    for pos_path, configs_path in shard_pairs:
        blk, n_sites, n_meas = shard_to_block(pos_path, configs_path, args.l_max, args.max_configs_per_shard)
        if n_sites_seen is None:
            n_sites_seen = n_sites
        elif n_sites != n_sites_seen:
            raise ValueError(f"{pos_path}: n_sites={n_sites} mismatches earlier shard's {n_sites_seen}")
        blocks.append(blk)
        print(f"# {pos_path} + {configs_path}: n_meas={n_meas}, n_sites={n_sites}", file=sys.stderr)

    n_blocks, total_n, results = analyze_blocks(args.l_max, blocks, args.min_blocks_factor)

    print(f"# {len(shard_pairs)} shard(s), n_sites={n_sites_seen}")
    print(f"# l_max={args.l_max} n_blocks(shards)={n_blocks} total_samples={total_n}")
    print(f"# R_l/chi2 are basis-dependent (mesh's fixed quantization axis); "
          f"kappa{{2,3,4}}_r are basis-independent spectral cumulants of M_l's "
          f"eigenvalues (0 iff M_l is proportional to the identity)")
    print(f"{'l':>3} {'dim':>4} {'R_l':>14} {'R_l_err':>14} {'chi2':>10} {'dof':>4} "
          f"{'chi2/dof':>9} {'c_hat':>14} {'kappa2_r':>14} {'err':>14} "
          f"{'kappa3_r':>14} {'err':>14} {'kappa4_r':>14} {'err':>14}")
    for r in results:
        chi2_str = f"{r['chi2']:10.3f}" if r["chi2"] is not None else f"{'--':>10}"
        rchi2_str = (
            f"{r['chi2'] / r['dof']:9.3f}" if r["chi2"] is not None and r["dof"] > 0 else f"{'--':>9}"
        )
        c_hat_str = f"{r['c_hat']:14.6e}" if r["c_hat"] is not None else f"{'--':>14}"
        k3_str = f"{r['kappa3_r']:14.6e}" if r["kappa3_r"] is not None else f"{'--':>14}"
        k3e_str = f"{r['kappa3_r_err']:14.6e}" if r["kappa3_r_err"] is not None else f"{'--':>14}"
        k4_str = f"{r['kappa4_r']:14.6e}" if r["kappa4_r"] is not None else f"{'--':>14}"
        k4e_str = f"{r['kappa4_r_err']:14.6e}" if r["kappa4_r_err"] is not None else f"{'--':>14}"
        print(f"{r['l']:>3} {r['dim']:>4} {r['R_l']:14.6e} {r['R_l_err']:14.6e} "
              f"{chi2_str} {r['dof']:>4} {rchi2_str} {c_hat_str} "
              f"{r['kappa2_r']:14.6e} {r['kappa2_r_err']:14.6e} {k3_str} {k3e_str} {k4_str} {k4e_str}")

    # F_l / Delta_s / Delta_l_pair -- l-space Delta estimators, reusing
    # cft_symmetry_test.py's formulas against the same `blocks` (see
    # module docstring: expected to reproduce the known ~0.33 l-space
    # anomaly, not confirm 1/8 -- included as a direct requested
    # cross-check, not a validation attempt).
    F_full, F_jk = compute_F_l(args.l_max, blocks)
    F_err = jackknife_err(F_full, F_jk)
    ds_full = delta_s_formula(F_full[0], F_full[1])
    ds_jk = delta_s_formula(F_jk[:, 0], F_jk[:, 1])
    ds_err = jackknife_err(ds_full, ds_jk)

    print(f"\n# Delta_s(a) = 2*F_1/(F_1+F_0) [exact value 1/8 = 0.125; F_l now divided "
          f"by (2l+1) -- see cft_symmetry_test.py's compute_F_l docstring, "
          f"2026-08-26 fix. Matched the real-space method to 3 decimal places at "
          f"n_refine=4 after this fix -- no longer expected to be anomalous]:")
    print(f"#   {ds_full:.8f} +/- {ds_err:.8f}")
    print(f"# F_l(a):")
    for l in range(args.l_max + 1):
        print(f"#   l={l:2d}  F_l={F_full[l]:14.6e} +/- {F_err[l]:14.6e}")
    print(f"# Delta_l_pair(a) [single-step recursion inversion at each l -- should all "
          f"agree with each other AND with Delta_s=Delta_1_pair if this is one CFT "
          f"primary; l=6 is expected to be icosahedral-invariant-contaminated "
          f"regardless of continuum limit -- see journal.md]:")
    print(f"{'l':>3} {'Delta_l_pair':>16} {'err':>16}")
    for l in range(1, args.l_max + 1):
        dp_full = delta_l_pair_formula(F_full[l - 1], F_full[l], l)
        dp_jk = delta_l_pair_formula(F_jk[:, l - 1], F_jk[:, l], l)
        dp_err = jackknife_err(dp_full, dp_jk)
        print(f"{l:3d} {dp_full:16.6e} {dp_err:16.6e}")


if __name__ == "__main__":
    main()
