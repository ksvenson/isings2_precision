#!/usr/bin/env python3
"""Build the M_l[m,m'] = <S_lm S*_lm'> spherical-symmetry matrices OFFLINE
from bin/full_corr_test's MC-averaged, orbit-symmetrized all-pairs
correlator (*_full_corr_<seed>.dat) -- rather than accumulating them inside
the simulator per raw configuration, which is what ising_s2_crit.cc's
--l_max path did. Built 2026-08-26 per explicit user direction: "we will
build the matrices from the MC averaged two point functions. Building them
in the simulator was a mistake since the two point functions were not
smoothed by high stats yet."

Math: S_lm(c) = sum_i w_i Y*_lm(x_i) s_i(c) (same convention as
analyze_saved_configs.py's compute_Ml_full). By linearity of expectation,
<S_lm S*_lm'> = sum_i sum_j w_i w_j Y*_lm(x_i) Y_lm'(x_j) <s_i s_j>, i.e.
M_l = A^H @ corr @ A where A[i,k] = w_i * conj(Y_{l,m}(x_i)) and corr is
full_corr_test.cc's already-converged, orbit-symmetrized <s_i s_j> matrix
-- no raw per-configuration spins needed at all here, only the smoothed
correlator + site positions/weights it already wrote to disk.

Each *_full_corr_<seed>.dat file is one independent-seed shard (same
"each shard is its own jackknife unit" pattern as
fit_real_space_ladder.py). This script turns a set of shards at the same
(n_refine, coupling_rule) point into exactly the `blocks` data structure
scripts/symmetry_test.py's analyze_blocks() already expects (one block per
shard: blk['n'] = n_meas, blk['sum'][(l,m,mp)] = n_meas * M_l[m,mp] for
m<=mp) -- so R_l/chi2/kappa{2,3,4}_r and their delete-one-shard jackknife
errors are computed by reusing that module's existing, already-verified
code, not reimplemented here.

Usage:
  scripts/build_ylm_matrix_from_corr.py <full_corr_dat_file> [more files...]
    --l_max 8 [--min-blocks-factor 5]
"""
import argparse
import os
import struct
import sys

import numpy as np
from scipy.special import sph_harm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symmetry_test import analyze_blocks, full_pairs_for_l  # noqa: E402


def load_full_corr(path):
    """Returns (n_sites, n_meas, x, y, z, w, corr) -- corr is the full
    symmetric n_sites x n_sites matrix, mirrored from the upper-triangle
    storage full_corr_test.cc writes."""
    with open(path, "rb") as f:
        n = struct.unpack("<i", f.read(4))[0]
        n_meas = struct.unpack("<q", f.read(8))[0]
        n2 = struct.unpack("<i", f.read(4))[0]
        if n2 != n:
            raise ValueError(f"{path}: header n_sites={n} but WritePositions n_sites={n2}")
        pos = np.frombuffer(f.read(8 * 3 * n), dtype="<f8").reshape(n, 3)
        w = np.frombuffer(f.read(8 * n), dtype="<f8").copy()
        tri = np.frombuffer(f.read(), dtype="<f8")

    expected = n * (n + 1) // 2
    if len(tri) != expected:
        raise ValueError(f"{path}: got {len(tri)} triangle entries, expected {expected}")

    corr = np.zeros((n, n))
    idx = 0
    for i in range(n):
        k = n - i
        corr[i, i:] = tri[idx : idx + k]
        idx += k
    corr = corr + corr.T - np.diag(np.diag(corr))
    return n, n_meas, pos[:, 0].copy(), pos[:, 1].copy(), pos[:, 2].copy(), w, corr


def build_ylm_matrix(x, y, z, l_max):
    """Y[i,k] = Y_{l,m}(site i), matching QfeLatticeS2::CalcYlm's convention
    (boost::math::spherical_harmonic(l, m, polar, azimuth)) -- verbatim
    copy of analyze_saved_configs.py's helper of the same name (kept in
    sync by hand, not shared code, same convention as the rest of this
    campaign's standalone diagnostics)."""
    theta = np.arccos(np.clip(z, -1, 1))  # polar
    phi = np.arctan2(y, x)  # azimuthal
    lm_list = [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]
    Y = np.zeros((len(x), len(lm_list)), dtype=complex)
    for k, (l, m) in enumerate(lm_list):
        Y[:, k] = sph_harm(m, l, phi, theta)
    return lm_list, Y


def compute_Ml_from_corr(corr, w, lm_list, Y):
    """M_l = A^H @ corr @ A per harmonic level l, A[i,k]=w_i*conj(Y_lm(i)).
    Returns dict l -> (2l+1)x(2l+1) complex Hermitian matrix, indexed
    [m+l, mp+l]."""
    A = w[:, None] * np.conj(Y)  # (n_sites, n_lm)
    M_full = A.conj().T @ corr @ A  # (n_lm, n_lm), one big matmul for all l at once
    l_max = max(l for l, m in lm_list)
    out = {}
    for l in range(l_max + 1):
        idxs = [k for k, (ll, m) in enumerate(lm_list) if ll == l]
        out[l] = M_full[np.ix_(idxs, idxs)]
    return out


def shard_to_block(path, l_max):
    n, n_meas, x, y, z, w, corr = load_full_corr(path)
    lm_list, Y = build_ylm_matrix(x, y, z, l_max)
    Ml = compute_Ml_from_corr(corr, w, lm_list, Y)
    block_sum = {}
    for l in range(l_max + 1):
        for mm, mp in full_pairs_for_l(l):
            block_sum[(l, mm, mp)] = n_meas * Ml[l][mm + l, mp + l]
    return dict(n=n_meas, sum=block_sum), n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="*_full_corr_<seed>.dat shard file(s)")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--min-blocks-factor", type=int, default=5)
    args = ap.parse_args()

    blocks = []
    n_sites_seen = None
    for p in args.paths:
        blk, n_sites = shard_to_block(p, args.l_max)
        if n_sites_seen is None:
            n_sites_seen = n_sites
        elif n_sites != n_sites_seen:
            raise ValueError(f"{p}: n_sites={n_sites} mismatches earlier shard's {n_sites_seen}")
        blocks.append(blk)
        print(f"# {p}: n_meas={blk['n']}, n_sites={n_sites}", file=sys.stderr)

    n_blocks, total_n, results = analyze_blocks(args.l_max, blocks, args.min_blocks_factor)

    print(f"# {len(args.paths)} shard(s), n_sites={n_sites_seen}")
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


if __name__ == "__main__":
    main()
