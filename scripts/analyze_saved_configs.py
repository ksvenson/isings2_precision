#!/usr/bin/env python3
"""Post-process raw configs saved by bin/save_configs: compute the exact
S_lm/M_l[m,m']/F_l harmonic decomposition (same math as ising_s2_crit.cc,
done here in Python for full analysis flexibility) directly from saved
site positions/weights and packed spin configurations -- and, unlike the
production C++ pipeline, can also compute a UV-pair-excluded M_l (exclude
site pairs within an angular cutoff from the double sum defining M_l,
using a precomputed near-pair list -- something the factorized
S_lm=sum_i w_i Y*_lm(x_i) s_i trick cannot do on its own, since pair
exclusion is inherently bilinear). Built 2026-08-25 per user direction
("save the real-space correlators and do the spherical symmetry tests by
projecting afterwards") -- see IsingS2_precision/journal.md.

Usage:
  scripts/analyze_saved_configs.py <positions.dat> <configs.bin> [--l_max 8]
    [--uv-cut-theta 0.3] [--max-configs N]
"""
import argparse
import sys

import numpy as np
from scipy.special import sph_harm
from scipy.spatial import cKDTree


def load_positions(path):
    s, x, y, z, wt = [], [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            x.append(float(p[1])); y.append(float(p[2])); z.append(float(p[3])); wt.append(float(p[4]))
    return np.array(x), np.array(y), np.array(z), np.array(wt)


def load_configs(path, n_sites, max_configs=None):
    buf_size = (n_sites + 7) // 8
    raw = np.fromfile(path, dtype=np.uint8)
    n_configs = len(raw) // buf_size
    raw = raw[: n_configs * buf_size].reshape(n_configs, buf_size)
    bits = np.unpackbits(raw, axis=1, bitorder="little")[:, :n_sites]
    spins = np.where(bits == 1, -1.0, 1.0)
    if max_configs is not None:
        spins = spins[:max_configs]
    return spins


def build_ylm_matrix(x, y, z, l_max):
    """Returns (lm_list, Y) where Y[i, k] = Y_{l,m}(site i) for the k-th
    (l,m) pair in lm_list (m ranges over the full -l..l here, unlike the
    C++ driver's m>=0-only storage, since we're not trying to match its
    on-disk format -- just recomputing everything fresh)."""
    theta = np.arccos(np.clip(z, -1, 1))  # polar
    phi = np.arctan2(y, x)  # azimuthal
    lm_list = [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]
    Y = np.zeros((len(x), len(lm_list)), dtype=complex)
    for k, (l, m) in enumerate(lm_list):
        # scipy sph_harm(m, l, azimuth, polar) matches boost's Y_lm(l,m,polar,azimuth)
        Y[:, k] = sph_harm(m, l, phi, theta)
    return lm_list, Y


def compute_Ml_full(spins, w, lm_list, Y):
    """M_l[m,m'] = <S_lm S*_lm'>, S_lm = sum_i w_i Y*_lm(x_i) s_i -- the
    standard (unrestricted) harmonic projection, vectorized as one matmul.

    S = spins @ A is computed as two real matmuls (spins @ A.real,
    spins @ A.imag) rather than one real-times-complex matmul -- numpy
    upcasts the real `spins` operand to complex128 for a mixed-dtype
    matmul, which at production n_refine (e.g. n_refine=128,
    n_meas=10000: spins is 13 GB real) means allocating a ~26 GB complex
    copy of `spins` just to multiply it by an (n_sites, n_lm) matrix whose
    complex part is needed for a result 1000x smaller. Splitting into two
    real matmuls avoids ever materializing that copy."""
    A = w[:, None] * np.conj(Y)  # (n_sites, n_lm)
    S = spins @ A.real + 1j * (spins @ A.imag)  # (n_configs, n_lm)
    M = (S.conj().T @ S) / spins.shape[0]  # (n_lm, n_lm)
    return M, S


def compute_Ml_near(spins, w, x, y, z, lm_list, Y, theta_cut):
    """UV (short-distance) contribution to M_l[m,m']: same double sum as
    the full M_l, restricted to site pairs (including i=j) within
    theta_cut of each other. Uses a KD-tree on chord distance
    (chord = sqrt(2-2cos(theta_cut))) to find the near-pair list once,
    then computes each pair's config-averaged <s_i s_j> and forms the
    weighted Y*_lm(i)Y_lm'(j) sum only over those pairs."""
    pts = np.stack([x, y, z], axis=1)
    tree = cKDTree(pts)
    chord_cut = np.sqrt(max(0.0, 2 - 2 * np.cos(theta_cut)))
    pairs = tree.query_pairs(r=chord_cut, output_type="ndarray")  # i<j only
    print(f"  near-pair count (theta_cut={theta_cut}): {len(pairs)} (+ {len(x)} self-pairs)")

    n_configs = spins.shape[0]
    I, J = pairs[:, 0], pairs[:, 1]
    # <s_i s_j> for i<j pairs, chunked over pairs to bound memory
    corr_offdiag = np.empty(len(pairs))
    # keep chunk*n_configs*8 bytes bounded (~1.5GB) regardless of dataset size
    chunk = max(1000, int(1.5e8 / max(1, n_configs)))
    for start in range(0, len(pairs), chunk):
        end = min(start + chunk, len(pairs))
        corr_offdiag[start:end] = (spins[:, I[start:end]] * spins[:, J[start:end]]).mean(axis=0)
    # self-pairs (i=j): <s_i^2> = 1 always
    corr_diag = np.ones(len(x))

    n_lm = len(lm_list)
    M_near = np.zeros((n_lm, n_lm), dtype=complex)
    wI, wJ = w[I], w[J]
    for k in range(n_lm):
        YkI = np.conj(Y[I, k])
        # off-diagonal pairs contribute symmetrically (i,j) and (j,i)
        for kp in range(n_lm):
            YkpJ = Y[J, kp]
            term_off = np.sum(wI * wJ * YkI * YkpJ * corr_offdiag)
            YkJ = np.conj(Y[J, k])
            YkpI = Y[I, kp]
            term_off2 = np.sum(wJ * wI * YkJ * YkpI * corr_offdiag)
            term_self = np.sum(w * w * np.conj(Y[:, k]) * Y[:, kp] * corr_diag)
            M_near[k, kp] = term_off + term_off2 + term_self
    return M_near


def delta_l_pair_from_diag(F_l, l_max):
    out = {}
    for l in range(1, l_max + 1):
        rho = F_l[l] / F_l[l - 1]
        out[l] = (rho * (l + 1) - (l - 1)) / (1.0 + rho)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("positions")
    ap.add_argument("configs")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--uv-cut-theta", type=float, default=None)
    ap.add_argument("--max-configs", type=int, default=None)
    args = ap.parse_args()

    x, y, z, w = load_positions(args.positions)
    n_sites = len(x)
    spins = load_configs(args.configs, n_sites, args.max_configs)
    print(f"loaded {n_sites} sites, {spins.shape[0]} configs")

    lm_list, Y = build_ylm_matrix(x, y, z, args.l_max)
    M_full, S = compute_Ml_full(spins, w, lm_list, Y)

    # diagonal trace per l (sum over m of M[l,m; l,m]) -> F_l, same
    # convention as cft_symmetry_test.py (F_l = Trace(M_l)/(2*pi))
    l_index = {}
    for k, (l, m) in enumerate(lm_list):
        l_index.setdefault(l, []).append(k)

    def trace_per_l(M):
        F = np.zeros(args.l_max + 1)
        for l in range(args.l_max + 1):
            idx = l_index[l]
            F[l] = np.real(np.sum(M[idx, idx])) / (2 * np.pi)
        return F

    F_full = trace_per_l(M_full)
    print("\nF_l (standard, unrestricted, from saved raw configs):")
    for l in range(args.l_max + 1):
        print(f"  l={l}: {F_full[l]:.6e}")
    print("\nDelta_l_pair (standard):")
    for l, d in delta_l_pair_from_diag(F_full, args.l_max).items():
        print(f"  l={l}: {d:.5f}")

    if args.uv_cut_theta is not None:
        M_near = compute_Ml_near(spins, w, x, y, z, lm_list, Y, args.uv_cut_theta)
        M_far = M_full - M_near
        F_far = trace_per_l(M_far)
        print(f"\nF_l (UV-excluded, theta_cut={args.uv_cut_theta}, exact pair-restricted M_l):")
        for l in range(args.l_max + 1):
            print(f"  l={l}: {F_far[l]:.6e}")
        print("\nDelta_l_pair (UV-excluded):")
        for l, d in delta_l_pair_from_diag(F_far, args.l_max).items():
            print(f"  l={l}: {d:.5f}")


if __name__ == "__main__":
    main()
