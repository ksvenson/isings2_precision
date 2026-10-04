#!/usr/bin/env python3
"""Spherical-symmetry test from a *_ylm_2pt_full_jackblocks_*.dat file.

Reads the block-jackknife checkpoint written by
../IsingS2/src/ising_s2_crit.cc (--jack_block_size), reconstructs the full
Hermitian M_l[m,m'] = <S_lm S*_lm'> matrix per harmonic level l via
leave-one-block-out jackknife, and reports two families of SO(3)-breaking
diagnostics per l:

  R_l    = sum_{m!=m'} |M_l[m,m']|^2 / sum_m |M_l[m,m]|^2
           (relative off-diagonal power; ->0 for an SO(3)-invariant theory)
  chi2_l = generalized-least-squares chi-square of the (2l+1) diagonal
           entries against a constant, using the jackknife covariance
           matrix across m (accounts for correlations between different m
           at fixed l, unlike a plain coefficient of variation)

  kappa{2,3,4}_r = reduced spectral cumulants of M_l's own eigenvalue
           distribution, built from Tr(M_l^k) (k=1..4) via the standard
           moment-to-cumulant formulas and normalized by kappa1^k. Under
           exact SO(3) symmetry M_l is proportional to the identity, so
           *every* eigenvalue is equal and every spectral cumulant of
           order >=2 vanishes exactly -- not just the variance. Unlike
           R_l/chi2_l, these are basis-independent: Tr(M^k) doesn't care
           which quantization axis the mesh happened to fix, so they can't
           be fooled by that axis being accidentally aligned (or
           misaligned) with a mesh symmetry direction. R_l/chi2_l are kept
           because they test something different -- whether the
           *specific* mesh-fixed axis looks anomalous.

All point estimates are computed from the full dataset; all errors come
from the jackknife blocks (leave-one-block-out). See
IsingS2_precision/journal.md 2026-08-22 for the design rationale (block
size is not assumed small a priori -- it is whatever --jack_block_size the
run used; kappa_r motivation is the same date's later entry).
"""

import argparse
import re
import sys
from collections import defaultdict

import numpy as np

HEADER_RE = re.compile(r"#\s*n_pairs=(\d+)\s+jack_block_size=(\d+)\s+n_blocks=(\d+)")


def read_jackblocks(path):
    """Returns (l_max, blocks) where blocks is a list of dicts:
    {'n': int, 'sum': {(l, m, mp): complex}}
    """
    with open(path) as f:
        header = f.readline()
        m = HEADER_RE.match(header)
        if not m:
            raise ValueError(f"{path}: unrecognized header {header!r}")
        n_pairs, jack_block_size, n_blocks = (int(x) for x in m.groups())

        blocks = [dict(n=None, sum={}) for _ in range(n_blocks)]
        l_max = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            k, l, mm, mp, sum_re, sum_im, n_k = line.split()
            k, l, mm, mp, n_k = int(k), int(l), int(mm), int(mp), int(n_k)
            sum_re, sum_im = float(sum_re), float(sum_im)
            l_max = max(l_max, l)
            blk = blocks[k]
            if blk["n"] is None:
                blk["n"] = n_k
            elif blk["n"] != n_k:
                raise ValueError(f"{path}: block {k} has inconsistent n ({blk['n']} vs {n_k})")
            blk["sum"][(l, mm, mp)] = complex(sum_re, sum_im)

    for k, blk in enumerate(blocks):
        if len(blk["sum"]) != n_pairs:
            raise ValueError(
                f"{path}: block {k} has {len(blk['sum'])} pairs, expected {n_pairs}"
            )
    return l_max, jack_block_size, blocks


def full_pairs_for_l(l):
    """(m, mp) with m <= mp, matching the C++ generation order."""
    return [(mm, mp) for mm in range(-l, l + 1) for mp in range(mm, l + 1)]


def hermitian_matrix(l, mean_upper):
    """Build the full (2l+1)x(2l+1) M_l matrix from upper-triangle means.

    mean_upper: dict {(m, mp): complex} for m <= mp.
    Returns a numpy array indexed [m+l, mp+l].
    """
    dim = 2 * l + 1
    mat = np.zeros((dim, dim), dtype=complex)
    for (mm, mp), val in mean_upper.items():
        mat[mm + l, mp + l] = val
        if mm != mp:
            mat[mp + l, mm + l] = np.conj(val)
    return mat


def off_diag_ratio(mat):
    diag = np.real(np.diag(mat))
    denom = np.sum(diag ** 2)
    off = mat - np.diag(np.diag(mat))
    numer = np.sum(np.abs(off) ** 2)
    if denom == 0.0:
        return np.nan
    return numer / denom


def matrix_trace_powers(mat, k_max=4):
    """[Tr(M), Tr(M^2), ..., Tr(M^k_max)], via repeated matmul."""
    powers = []
    mp = mat
    for k in range(1, k_max + 1):
        powers.append(np.real(np.trace(mp)))
        if k < k_max:
            mp = mp @ mat
    return powers


def moments_to_cumulants(p1, p2, p3, p4):
    """Cumulants kappa1..4 of a distribution from its raw moments p1..4."""
    mu2 = p2 - p1 ** 2
    mu3 = p3 - 3 * p1 * p2 + 2 * p1 ** 3
    mu4 = p4 - 4 * p1 * p3 + 6 * p1 ** 2 * p2 - 3 * p1 ** 4
    kappa1 = p1
    kappa2 = mu2
    kappa3 = mu3
    kappa4 = mu4 - 3 * mu2 ** 2
    return kappa1, kappa2, kappa3, kappa4


def reduced_spectral_cumulants(mat):
    """(kappa2_r, kappa3_r, kappa4_r) of mat's eigenvalue distribution.

    Computed from Tr(M^k)/dim (k=1..4) without diagonalizing, then
    kappa_j normalized by kappa1^j so different l are comparable. All
    exactly zero iff mat is proportional to the identity.
    """
    dim = mat.shape[0]
    t1, t2, t3, t4 = matrix_trace_powers(mat, k_max=4)
    p1, p2, p3, p4 = t1 / dim, t2 / dim, t3 / dim, t4 / dim
    kappa1, kappa2, kappa3, kappa4 = moments_to_cumulants(p1, p2, p3, p4)
    if kappa1 == 0.0:
        return np.nan, np.nan, np.nan
    return kappa2 / kappa1 ** 2, kappa3 / kappa1 ** 3, kappa4 / kappa1 ** 4


def jackknife_error(theta_full, theta_samples):
    n = len(theta_samples)
    theta_bar = np.mean(theta_samples)
    var = (n - 1) / n * np.sum((theta_samples - theta_bar) ** 2)
    return theta_full, np.sqrt(var)


def diag_chi2_const(d_full, cov, min_blocks_factor=5):
    """GLS chi-square of d_full (real vector) against a constant, given the
    jackknife covariance matrix cov. Returns (chi2, dof, c_hat) or
    (None, dof, None) if cov is judged too poorly-conditioned to trust
    (fewer than min_blocks_factor * dim jackknife samples went into it --
    a rule of thumb, not a derived bound; flagged rather than silently
    trusted, same spirit as the rest of this campaign's statistics policy).
    """
    dim = len(d_full)
    dof = dim - 1
    try:
        w = np.linalg.pinv(cov)
    except np.linalg.LinAlgError:
        return None, dof, None
    ones = np.ones(dim)
    denom = ones @ w @ ones
    if denom == 0.0 or not np.isfinite(denom):
        return None, dof, None
    c_hat = (ones @ w @ d_full) / denom
    resid = d_full - c_hat * ones
    chi2 = resid @ w @ resid
    return chi2, dof, c_hat


def read_jackblocks_multi(paths):
    """Pool multiple shard jackblocks files (same n_pairs/l_max/jack_block_size)
    from independent seeded runs of the same (n_refine, coupling_rule) point
    into one combined block list, by simple concatenation -- each shard's
    per-block sums are already independent, so the concatenation is itself a
    valid (larger) set of leave-one-block-out jackknife blocks.
    """
    l_max = jack_block_size = None
    all_blocks = []
    for path in paths:
        l_max_i, jbs_i, blocks_i = read_jackblocks(path)
        if l_max is None:
            l_max, jack_block_size = l_max_i, jbs_i
        elif l_max_i != l_max:
            raise ValueError(f"{path}: l_max={l_max_i} != {l_max} from earlier shard")
        elif jbs_i != jack_block_size:
            raise ValueError(
                f"{path}: jack_block_size={jbs_i} != {jack_block_size} from earlier shard"
            )
        all_blocks.extend(blocks_i)
    return l_max, jack_block_size, all_blocks


def analyze_blocks(l_max, blocks, min_blocks_factor=5):
    n_blocks = len(blocks)
    total_n = sum(blk["n"] for blk in blocks)

    results = []
    for l in range(l_max + 1):
        pairs = full_pairs_for_l(l)
        dim = 2 * l + 1

        total_sum = {p: sum(blk["sum"][(l,) + p] for blk in blocks) for p in pairs}
        mean_full = {p: v / total_n for p, v in total_sum.items()}
        mat_full = hermitian_matrix(l, mean_full)
        r_full = off_diag_ratio(mat_full)
        d_full = np.real(np.diag(mat_full))
        k2_full, k3_full, k4_full = reduced_spectral_cumulants(mat_full)

        r_jk = np.empty(n_blocks)
        d_jk = np.empty((n_blocks, dim))
        k2_jk = np.empty(n_blocks)
        k3_jk = np.empty(n_blocks)
        k4_jk = np.empty(n_blocks)
        for k, blk in enumerate(blocks):
            n_k = blk["n"]
            loo_n = total_n - n_k
            mean_loo = {p: (total_sum[p] - blk["sum"][(l,) + p]) / loo_n for p in pairs}
            mat_loo = hermitian_matrix(l, mean_loo)
            r_jk[k] = off_diag_ratio(mat_loo)
            d_jk[k] = np.real(np.diag(mat_loo))
            k2_jk[k], k3_jk[k], k4_jk[k] = reduced_spectral_cumulants(mat_loo)

        _, r_err = jackknife_error(r_full, r_jk)
        _, k2_err = jackknife_error(k2_full, k2_jk)
        _, k3_err = jackknife_error(k3_full, k3_jk)
        _, k4_err = jackknife_error(k4_full, k4_jk)

        # kappa3/kappa4 need enough eigenvalues to be distinct from lower
        # moments -- below that, report as unavailable rather than a
        # numerically noisy near-zero.
        if dim < 3:
            k3_full = k3_err = None
        if dim < 4:
            k4_full = k4_err = None

        d_bar = np.mean(d_jk, axis=0)
        cov = (n_blocks - 1) / n_blocks * (d_jk - d_bar).T @ (d_jk - d_bar)

        if n_blocks >= min_blocks_factor * dim:
            chi2, dof, c_hat = diag_chi2_const(d_full, cov)
        else:
            chi2, dof, c_hat = None, dim - 1, None

        results.append(
            dict(
                l=l, dim=dim, R_l=r_full, R_l_err=r_err, chi2=chi2, dof=dof, c_hat=c_hat,
                kappa2_r=k2_full, kappa2_r_err=k2_err,
                kappa3_r=k3_full, kappa3_r_err=k3_err,
                kappa4_r=k4_full, kappa4_r_err=k4_err,
            )
        )

    return n_blocks, total_n, results


def analyze(path, min_blocks_factor=5):
    l_max, jack_block_size, blocks = read_jackblocks(path)
    n_blocks, total_n, results = analyze_blocks(l_max, blocks, min_blocks_factor)
    return l_max, jack_block_size, n_blocks, total_n, results


def analyze_multi(paths, min_blocks_factor=5):
    l_max, jack_block_size, blocks = read_jackblocks_multi(paths)
    n_blocks, total_n, results = analyze_blocks(l_max, blocks, min_blocks_factor)
    return l_max, jack_block_size, n_blocks, total_n, results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "jackblocks_path", nargs="+",
        help="*_ylm_2pt_full_jackblocks_*.dat file(s); multiple shard files "
        "for the same (n_refine, coupling_rule) point are pooled",
    )
    ap.add_argument(
        "--min-blocks-factor",
        type=int,
        default=5,
        help="skip the chi2 test at l if n_blocks < this * (2l+1) (default 5, a rule of "
        "thumb for covariance-matrix stability, not a derived bound)",
    )
    args = ap.parse_args()

    if len(args.jackblocks_path) == 1:
        l_max, jack_block_size, n_blocks, total_n, results = analyze(
            args.jackblocks_path[0], args.min_blocks_factor
        )
    else:
        l_max, jack_block_size, n_blocks, total_n, results = analyze_multi(
            args.jackblocks_path, args.min_blocks_factor
        )

    print(f"# {args.jackblocks_path}")
    print(f"# l_max={l_max} jack_block_size={jack_block_size} n_blocks={n_blocks} "
          f"total_samples={total_n}")
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

        def fmt(key, err_key):
            v, e = r[key], r[err_key]
            if v is None or e is None:
                return f"{'--':>14} {'--':>14}"
            return f"{v:14.6e} {e:14.6e}"

        print(
            f"{r['l']:3d} {r['dim']:4d} {r['R_l']:14.6e} {r['R_l_err']:14.6e} "
            f"{chi2_str} {r['dof']:4d} {rchi2_str} {c_hat_str} "
            f"{fmt('kappa2_r', 'kappa2_r_err')} "
            f"{fmt('kappa3_r', 'kappa3_r_err')} "
            f"{fmt('kappa4_r', 'kappa4_r_err')}"
        )

    skipped = [r["l"] for r in results if r["chi2"] is None]
    if skipped:
        print(
            f"# chi2 skipped for l={skipped}: n_blocks < {args.min_blocks_factor} * (2l+1) "
            "-- covariance matrix not trusted at this block count",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
