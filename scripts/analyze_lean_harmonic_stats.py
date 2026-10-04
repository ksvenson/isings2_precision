#!/usr/bin/env python3
"""Analyze bin/lean_harmonic_stats output: F_l/Delta_s/Delta_l_pair from the
online-accumulated diagonal |S_lm|^2 sums, plus magnetization moments
(<M>, <M^2>, <M^3>, <M^4>, susceptibility chi = (<M^2>-<M>^2)/n_sites,
Binder U4 = 1 - <M^4>/(3<M^2>^2)) -- all from per-jackknife-block
sums, no raw configs or full M_l matrix involved. See
src/lean_harmonic_stats.cc's header for what this driver accumulates.

Usage:
  scripts/analyze_lean_harmonic_stats.py <lean1.dat> [<lean2.dat> ...] [--l_max 8]
"""
import argparse
import re
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from cft_symmetry_test import delta_l_pair_formula, delta_s_formula  # noqa: E402


def load(path):
    """Parses both the old (diagonal-only) and new (2026-08-27, diagonal +
    off-diagonal) lean_harmonic_stats formats. `offdiag_list`/blocks'
    `offdiag_re`/`offdiag_im` are None/absent for old-format files (no
    `# offdiag_order:` header line, `n_offdiag` unset) -- callers that only
    need F_l/Delta_s (the original use case) are unaffected either way.
    """
    header = {}
    lm_list = None
    offdiag_list = None
    blocks = []
    with open(path) as f:
        for line in f:
            if line.startswith("# lm_order:"):
                lm_list = [tuple(int(v) for v in tok.strip("()").split(","))
                           for tok in re.findall(r"\(-?\d+,-?\d+\)", line)]
                continue
            if line.startswith("# offdiag_order:"):
                offdiag_list = [tuple(int(v) for v in tok.strip("()").split(","))
                                 for tok in re.findall(r"\(-?\d+,-?\d+,-?\d+\)", line)]
                continue
            if line.startswith("#"):
                for kv in re.findall(r"(\w+)=(\S+)", line):
                    header[kv[0]] = kv[1]
                continue
            if not line.strip():
                continue
            p = line.split()
            try:
                n_meas = int(p[1])
                M_sum, M2_sum, M3_sum, M4_sum = (float(x) for x in p[2:6])
                n_lm = len(lm_list)
                rest = [float(x) for x in p[6:]]
            except (ValueError, IndexError):
                # Torn trailing write from a still-running production task
                # appending to this file concurrently with this read --
                # drop the incomplete final block rather than crash (seen
                # repeatedly reading eqarea/eqrp dirs while job 7341060 was
                # still writing, 2026-08-28).
                break
            abs2 = np.array(rest[:n_lm])
            blk = dict(n=n_meas, M=M_sum, M2=M2_sum, M3=M3_sum, M4=M4_sum, abs2=abs2)
            if offdiag_list is not None:
                n_off = len(offdiag_list)
                offdiag_re = np.array(rest[n_lm:n_lm + n_off])
                offdiag_im = np.array(rest[n_lm + n_off:n_lm + 2 * n_off])
                blk["offdiag_re"] = offdiag_re
                blk["offdiag_im"] = offdiag_im
            blocks.append(blk)
    return header, lm_list, blocks, offdiag_list


def to_analyze_blocks_format(lm_list, offdiag_list, blocks):
    """Converts lean_harmonic_stats blocks (abs2 + offdiag_re/im arrays) into
    the {n, sum: {(l,mm,mp): complex}} structure scripts/symmetry_test.py's
    analyze_blocks() expects, so the REAL kappa2_r/R_l (needs the full
    M_l[m,m'] matrix, not just the diagonal) reuses that already-verified
    code instead of being reimplemented. Requires offdiag_list (i.e. the
    2026-08-27+ format) -- raises if the file predates that.
    """
    if offdiag_list is None:
        raise ValueError("this file has no off-diagonal data (old lean_harmonic_stats format) "
                          "-- real kappa2_r needs a file written after the 2026-08-27 mesh/offdiag update")
    out_blocks = []
    for blk in blocks:
        s = {}
        for (l, m), val in zip(lm_list, blk["abs2"]):
            s[(l, m, m)] = complex(val, 0.0)
        for (l, m, mp), re_v, im_v in zip(offdiag_list, blk["offdiag_re"], blk["offdiag_im"]):
            s[(l, m, mp)] = complex(re_v, im_v)
        out_blocks.append(dict(n=blk["n"], sum=s))
    return out_blocks


def combine(paths):
    all_blocks = []
    header = lm_list = offdiag_list = None
    for path in paths:
        h, lm, blocks, offdiag = load(path)
        if lm_list is None:
            header, lm_list, offdiag_list = h, lm, offdiag
        elif lm != lm_list:
            raise ValueError(f"{path}: lm_order mismatch across files")
        elif offdiag != offdiag_list:
            raise ValueError(f"{path}: offdiag_order mismatch across files")
        all_blocks.extend(blocks)
    return header, lm_list, all_blocks, offdiag_list


def jackknife_err(full, jk):
    n = len(jk)
    bar = np.mean(jk, axis=0)
    return np.sqrt((n - 1) / n * np.sum((jk - bar) ** 2, axis=0))


def analyze(paths, l_max=None):
    header, lm_list, blocks, _offdiag_list = combine(paths)
    n_blocks = len(blocks)
    total_n = sum(b["n"] for b in blocks)
    l_of_k = np.array([lm[0] for lm in lm_list])
    l_max_data = l_of_k.max()
    if l_max is None:
        l_max = l_max_data

    abs2_stack = np.array([b["abs2"] for b in blocks])  # (n_blocks, n_lm)
    n_stack = np.array([b["n"] for b in blocks])

    total_abs2 = abs2_stack.sum(axis=0)
    # F_l = Trace(M_l)/(2*pi)/(2l+1), Trace(M_l) = sum_m |S_lm|^2 summed
    # over configs / n_configs -- the /(2l+1) is REQUIRED to match
    # F_l^cont (see cft_symmetry_test.py's compute_F_l docstring, 2026-08-26
    # "F_l(a) missing (2l+1) division" fix, for the full derivation:
    # Trace(M_l)/(2*pi) equals (2l+1)*F_l^cont in the continuum limit, not
    # F_l^cont itself, since this is a double (bilinear) sum over site
    # pairs, not a single-reference-row projection).
    def F_l_from_abs2(abs2_sum, n):
        F = np.zeros(l_max + 1)
        for l in range(l_max + 1):
            mask = l_of_k == l
            F[l] = abs2_sum[mask].sum() / n / (2.0 * np.pi) / (2 * l + 1)
        return F

    F_full = F_l_from_abs2(total_abs2, total_n)
    F_jk = np.empty((n_blocks, l_max + 1))
    for k in range(n_blocks):
        loo_abs2 = total_abs2 - abs2_stack[k]
        loo_n = total_n - n_stack[k]
        F_jk[k] = F_l_from_abs2(loo_abs2, loo_n)
    F_err = jackknife_err(F_full, F_jk)

    ds_full = delta_s_formula(F_full[0], F_full[1])
    ds_jk = delta_s_formula(F_jk[:, 0], F_jk[:, 1])
    ds_err = jackknife_err(ds_full, ds_jk)

    # magnetization moments + Binder cumulant, same jackknife-over-blocks pattern
    n_sites = int(header.get("n_sites", 0))
    M_tot = sum(b["M"] for b in blocks)
    M2_tot = sum(b["M2"] for b in blocks)
    M3_tot = sum(b["M3"] for b in blocks)
    M4_tot = sum(b["M4"] for b in blocks)
    m1_full = M_tot / total_n
    m2_full = M2_tot / total_n
    m3_full = M3_tot / total_n
    m4_full = M4_tot / total_n
    u4_full = 1.0 - m4_full / (3.0 * m2_full ** 2)
    # susceptibility per site: chi = (<M^2> - <M>^2) / n_sites (standard
    # finite-size-scaling normalization -- M is extensive so <M^2> ~ n_sites,
    # dividing by n_sites makes chi intensive)
    chi_full = (m2_full - m1_full ** 2) / n_sites if n_sites else np.nan

    m1_jk = np.array([(M_tot - b["M"]) / (total_n - b["n"]) for b in blocks])
    m2_jk = np.array([(M2_tot - b["M2"]) / (total_n - b["n"]) for b in blocks])
    m3_jk = np.array([(M3_tot - b["M3"]) / (total_n - b["n"]) for b in blocks])
    m4_jk = np.array([(M4_tot - b["M4"]) / (total_n - b["n"]) for b in blocks])
    u4_jk = 1.0 - m4_jk / (3.0 * m2_jk ** 2)
    chi_jk = (m2_jk - m1_jk ** 2) / n_sites if n_sites else np.full(n_blocks, np.nan)
    m1_err = jackknife_err(m1_full, m1_jk)
    m2_err = jackknife_err(m2_full, m2_jk)
    m3_err = jackknife_err(m3_full, m3_jk)
    m4_err = jackknife_err(m4_full, m4_jk)
    u4_err = jackknife_err(u4_full, u4_jk)
    chi_err = jackknife_err(chi_full, chi_jk)

    print(f"# {paths}")
    print(f"# n_sites={header.get('n_sites')} l_max={l_max} n_blocks={n_blocks} total_meas={total_n}")
    print(f"# <M>   = {m1_full:.6e} +/- {m1_err:.6e}")
    print(f"# <M^2> = {m2_full:.6e} +/- {m2_err:.6e}")
    print(f"# <M^3> = {m3_full:.6e} +/- {m3_err:.6e}")
    print(f"# <M^4> = {m4_full:.6e} +/- {m4_err:.6e}")
    print(f"# chi = (<M^2>-<M>^2)/n_sites = {chi_full:.6e} +/- {chi_err:.6e}")
    print(f"# Binder U4 = {u4_full:.6f} +/- {u4_err:.6f}")
    print(f"# Delta_s(a) = 2*F_1/(F_1+F_0) [exact value 1/8 = 0.125]:")
    print(f"#   {ds_full:.6f} +/- {ds_err:.6f}")
    print(f"# F_l(a):")
    for l in range(l_max + 1):
        print(f"#   l={l:2d}  F_l={F_full[l]:14.6e} +/- {F_err[l]:14.6e}")
    print(f"{'l':>3} {'Delta_l_pair':>16} {'err':>16}")
    for l in range(1, l_max + 1):
        dp_full = delta_l_pair_formula(F_full[l - 1], F_full[l], l)
        dp_jk = delta_l_pair_formula(F_jk[:, l - 1], F_jk[:, l], l)
        dp_err = jackknife_err(dp_full, dp_jk)
        print(f"{l:3d} {dp_full:16.6e} {dp_err:16.6e}")

    dp_full_arr = np.empty(l_max)
    dp_err_arr = np.empty(l_max)
    for l in range(1, l_max + 1):
        dp_full_arr[l - 1] = delta_l_pair_formula(F_full[l - 1], F_full[l], l)
        dp_jk = delta_l_pair_formula(F_jk[:, l - 1], F_jk[:, l], l)
        dp_err_arr[l - 1] = jackknife_err(dp_full_arr[l - 1], dp_jk)

    return dict(
        n_sites=int(header.get("n_sites", 0)), n_refine=int(header.get("n_refine", 0)),
        l_max=l_max, n_blocks=n_blocks, total_n=total_n,
        F_full=F_full, F_err=F_err,
        delta_s=ds_full, delta_s_err=ds_err,
        delta_l_pair=dp_full_arr, delta_l_pair_err=dp_err_arr,  # index l-1 -> l=1..l_max
        m1=m1_full, m1_err=m1_err,
        m2=m2_full, m2_err=m2_err,
        m3=m3_full, m3_err=m3_err,
        m4=m4_full, m4_err=m4_err,
        chi=chi_full, chi_err=chi_err,
        u4=u4_full, u4_err=u4_err,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--l_max", type=int, default=None)
    args = ap.parse_args()
    analyze(args.paths, args.l_max)


if __name__ == "__main__":
    main()
