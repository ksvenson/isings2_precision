#!/usr/bin/env python3
"""Two-operator fit of F_l(a): F_l = A*F_l^cont(Delta_1) + B*F_l^cont(Delta_2),
a linear superposition of two copies of Owen_section_D's exact recursion
(the Legendre decomposition is linear in the correlator, so a correlator
that is itself an admixture of two power laws (2-2z)^{-Delta_1} and
(2-2z)^{-Delta_2} decomposes into this sum term by term). Built 2026-08-25
after a single-Delta fit was shown to catastrophically fail
(chi2/dof ~7500) while this two-term model fits excellently
(chi2/dof < 1 from n_refine>=16) -- see journal.md's "two-operator fit"
entry for the full writeup and important caveats (diagonal chi-square,
not yet the full jackknife covariance across l; naive a_inf+b/n_refine
continuum extrapolation, not yet validated for these fitted parameters).

Usage:
  scripts/two_operator_fit.py <run_dir> --ladder 8,12,16,24,32,48,64,96,128 --mode naive
"""
import argparse
import glob
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(__file__))
from cft_symmetry_test import analyze

MODE_SUFFIX = {"naive": "", "equal_area": "_eqarea", "equal_rp": "_eqrp"}


def Fcont(Delta, l_max=8):
    F0 = 2.0 ** (1 - 2 * Delta) / (1 - Delta)
    out = [F0]
    for ll in range(1, l_max + 1):
        out.append((ll - 1 + Delta) / (ll + 1 - Delta) * out[-1])
    return np.array(out)


def fit_two_op(F, Ferr):
    l_max = len(F) - 1

    def resid(p):
        A, B, D1, D2 = p
        return (A * Fcont(D1, l_max) + B * Fcont(D2, l_max) - F) / Ferr

    guesses = [
        (F[0] * 0.5, F[0] * 0.5, 0.125, 0.6),
        (F[0] * 0.2, F[0] * 0.8, 0.1, 0.65),
        (F[0] * 0.8, F[0] * 0.2, 0.15, 0.55),
    ]
    best = None
    for g in guesses:
        sol = least_squares(resid, g, bounds=([-np.inf, -np.inf, 0.001, 0.001], [np.inf, np.inf, 0.99, 2.5]))
        chi2 = np.sum(sol.fun ** 2)
        if best is None or chi2 < best[0]:
            best = (chi2, sol)
    chi2, sol = best
    dof = len(F) - 4
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
        errs = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        errs = np.full(4, np.nan)
    return sol.x, errs, chi2 / dof


def pooled_paths(run_dir, k, mode):
    suffix = MODE_SUFFIX[mode]
    return sorted(glob.glob(os.path.join(
        run_dir, f"q5k{k}{suffix}", "shard_*", f"q5k{k}{suffix}", "*_ylm_2pt_full_jackblocks_*.dat")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--ladder", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODE_SUFFIX))
    args = ap.parse_args()

    ladder = [int(x) for x in args.ladder.split(",")]

    print(f"{'k':>5} {'A':>11} {'B':>11} {'Delta_1':>10} {'D1err':>8} {'Delta_2':>10} {'D2err':>8} {'chi2/dof':>9}")
    rows = []
    for k in ladder:
        paths = pooled_paths(args.run_dir, k, args.mode)
        if not paths:
            print(f"WARNING: k={k} no jackblocks files found, skipping", file=sys.stderr)
            continue
        r = analyze(paths)
        x, errs, chi2dof = fit_two_op(r["F_full"], r["F_err"])
        A, B, D1, D2 = x
        print(f"{k:5d} {A:11.4e} {B:11.4e} {D1:10.5f} {errs[2]:8.5f} {D2:10.5f} {errs[3]:8.5f} {chi2dof:9.3f}")
        rows.append((k, D1, errs[2], D2, errs[3], chi2dof))

    arr = np.array([r[:5] for r in rows])
    mask = arr[:, 0] >= 16
    if mask.sum() >= 3:
        k = arr[mask, 0]
        print()
        print("# continuum extrapolation (a_inf + b/n_refine, n_refine>=16 only -- naive ansatz, see journal.md caveats):")
        for name, y, yerr in [("Delta_1", arr[mask, 1], arr[mask, 2]), ("Delta_2", arr[mask, 3], arr[mask, 4])]:
            X = np.vstack([np.ones_like(k), 1.0 / k]).T
            W = np.diag(1.0 / yerr ** 2)
            cov = np.linalg.inv(X.T @ W @ X)
            beta = cov @ X.T @ W @ y
            resid = y - X @ beta
            chi2 = resid @ W @ resid
            dof = len(k) - 2
            err = np.sqrt(np.diag(cov))
            print(f"#   {name}(a->0) = {beta[0]:.5f} +/- {err[0]:.5f}  slope={beta[1]:+.3f}+/-{err[1]:.3f}  chi2/dof={chi2/dof:.2f}")


if __name__ == "__main__":
    main()
