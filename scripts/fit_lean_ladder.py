#!/usr/bin/env python3
"""Continuum-extrapolation fit of Delta_s(a) from bin/lean_harmonic_stats
ladder output, same pattern as fit_real_space_ladder.py but reusing
analyze_lean_harmonic_stats.analyze()'s already-correct per-point
jackknife error (no separate shard-level jackknife needed here -- the
driver's own per-block accumulation already gives one).

Fits Delta(a) = Delta_inf + c * n_refine^(-power). Per this campaign's
established finding (journal.md "CLEAN CONFIRMATION", 2026-08-25/26),
--power's correct value is NOT assumed -- both p=1 (linear) and p=2
(quadratic) are fit and reported side by side; prefer whichever has the
better chi2/dof, don't just trust the default.

Usage:
  scripts/fit_lean_ladder.py <dir_with_*_lean_*.dat files> [--l_max 8]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_lean_harmonic_stats import analyze  # noqa: E402


def fit_power(n_refines, deltas, errs, power):
    def resid(p):
        d_inf, c = p
        return (d_inf + c * n_refines ** (-power) - deltas) / errs

    sol = least_squares(resid, (0.125, 0.0))
    chi2 = np.sum(sol.fun ** 2)
    dof = len(n_refines) - 2
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
        d_inf_err = np.sqrt(cov[0, 0])
    except np.linalg.LinAlgError:
        d_inf_err = float("nan")
    return sol.x[0], d_inf_err, sol.x[1], chi2 / dof


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--min_n_refine", type=int, default=None,
                     help="drop ladder points below this n_refine (exclude small-L aliasing/finite-size outliers from the fit)")
    args = ap.parse_args()

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")

    by_n_refine = {}
    for f in files:
        n_refine = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n_refine.setdefault(n_refine, []).append(f)

    n_refines, deltas, errs = [], [], []
    print(f"{'n_refine':>9} {'n_shards':>9} {'Delta_s':>12} {'err':>10}")
    for n_refine in sorted(by_n_refine):
        if args.min_n_refine and n_refine < args.min_n_refine:
            continue
        shard_files = sorted(by_n_refine[n_refine])
        r = analyze(shard_files, args.l_max)
        print(f"{n_refine:>9} {len(shard_files):>9} {r['delta_s']:>12.5f} {r['delta_s_err']:>10.5f}")
        n_refines.append(n_refine)
        deltas.append(r["delta_s"])
        errs.append(r["delta_s_err"])

    n_refines = np.array(n_refines, dtype=float)
    deltas = np.array(deltas)
    errs = np.array(errs)

    print(f"\nextrapolation model: Delta(a) = Delta_inf + c * n_refine^(-power)")
    for power in (1.0, 2.0):
        d_inf, d_inf_err, c, chi2_dof = fit_power(n_refines, deltas, errs, power)
        sigma = (d_inf - 0.125) / d_inf_err if d_inf_err == d_inf_err else float("nan")
        print(f"power={power:.0f}: Delta_inf = {d_inf:.5f} +/- {d_inf_err:.5f}   c = {c:.5f}   "
              f"chi2/dof = {chi2_dof:.3f}   deviation = {sigma:.2f} sigma")


if __name__ == "__main__":
    main()
