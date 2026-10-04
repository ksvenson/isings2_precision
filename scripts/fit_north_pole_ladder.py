#!/usr/bin/env python3
"""Fit Delta_sigma(a) from bin/analyze_north_pole_configs's output across a
resolution ladder, then extrapolate to the continuum (Delta_inf).

Each <out>.dat file already carries a per-class standard error (Welford
across independent shards, computed inside analyze_north_pole_configs.cc)
-- so unlike fit_real_space_ladder.py, no further shard-level jackknife is
needed here; this script does one err-weighted nonlinear least-squares fit
per ladder point directly.

Usage:
  scripts/fit_north_pole_ladder.py <northpole_q5k4.dat> <northpole_q5k8.dat> ... \\
    [--theta-min 0.3] [--theta-max 2.8] [--power 1]
"""
import argparse
import re

import numpy as np
from scipy.optimize import least_squares


def load(path):
    theta, mean, err = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            theta.append(float(p[1]))
            mean.append(float(p[2]))
            err.append(float(p[3]))
    return np.array(theta), np.array(mean), np.array(err)


def model1(theta, A, Delta):
    return A * (2 - 2 * np.cos(theta)) ** (-Delta)


def fit_delta(theta, mean, err, theta_min, theta_max):
    mask = (theta >= theta_min) & (theta <= theta_max)
    th, y, ye = theta[mask], mean[mask], err[mask]

    def resid(p):
        A, D = p
        return (model1(th, A, D) - y) / ye

    sol = least_squares(resid, (y[len(y) // 2], 0.15), bounds=([0, 0.001], [np.inf, 3.0]))
    chi2 = np.sum(sol.fun ** 2)
    dof = len(th) - 2
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
        d_err = np.sqrt(cov[1, 1])
    except np.linalg.LinAlgError:
        d_err = float("nan")
    return sol.x[1], d_err, chi2 / dof if dof > 0 else float("nan"), mask.sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--theta-min", type=float, default=0.3)
    ap.add_argument("--theta-max", type=float, default=2.8)
    ap.add_argument("--power", type=float, default=1.0)
    args = ap.parse_args()

    points = []
    for path in sorted(args.paths, key=lambda p: int(re.search(r"q5k(\d+)", p).group(1))):
        n_refine = int(re.search(r"q5k(\d+)", path).group(1))
        theta, mean, err = load(path)
        delta, delta_err, rchi2, n_used = fit_delta(theta, mean, err, args.theta_min, args.theta_max)
        points.append((n_refine, delta, delta_err, rchi2, n_used))
        print(f"n_refine={n_refine:>4} Delta={delta:.5f} +/- {delta_err:.5f}  chi2/dof={rchi2:.3f}  n_used={n_used}")

    n_refines = np.array([p[0] for p in points], dtype=float)
    deltas = np.array([p[1] for p in points])
    errs = np.array([p[2] for p in points])

    if len(points) < 3:
        print("\n(need >=3 ladder points for a continuum extrapolation fit)")
        return

    def resid(p):
        d_inf, c = p
        return (d_inf + c * n_refines ** (-args.power) - deltas) / errs

    sol = least_squares(resid, (0.125, 0.0))
    chi2 = np.sum(sol.fun ** 2)
    dof = len(n_refines) - 2
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
        d_inf_err = np.sqrt(cov[0, 0])
    except np.linalg.LinAlgError:
        d_inf_err = float("nan")

    print(f"\nextrapolation model: Delta(a) = Delta_inf + c * n_refine^(-{args.power})")
    print(f"Delta_inf = {sol.x[0]:.5f} +/- {d_inf_err:.5f}   c = {sol.x[1]:.5f}   chi2/dof = {chi2/dof:.3f}")
    print(f"exact CFT value: 1/8 = 0.125   deviation = {(sol.x[0]-0.125):.5f} = {(sol.x[0]-0.125)/d_inf_err:.2f} sigma")


if __name__ == "__main__":
    main()
