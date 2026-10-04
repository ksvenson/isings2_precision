#!/usr/bin/env python3
"""Continuum-extrapolate Delta_sigma(a) from a resolution ladder of
real_space_2pt_test production shards (cluster/sge/submit_real_space_ladder.sh).

Each ladder point (n_refine) has N_SHARDS independent-seed shards, each
producing one <run_id>_real_space_2pt_<seed>.dat file (real_space_2pt_test
has no in-C++ jackknife -- each shard is an independent MC chain, used
directly as the jackknife unit here via delete-one-shard resampling).

**Updated 2026-08-26** for real_space_2pt_test.cc's rewrite that dropped
z=cos(theta) histogram binning for exact icosahedral equivalence classes
(see that file's header comment). Each shard now picks its own random
reference sites (seed-dependent), so shards no longer share a common theta
grid the way fixed z-bins did -- there is nothing to average bin-by-bin
across shards any more. Instead, each fold's fit is a single weighted
nonlinear least-squares fit over the UNION of every included shard's
(theta, mean, err) rows, weighted by each row's own online-accumulator
error (real_space_2pt_test.cc's QfeMeasReal.Error()) -- concatenation, not
averaging, since different shards' classes are (mostly) different pairs,
not repeated measurements of the same pair.

For each n_refine:
  1. Load every shard's classes.
  2. Delete-one-shard jackknife: for k = 0..n_shards-1, weighted-fit
     Delta_k over the concatenated rows of every shard except k.
  3. Jackknife mean/error of Delta from the n_shards delete-one fits.

Then fits Delta(a) = Delta_inf + c * n_refine**(-p) (p fixed, default 2)
across the ladder, weighted by the jackknife errors, to extract the
continuum-limit Delta_inf.

Usage:
  scripts/fit_real_space_ladder.py <campaign_runs/real_space_ladder_.../>
    [--theta-min 0.3] [--theta-max 2.8] [--power 2]
"""
import argparse
import glob
import os
import re

import numpy as np
from scipy.optimize import least_squares


def load_classes(path):
    theta, mean, err = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            theta.append(float(p[2]))
            mean.append(float(p[3]))
            err.append(float(p[4]))
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
    return sol.x[1]  # Delta


def analyze_point(run_dir, theta_min, theta_max):
    shard_dirs = sorted(glob.glob(os.path.join(run_dir, "shard_*")))
    shard_data = []
    for sd in shard_dirs:
        dat_files = glob.glob(os.path.join(sd, "*_real_space_2pt_*.dat"))
        if not dat_files:
            continue
        shard_data.append(load_classes(dat_files[0]))
    n_shards = len(shard_data)
    if n_shards < 2:
        raise RuntimeError(f"{run_dir}: need >=2 shards, found {n_shards}")

    def pooled(indices):
        theta = np.concatenate([shard_data[k][0] for k in indices])
        mean = np.concatenate([shard_data[k][1] for k in indices])
        err = np.concatenate([shard_data[k][2] for k in indices])
        return theta, mean, err

    all_idx = list(range(n_shards))
    delta_jk = np.empty(n_shards)
    for k in range(n_shards):
        loo_theta, loo_mean, loo_err = pooled([i for i in all_idx if i != k])
        delta_jk[k] = fit_delta(loo_theta, loo_mean, loo_err, theta_min, theta_max)

    delta_bar = delta_jk.mean()
    delta_err = np.sqrt((n_shards - 1) / n_shards * np.sum((delta_jk - delta_bar) ** 2))
    # full-statistics point estimate (not the jackknife mean, for reference)
    full_theta, full_mean, full_err = pooled(all_idx)
    delta_full = fit_delta(full_theta, full_mean, full_err, theta_min, theta_max)
    return n_shards, delta_full, delta_bar, delta_err


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_root")
    ap.add_argument("--theta-min", type=float, default=0.3)
    ap.add_argument("--theta-max", type=float, default=2.8)
    ap.add_argument("--power", type=float, default=2.0)
    args = ap.parse_args()

    point_dirs = sorted(
        glob.glob(os.path.join(args.run_root, "q5k*")),
        key=lambda p: int(re.search(r"q5k(\d+)", p).group(1)),
    )

    n_refines, deltas, errs = [], [], []
    print(f"{'n_refine':>9} {'n_shards':>9} {'Delta_full':>12} {'Delta_jk':>12} {'err':>10}")
    for pd in point_dirs:
        n_refine = int(re.search(r"q5k(\d+)", pd).group(1))
        n_shards, delta_full, delta_bar, delta_err = analyze_point(pd, args.theta_min, args.theta_max)
        print(f"{n_refine:>9} {n_shards:>9} {delta_full:>12.5f} {delta_bar:>12.5f} {delta_err:>10.5f}")
        n_refines.append(n_refine)
        deltas.append(delta_bar)
        errs.append(delta_err)

    n_refines = np.array(n_refines, dtype=float)
    deltas = np.array(deltas)
    errs = np.array(errs)

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
