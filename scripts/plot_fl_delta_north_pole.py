#!/usr/bin/env python3
"""Plots Delta_s(a) from fl_delta_from_north_pole.py (the corrected,
no-(2l+1) l-space F_l projection of the north-pole two-point function)
vs 1/n_refine, with a continuum extrapolation fit, alongside the
real-space (fit_north_pole_ladder.py) Delta(a) trend for comparison.

Usage:
  scripts/plot_fl_delta_north_pole.py <analysis_dir> -o <out_path> [--power 2]
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cft_symmetry_test import delta_s_formula  # noqa: E402
from fit_north_pole_ladder import fit_delta as fit_delta_realspace  # noqa: E402
from fit_north_pole_ladder import load as load_northpole  # noqa: E402
from fl_delta_from_north_pole import compute_Fl, compute_Fl_jk  # noqa: E402
from fl_delta_from_north_pole import load as load_northpole_full  # noqa: E402


def fl_series(paths):
    rng = np.random.default_rng(0)
    n_refines, deltas, errs = [], [], []
    for path in sorted(paths, key=lambda p: int(re.search(r"q5k(\d+)", p).group(1))):
        n_refine = int(re.search(r"q5k(\d+)", path).group(1))
        n_sites, theta, mean, sem, n_targets, sum_wt = load_northpole_full(path)
        Fl = compute_Fl(n_sites, theta, mean, sum_wt, 1)
        ds = delta_s_formula(Fl[0], Fl[1])
        Fl_samples = compute_Fl_jk(n_sites, theta, mean, sem, sum_wt, 1, 200, rng)
        ds_samples = delta_s_formula(Fl_samples[:, 0], Fl_samples[:, 1])
        n_refines.append(n_refine)
        deltas.append(ds)
        errs.append(ds_samples.std())
    return np.array(n_refines, dtype=float), np.array(deltas), np.array(errs)


def realspace_series(paths, theta_min=0.3, theta_max=2.8):
    n_refines, deltas, errs = [], [], []
    for path in sorted(paths, key=lambda p: int(re.search(r"q5k(\d+)", p).group(1))):
        n_refine = int(re.search(r"q5k(\d+)", path).group(1))
        theta, mean, err = load_northpole(path)
        d, de, _, _ = fit_delta_realspace(theta, mean, err, theta_min, theta_max)
        n_refines.append(n_refine)
        deltas.append(d)
        errs.append(de)
    return np.array(n_refines, dtype=float), np.array(deltas), np.array(errs)


def extrapolate(n_refines, deltas, errs, power):
    def resid(p):
        d_inf, c = p
        return (d_inf + c * n_refines ** (-power) - deltas) / errs

    sol = least_squares(resid, (0.125, 0.0))
    chi2 = np.sum(sol.fun ** 2)
    dof = max(len(n_refines) - 2, 1)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
        d_inf_err = np.sqrt(cov[0, 0])
    except np.linalg.LinAlgError:
        d_inf_err = float("nan")
    return sol.x[0], d_inf_err, sol.x[1], chi2 / dof


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("analysis_dir")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--power", type=float, default=2.0)
    ap.add_argument("--min-n-refine-fl", type=int, default=8,
                     help="drop ladder points below this for the F_l fit (n_refine=4 is too coarse)")
    args = ap.parse_args()

    paths = [p for p in glob.glob(os.path.join(args.analysis_dir, "northpole_q5k*.dat")) if os.path.getsize(p) > 0]

    nr_fl, d_fl, e_fl = fl_series(paths)
    mask = nr_fl >= args.min_n_refine_fl
    nr_fl_fit, d_fl_fit, e_fl_fit = nr_fl[mask], d_fl[mask], e_fl[mask]

    nr_rs, d_rs, e_rs = realspace_series(paths)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.axhline(0.125, color="#4d4d4d", linestyle="--", linewidth=1, label="exact 1/8", zorder=1)

    ax.errorbar(1.0 / nr_rs, d_rs, yerr=e_rs, fmt="s", color="#762a83", capsize=3,
                label="real-space windowed fit", markersize=6, zorder=3)
    if len(nr_rs) >= 3:
        d_inf, d_inf_err, c, rchi2 = extrapolate(nr_rs, d_rs, e_rs, 1.0)
        xx = np.linspace(0, (1.0 / nr_rs).max() * 1.1, 100)
        ax.plot(xx, d_inf + c * xx, color="#762a83", alpha=0.5, linewidth=1)
        ax.plot(0, d_inf, marker="*", markersize=16, color="#762a83", zorder=5,
                label=f"real-space $\\Delta_\\infty$={d_inf:.5f}$\\pm${d_inf_err:.5f} (p=1)")

    ax.errorbar(1.0 / nr_fl, d_fl, yerr=e_fl, fmt="o", color="#1b7837", capsize=3,
                label="F_l projection (corrected, no 2l+1)", markersize=6, zorder=3)
    if len(nr_fl_fit) >= 3:
        d_inf, d_inf_err, c, rchi2 = extrapolate(nr_fl_fit, d_fl_fit, e_fl_fit, args.power)
        xx = np.linspace(0, (1.0 / nr_fl).max() * 1.1, 100)
        ax.plot(xx, d_inf + c * xx ** args.power, color="#1b7837", alpha=0.5, linewidth=1)
        ax.plot(0, d_inf, marker="*", markersize=16, color="#1b7837", zorder=5,
                label=f"F_l $\\Delta_\\infty$={d_inf:.5f}$\\pm${d_inf_err:.5f} (p={args.power:g}, "
                      f"n_refine$\\geq${args.min_n_refine_fl}, $\\chi^2$/dof={rchi2:.0f})")

    for nr, x, y in zip(nr_fl.astype(int), 1.0 / nr_fl, d_fl):
        ax.annotate(f"{nr}", (x, y), textcoords="offset points", xytext=(5, -10), fontsize=7, color="#1b7837")

    ax.set_xlabel("1 / n_refine  (continuum $\\rightarrow$ 0)")
    ax.set_ylabel(r"$\Delta_\sigma(a)$")
    ax.set_title("$\\Delta_\\sigma$: F_l-projection vs. real-space window fit")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

    print("\nF_l method:")
    for nr, d, e in zip(nr_fl, d_fl, e_fl):
        print(f"  n_refine={int(nr):3d}  Delta_s={d:.6f} +/- {e:.6f}")
    print("real-space method:")
    for nr, d, e in zip(nr_rs, d_rs, e_rs):
        print(f"  n_refine={int(nr):3d}  Delta={d:.6f} +/- {e:.6f}")


if __name__ == "__main__":
    main()
