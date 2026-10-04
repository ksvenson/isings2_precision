#!/usr/bin/env python3
"""Plot Delta_sigma(a) vs 1/n_refine (~1/L) for one or more real-space
ladder runs (fit_real_space_ladder.py's analyze_point, reused directly),
overlaid with the Delta_inf + c*n_refine^-p extrapolation fit.

Usage:
  scripts/plot_real_space_ladder.py <run_root> [<run_root> ...] \\
    [--power 1] [--theta-min 0.3] [--theta-max 2.8] [--out PATH]
"""
import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from fit_real_space_ladder import analyze_point


def ladder(run_root, theta_min, theta_max):
    point_dirs = sorted(
        glob.glob(os.path.join(run_root, "q5k*")),
        key=lambda p: int(re.search(r"q5k(\d+)", p).group(1)),
    )
    n_refines, deltas, errs = [], [], []
    for pd in point_dirs:
        n_refine = int(re.search(r"q5k(\d+)", pd).group(1))
        _, _, delta_bar, delta_err = analyze_point(pd, theta_min, theta_max)
        n_refines.append(n_refine)
        deltas.append(delta_bar)
        errs.append(delta_err)
    return np.array(n_refines, dtype=float), np.array(deltas), np.array(errs)


def fit(n_refines, deltas, errs, power):
    def resid(p):
        d_inf, c = p
        return (d_inf + c * n_refines ** (-power) - deltas) / errs

    sol = least_squares(resid, (0.125, 0.0))
    chi2 = np.sum(sol.fun ** 2)
    dof = len(n_refines) - 2
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
    return sol.x[0], np.sqrt(cov[0, 0]), sol.x[1], chi2 / dof


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_roots", nargs="+")
    ap.add_argument("--power", type=float, default=1.0)
    ap.add_argument("--theta-min", type=float, default=0.3)
    ap.add_argument("--theta-max", type=float, default=2.8)
    ap.add_argument("--out", default="campaign_runs/real_space_ladder_2026-08-25/plots/delta_vs_invL.png")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    colors = plt.cm.tab10.colors

    for i, run_root in enumerate(args.run_roots):
        label = os.path.basename(os.path.normpath(run_root))
        n_refines, deltas, errs = ladder(run_root, args.theta_min, args.theta_max)
        inv_l = 1.0 / n_refines
        d_inf, d_inf_err, c, chi2dof = fit(n_refines, deltas, errs, args.power)

        color = colors[i % len(colors)]
        ax.errorbar(inv_l, deltas, yerr=errs, fmt="o", color=color, label=f"{label} data")

        x_fit = np.linspace(0, inv_l.max() * 1.05, 200)
        y_fit = d_inf + c * x_fit ** args.power
        ax.plot(x_fit, y_fit, "-", color=color, alpha=0.7,
                 label=f"{label} fit: Delta_inf={d_inf:.5f}+/-{d_inf_err:.5f} "
                       f"({(d_inf-0.125)/d_inf_err:+.2f}sig, chi2/dof={chi2dof:.2f})")

    ax.axhline(0.125, color="black", linestyle="--", linewidth=1, label="exact 1/8")
    ax.set_xlabel(r"$1/n_{\rm refine} \approx 1/L$")
    ax.set_ylabel(r"$\Delta_\sigma(a)$")
    ax.set_title(rf"Real-space $\Delta_\sigma$ continuum extrapolation ($\Delta = \Delta_\infty + c\, n_{{\rm refine}}^{{-{args.power:g}}}$)")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
