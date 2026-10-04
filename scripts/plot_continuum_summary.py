#!/usr/bin/env python3
"""Plots Delta_sigma(a) (from bin/analyze_north_pole_configs, via
fit_north_pole_ladder.py's fit_delta) and kappa2_r(a) per harmonic level l
(from build_ylm_matrix_from_configs.py's output) vs 1/n_refine, for
whatever ladder points are currently available -- built 2026-08-26 to give
a running view of the continuum behavior of both objective (2) (Delta
extraction) and objective (1) (spherical symmetry) while the rest of the
save_configs_ladder_2026-08-26 production analysis is still running.

Usage:
  scripts/plot_continuum_summary.py <analysis_dir> -o <out_prefix>
    (reads <analysis_dir>/northpole_q5k*.dat and <analysis_dir>/ml_q5k*.txt,
    whichever are present and non-empty)
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fit_north_pole_ladder import fit_delta, load as load_northpole  # noqa: E402


def plot_delta(analysis_dir, out_path, theta_min=0.3, theta_max=2.8, power=1.0):
    paths = sorted(
        glob.glob(os.path.join(analysis_dir, "northpole_q5k*.dat")),
        key=lambda p: int(re.search(r"q5k(\d+)", p).group(1)),
    )
    paths = [p for p in paths if os.path.getsize(p) > 0]
    if not paths:
        print("no north-pole data yet, skipping Delta plot")
        return

    n_refines, deltas, errs = [], [], []
    for p in paths:
        n_refine = int(re.search(r"q5k(\d+)", p).group(1))
        theta, mean, err = load_northpole(p)
        d, de, rchi2, n_used = fit_delta(theta, mean, err, theta_min, theta_max)
        n_refines.append(n_refine)
        deltas.append(d)
        errs.append(de)
        print(f"n_refine={n_refine}: Delta={d:.5f}+/-{de:.5f} chi2/dof={rchi2:.2f}")

    n_refines = np.array(n_refines, dtype=float)
    deltas = np.array(deltas)
    errs = np.array(errs)
    inv_l = 1.0 / n_refines

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(inv_l, deltas, yerr=errs, fmt="o", color="#2166ac", capsize=3, label="measured $\\Delta_\\sigma(a)$")
    ax.axhline(0.125, color="#b2182b", linestyle="--", linewidth=1, label="exact 1/8")

    if len(n_refines) >= 3:
        from scipy.optimize import least_squares

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
        xx = np.linspace(0, inv_l.max() * 1.05, 100)
        ax.plot(xx, sol.x[0] + sol.x[1] * xx ** power, color="#2166ac", linewidth=1, alpha=0.6,
                 label=f"fit: $\\Delta_\\infty$={sol.x[0]:.5f}$\\pm${d_inf_err:.5f}")
        ax.plot(0, sol.x[0], marker="*", markersize=14, color="#2166ac", zorder=5)
        sigma = (sol.x[0] - 0.125) / d_inf_err if d_inf_err > 0 else float("nan")
        print(f"Delta_inf = {sol.x[0]:.5f} +/- {d_inf_err:.5f} ({sigma:+.2f} sigma from 1/8)")

    for nr, x, y in zip(n_refines.astype(int), inv_l, deltas):
        ax.annotate(f"n_refine={nr}", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)

    ax.set_xlabel("1 / n_refine  (continuum $\\rightarrow$ 0)")
    ax.set_ylabel(r"$\Delta_\sigma(a)$")
    ax.set_title("Real-space $\\Delta_\\sigma$ continuum extrapolation (north-pole method)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


# ml_q5k*.txt column indices (0-based, matches build_ylm_matrix_from_configs.py's
# print format: l dim R_l R_l_err chi2 dof chi2/dof c_hat kappa2_r err kappa3_r err kappa4_r err)
_KAPPA_COLUMNS = {
    # order -> (value_col, err_col, is_reduced_cumulant)
    # kappa1 has no direct column -- c_hat (GLS common-diagonal estimate)
    # is used as a proxy: c_hat == Tr(M_l)/dim == kappa1 exactly when the
    # diagonal is constant (i.e. when kappa2_r etc are small), which is
    # the same regime this whole diagnostic is meant to probe.
    1: (7, None, False),
    2: (8, 9, True),
    3: (10, 11, True),
    4: (12, 13, True),
}


def plot_kappa_order(analysis_dir, out_path, order, l_max_show=8, linear=False):
    value_col, err_col, is_reduced = _KAPPA_COLUMNS[order]
    paths = sorted(
        glob.glob(os.path.join(analysis_dir, "ml_q5k*.txt")),
        key=lambda p: int(re.search(r"q5k(\d+)", p).group(1)),
    )
    paths = [p for p in paths if os.path.getsize(p) > 0]
    if not paths:
        print(f"no M_l data yet, skipping kappa{order} plot")
        return

    data = {}  # l -> list of (n_refine, value, err)
    for p in paths:
        n_refine = int(re.search(r"q5k(\d+)", p).group(1))
        with open(p) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if parts[0] == "l":  # table header row (no leading '#')
                    continue
                l = int(parts[0])
                val_str = parts[value_col]
                if val_str == "--":
                    continue
                val = float(val_str)
                err = float(parts[err_col]) if err_col is not None and parts[err_col] != "--" else 0.0
                data.setdefault(l, []).append((n_refine, val, err))

    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")
    ls = sorted(l for l in data if l <= l_max_show)
    for i, l in enumerate(ls):
        pts = sorted(data[l])
        nr = np.array([p[0] for p in pts], dtype=float)
        vals = np.array([p[1] for p in pts])
        errs = np.array([p[2] for p in pts])
        color = cmap(i / max(1, len(ls) - 1))
        y = vals if (linear or not is_reduced) else np.abs(vals)
        ax.errorbar(1.0 / nr, y, yerr=errs, fmt="o-", color=color, capsize=3, label=f"l={l}", markersize=5)

    if linear or not is_reduced:
        ax.axhline(0, color="#4d4d4d", linewidth=0.8, zorder=1)

    if order == 1:
        ax.set_ylabel(r"$c\_hat(l) \approx \kappa_1(l) = \mathrm{Tr}(M_l)/\mathrm{dim}$")
        title = "M_l overall scale (kappa1 proxy via GLS c_hat) vs. resolution"
    else:
        ylabel = rf"$\kappa_{{{order},r}}(l)$"
        if not linear:
            ylabel = rf"$|\kappa_{{{order},r}}(l)|$"
        ax.set_ylabel(ylabel + "  (spherical-symmetry breaking)")
        title = f"SO(3)-breaking spectral cumulant (order {order}) vs. resolution, by harmonic level"

    ax.set_xlabel("1 / n_refine  (continuum $\\rightarrow$ 0)")
    ax.set_title(title)
    if not linear:
        ax.set_yscale("log")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("analysis_dir")
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--linear", action="store_true", help="linear y-axis for the kappa2_r plot (default: log, |kappa2_r|)")
    args = ap.parse_args()

    plot_delta(args.analysis_dir, f"{args.out_prefix}_delta_vs_invL.png")
    suffix = "_linear" if args.linear else ""
    for order in (1, 2, 3):
        plot_kappa_order(args.analysis_dir, f"{args.out_prefix}_kappa{order}_vs_invL{suffix}.png",
                          order, linear=args.linear)


if __name__ == "__main__":
    main()
