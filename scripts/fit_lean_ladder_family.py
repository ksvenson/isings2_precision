#!/usr/bin/env python3
"""Fit a family of continuum-extrapolation models to Delta_l_pair(a) (or
Delta_s(a) for l=1) independently for each harmonic level l=1..l_max of a
lean_harmonic_stats ladder, and report which model each level actually
prefers (by AIC, since some models have more free parameters than others
-- a lower chi2 alone is not a fair comparison across model complexity).

Motivated by the finding (this session, per user question "there should
be a quadratic scaling to the continuum") that l=1 prefers p=1 (linear)
while l>=2 prefer p=2 (quadratic) -- fit_lean_ladder.py only tried p=1/p=2
fixed-power fits. This script extends that with a floated-power fit and
two mixed two-term models, so the "family" isn't limited to the two
powers already known to disagree across l.

Models fit per l (all via scipy.optimize.least_squares, weighted by each
ladder point's own jackknife error from analyze_lean_harmonic_stats.py):
  p1       Delta_inf + c/n                      (2 params, fixed p=1)
  p2       Delta_inf + c/n^2                     (2 params, fixed p=2)
  p3       Delta_inf + c/n^3                     (2 params, fixed p=3)
  pfloat   Delta_inf + c/n^p                     (3 params, p free)
  mix12    Delta_inf + c1/n + c2/n^2             (3 params)
  mix13    Delta_inf + c1/n + c2/n^3             (3 params)

Model selection: AIC = chi2 + 2*k (k = number of free parameters), lowest
AIC wins per level -- penalizes the 3-parameter models for the extra
flexibility rather than letting them always win on raw chi2.

Usage:
  scripts/fit_lean_ladder_family.py <dir_with_*_lean_*.dat files> [--l_max 8] [-o out_prefix]
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
from analyze_lean_harmonic_stats import analyze  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]

EXACT = 0.125


def fit_single_power(n, y, e, p):
    def resid(params):
        d_inf, c = params
        return (d_inf + c * n ** (-p) - y) / e
    sol = least_squares(resid, (0.125, 0.0))
    chi2 = float(np.sum(sol.fun ** 2))
    dof = len(n) - 2
    cov = _cov(sol, chi2, dof)
    d_inf_err = np.sqrt(cov[0, 0]) if cov is not None else float("nan")
    return dict(name=f"p{p:g}", k=2, d_inf=sol.x[0], d_inf_err=d_inf_err,
                params=sol.x, chi2=chi2, dof=dof,
                model=lambda x, s=sol.x, p=p: s[0] + s[1] * x ** (-p))


def fit_floated_power(n, y, e):
    def resid(params):
        d_inf, c, p = params
        return (d_inf + c * n ** (-p) - y) / e
    sol = least_squares(resid, (0.125, 0.05, 1.5), bounds=([0, -np.inf, 0.1], [np.inf, np.inf, 5.0]))
    chi2 = float(np.sum(sol.fun ** 2))
    dof = len(n) - 3
    cov = _cov(sol, chi2, dof)
    d_inf_err = np.sqrt(cov[0, 0]) if cov is not None else float("nan")
    return dict(name="pfloat", k=3, d_inf=sol.x[0], d_inf_err=d_inf_err,
                params=sol.x, chi2=chi2, dof=dof,
                model=lambda x, s=sol.x: s[0] + s[1] * x ** (-s[2])), sol.x[2]


def fit_mix(n, y, e, p2):
    def resid(params):
        d_inf, c1, c2 = params
        return (d_inf + c1 / n + c2 * n ** (-p2) - y) / e
    sol = least_squares(resid, (0.125, 0.05, 0.0))
    chi2 = float(np.sum(sol.fun ** 2))
    dof = len(n) - 3
    cov = _cov(sol, chi2, dof)
    d_inf_err = np.sqrt(cov[0, 0]) if cov is not None else float("nan")
    return dict(name=f"mix1{p2:g}", k=3, d_inf=sol.x[0], d_inf_err=d_inf_err,
                params=sol.x, chi2=chi2, dof=dof,
                model=lambda x, s=sol.x, p2=p2: s[0] + s[1] / x + s[2] * x ** (-p2))


def _cov(sol, chi2, dof):
    if dof <= 0:
        return None
    try:
        return np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
    except np.linalg.LinAlgError:
        return None


def fit_all(n, y, e):
    results = [
        fit_single_power(n, y, e, 1.0),
        fit_single_power(n, y, e, 2.0),
        fit_single_power(n, y, e, 3.0),
        fit_mix(n, y, e, 2.0),
        fit_mix(n, y, e, 3.0),
    ]
    pfloat, p_val = fit_floated_power(n, y, e)
    pfloat["p_value"] = p_val
    results.append(pfloat)
    for r in results:
        r["aic"] = r["chi2"] + 2 * r["k"]
        r["chi2_dof"] = r["chi2"] / r["dof"] if r["dof"] > 0 else float("nan")
    results.sort(key=lambda r: r["aic"])
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out_prefix", default=None)
    args = ap.parse_args()
    out_prefix = args.out_prefix or os.path.join(args.data_dir, "plots", "delta_l_pair_family")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n.setdefault(n, []).append(f)

    n_refines = []
    dp_all, dp_err_all = [], []
    for n in sorted(by_n):
        r = analyze(sorted(by_n[n]), args.l_max)
        n_refines.append(n)
        dp_all.append(r["delta_l_pair"])
        dp_err_all.append(r["delta_l_pair_err"])
    n_refines = np.array(n_refines, dtype=float)
    dp_all = np.array(dp_all)
    dp_err_all = np.array(dp_err_all)

    ncols = 3
    nrows = int(np.ceil(args.l_max / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.6 * nrows), squeeze=False)
    inv_L = 1.0 / n_refines
    xx = np.linspace(1e-9, inv_L.max() * 1.05, 200)

    summary_lines = []
    header = f"{'l':>3} {'best_model':>10} {'Delta_inf':>11} {'err':>9} {'sigma':>7} {'chi2/dof':>9} {'AIC':>8}  {'p (if pfloat)':>14}"
    print(header)
    summary_lines.append(header)

    for l in range(1, args.l_max + 1):
        y = dp_all[:, l - 1]
        yerr = dp_err_all[:, l - 1]
        results = fit_all(n_refines, y, yerr)
        best = results[0]

        for r in results:
            extra = f"  p={r['p_value']:.2f}" if r["name"] == "pfloat" else ""
            sigma = (r["d_inf"] - EXACT) / r["d_inf_err"] if r["d_inf_err"] == r["d_inf_err"] else float("nan")
            line = (f"{l:3d} {r['name']:>10} {r['d_inf']:11.5f} {r['d_inf_err']:9.5f} "
                    f"{sigma:7.2f} {r['chi2_dof']:9.3f} {r['aic']:8.2f}{extra}"
                    + ("  <-- best" if r is best else ""))
            print(line)
            summary_lines.append(line)

        ax = axes[(l - 1) // ncols][(l - 1) % ncols]
        color = COLORS[(l - 1) % len(COLORS)]
        ax.errorbar(inv_L, y, yerr=yerr, fmt="o", ms=4, capsize=2, color=color, zorder=5, label="data")
        for r in results[:3]:
            n_x = xx ** -1  # placeholder unused
            yy = [r["model"](1.0 / iv) for iv in xx]
            style = "-" if r is best else "--"
            lw = 1.8 if r is best else 1.0
            alpha = 1.0 if r is best else 0.55
            ax.plot(xx, yy, style, linewidth=lw, alpha=alpha,
                    label=f"{r['name']} (AIC={r['aic']:.1f})" + (" *best*" if r is best else ""))
        ax.axhline(EXACT, color="gray", linestyle=":", linewidth=1.0)
        ax.set_title(f"l={l}  best={best['name']}  Delta_inf={best['d_inf']:.4f}", fontsize=9)
        ax.set_xlabel(r"$1/n_{\rm refine}$", fontsize=8)
        ax.set_ylabel(r"$\Delta_l^{\rm pair}(a)$", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=6, loc="best")

    for idx in range(args.l_max, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle(r"Model family per level: p1, p2, p3, pfloat, mix(1/n+1/n$^2$), mix(1/n+1/n$^3$) -- best by AIC shown solid")
    fig.tight_layout()
    plot_path = out_prefix + ".png"
    fig.savefig(plot_path, dpi=150)
    print(f"\nwrote {plot_path}", file=sys.stderr)

    txt_path = out_prefix + ".txt"
    with open(txt_path, "w") as f:
        f.write("\n".join(summary_lines) + "\n")
    print(f"wrote {txt_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
