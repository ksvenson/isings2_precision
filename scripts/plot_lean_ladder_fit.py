#!/usr/bin/env python3
"""Plot Delta_s(a) vs 1/n_refine for a lean_harmonic_stats ladder, overlaid
with the p=1 and p=2 continuum-extrapolation fits from fit_lean_ladder.py
(reuses its fit_power(), no reimplementation).

Usage:
  scripts/plot_lean_ladder_fit.py <dir_with_*_lean_*.dat files> [--l_max 8] [-o out.png]
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
from analyze_lean_harmonic_stats import analyze  # noqa: E402
from fit_lean_ladder import fit_power  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.data_dir, "plots", "delta_s_continuum_fit.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n.setdefault(n, []).append(f)

    n_refines, deltas, errs = [], [], []
    for n in sorted(by_n):
        r = analyze(sorted(by_n[n]), args.l_max)
        n_refines.append(n)
        deltas.append(r["delta_s"])
        errs.append(r["delta_s_err"])
    n_refines = np.array(n_refines, dtype=float)
    deltas = np.array(deltas)
    errs = np.array(errs)

    d1, e1, c1, chi1 = fit_power(n_refines, deltas, errs, 1.0)
    d2, e2, c2, chi2v = fit_power(n_refines, deltas, errs, 2.0)

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    inv_L = 1.0 / n_refines
    ax.errorbar(inv_L, deltas, yerr=errs, fmt="o", color="#0072B2", ms=5, capsize=3,
                label="lean_harmonic_stats data")
    xx = np.linspace(0, inv_L.max() * 1.05, 200)
    ax.plot(xx, d1 + c1 * xx, "-", color="#D55E00",
            label=rf"$p{{=}}1$: $\Delta_\infty$={d1:.5f}({e1*1e5:.0f}), $\chi^2$/dof={chi1:.2f}")
    ax.plot(xx, d2 + c2 * xx ** 2, "--", color="#009E73",
            label=rf"$p{{=}}2$: $\Delta_\infty$={d2:.5f}({e2*1e5:.0f}), $\chi^2$/dof={chi2v:.2f}")
    ax.axhline(0.125, color="gray", linestyle=":", linewidth=1.2, label="exact 1/8")
    ax.set_xlabel(r"$1/n_{\rm refine}$")
    ax.set_ylabel(r"$\Delta_s(a)$")
    ax.set_title("Continuum extrapolation")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
