#!/usr/bin/env python3
"""Diagonal-degeneracy spherical-symmetry check from bin/lean_harmonic_stats
data, reusing symmetry_test.py's diag_chi2_const (GLS chi-square of the
per-m diagonal against a constant) instead of the fuller R_l/kappa2_r
(which need the full off-diagonal M_l matrix -- NOT available from this
driver, which only accumulates diagonal |S_lm|^2 per m, see
src/lean_harmonic_stats.cc's header).

Why this is still a real, useful symmetry check despite missing the
off-diagonal terms: under exact SO(3) symmetry, <|S_lm|^2> must be equal
for every m at fixed l (independent of m) -- a rotation mixes different m
into each other and the ensemble average is rotation-invariant, so any
m-dependence in the diagonal itself is direct evidence of symmetry
breaking, entirely separate from the off-diagonal question. This is the
diagonal half of the same gate PLAN.md step 6 defines (chi2_l), just
computed from cheaper, already-collected data that goes all the way to
n_refine=128 -- unlike the full M_l/kappa2_r path (save_configs_ladder_
2026-08-26), which is only analyzed up to n_refine=32 so far (blocked on
a login-node compute bottleneck, see CLAUDE.md/journal.md 2026-08-26).
Does NOT replace kappa2_r/R_l -- a diagonal that passes this test can
still have nonzero off-diagonal terms breaking symmetry in a way this
check is blind to; treat this as a necessary-but-not-sufficient
complement, not a substitute.

Usage:
  scripts/symmetry_check_lean.py <dir_with_*_lean_*.dat files> [--l_max 8] [-o out_prefix]
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
from analyze_lean_harmonic_stats import load, combine  # noqa: E402
from symmetry_test import diag_chi2_const  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def diag_degeneracy(paths, l_max, min_blocks_factor=5):
    header, lm_list, blocks, _offdiag_list = combine(paths)
    n_blocks = len(blocks)
    l_of_k = np.array([lm[0] for lm in lm_list])
    total_n = sum(b["n"] for b in blocks)
    abs2_stack = np.array([b["abs2"] for b in blocks])  # (n_blocks, n_lm)
    n_stack = np.array([b["n"] for b in blocks])
    total_abs2 = abs2_stack.sum(axis=0)

    results = []
    for l in range(l_max + 1):
        mask = l_of_k == l
        dim = mask.sum()
        d_full = total_abs2[mask] / total_n  # per-m mean |S_lm|^2, diagonal
        d_jk = np.empty((n_blocks, dim))
        for k in range(n_blocks):
            loo_n = total_n - n_stack[k]
            d_jk[k] = (total_abs2[mask] - abs2_stack[k][mask]) / loo_n
        d_bar = d_jk.mean(axis=0)
        cov = (n_blocks - 1) / n_blocks * (d_jk - d_bar).T @ (d_jk - d_bar)

        if l == 0:
            results.append(dict(l=l, dim=dim, chi2=0.0, dof=0, c_hat=d_full[0], spread=0.0, spread_err=0.0))
            continue

        # jackknife error of the spread statistic (std/mean of the diagonal)
        spread = float(np.std(d_full) / np.mean(d_full)) if np.mean(d_full) != 0 else float("nan")
        spread_jk = np.array([
            np.std(d_jk[k]) / np.mean(d_jk[k]) if np.mean(d_jk[k]) != 0 else np.nan
            for k in range(n_blocks)
        ])
        spread_bar = np.nanmean(spread_jk)
        spread_err = float(np.sqrt((n_blocks - 1) / n_blocks * np.nansum((spread_jk - spread_bar) ** 2)))

        if n_blocks >= min_blocks_factor * dim:
            chi2, dof, c_hat = diag_chi2_const(d_full, cov)
            # jackknife error of chi2/dof: re-evaluate the GLS chi2 on each
            # leave-one-block-out diagonal sample against the SAME (full-
            # statistics) covariance matrix -- a fixed-cov jackknife of the
            # chi2 statistic's sensitivity to which block is omitted, not a
            # full double-jackknife of the covariance itself (that would need
            # re-estimating cov per leave-one-out sample from only n_blocks-1
            # points, too noisy at this dim/n_blocks ratio).
            chi2_dof_jk = np.array([
                diag_chi2_const(d_jk[k], cov)[0] / dof for k in range(n_blocks)
            ])
            chi2_dof_bar = chi2_dof_jk.mean()
            chi2_dof_err = float(np.sqrt((n_blocks - 1) / n_blocks * np.sum((chi2_dof_jk - chi2_dof_bar) ** 2)))
        else:
            chi2, dof, c_hat, chi2_dof_err = None, dim - 1, None, None
        results.append(dict(l=l, dim=dim, chi2=chi2, dof=dof, c_hat=c_hat, spread=spread,
                             spread_err=spread_err, chi2_dof_err=chi2_dof_err))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out_prefix", default=None)
    args = ap.parse_args()
    out_prefix = args.out_prefix or os.path.join(args.data_dir, "plots", "diag_symmetry")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    files = glob.glob(os.path.join(args.data_dir, "q5k*_lean_*.dat"))
    if not files:
        raise SystemExit(f"no *_lean_*.dat files found in {args.data_dir}")
    by_n = {}
    for f in files:
        n = int(re.search(r"q5k(\d+)_lean", f).group(1))
        by_n.setdefault(n, []).append(f)

    n_refines = sorted(by_n)
    chi2_dof = np.full((len(n_refines), args.l_max + 1), np.nan)
    chi2_dof_err = np.full((len(n_refines), args.l_max + 1), np.nan)
    spread = np.full((len(n_refines), args.l_max + 1), np.nan)
    spread_err = np.full((len(n_refines), args.l_max + 1), np.nan)

    print(f"{'n_refine':>9}  " + "  ".join(f"l={l}(chi2/dof)" for l in range(1, args.l_max + 1)))
    for i, n in enumerate(n_refines):
        res = diag_degeneracy(sorted(by_n[n]), args.l_max)
        row = []
        for r in res:
            l = r["l"]
            if l == 0:
                continue
            spread[i, l] = r["spread"]
            spread_err[i, l] = r["spread_err"]
            if r["chi2"] is not None and r["dof"] > 0:
                chi2_dof[i, l] = r["chi2"] / r["dof"]
                chi2_dof_err[i, l] = r["chi2_dof_err"]
                row.append(f"{chi2_dof[i, l]:.3f}+/-{chi2_dof_err[i, l]:.3f}")
            else:
                row.append(f"{'--':>10}")
        print(f"{n:>9}  " + "  ".join(row))

    n_refines_arr = np.array(n_refines, dtype=float)
    inv_L = 1.0 / n_refines_arr

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for l in range(1, args.l_max + 1):
        color = COLORS[(l - 1) % len(COLORS)]
        axes[0].errorbar(inv_L, chi2_dof[:, l], yerr=chi2_dof_err[:, l], fmt="o-", ms=4, capsize=2,
                          color=color, label=f"l={l}")
        axes[1].errorbar(inv_L, spread[:, l], yerr=spread_err[:, l], fmt="o-", ms=4, capsize=2,
                          color=color, label=f"l={l}")
    axes[0].axhline(1.0, color="gray", linestyle=":", linewidth=1.0, label="chi2/dof=1 (consistent)")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$1/n_{\rm refine}$")
    axes[0].set_ylabel(r"$\chi^2/{\rm dof}$ (diagonal vs constant)")
    axes[0].set_title("Diagonal m-degeneracy test (GLS chi2)")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$1/n_{\rm refine}$")
    axes[1].set_ylabel(r"$\mathrm{std}(|S_{lm}|^2)/\mathrm{mean}$")
    axes[1].set_title("Diagonal spread (basis-dependent, no jackknife)")
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle(r"Spherical-symmetry diagonal check from lean\_harmonic\_stats (l=1..%d, up to $n_{\rm refine}=%d$)"
                 % (args.l_max, int(n_refines_arr.max())))
    fig.tight_layout()
    out = out_prefix + ".png"
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
