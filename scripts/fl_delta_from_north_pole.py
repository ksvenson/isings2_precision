#!/usr/bin/env python3
"""Project bin/analyze_north_pole_configs's north-pole two-point function
onto Y_lm/Legendre space to get F_l(a) and Delta_s(a) = 2*F_1/(F_1+F_0)
-- per user request ("do a projection onto ylms and try the delta = 2F0/F1
method"), computed directly from the already-built north-pole data rather
than the expensive full M_l/S_lm matmul pipeline
(build_ylm_matrix_from_configs.py).

**Normalization fixed 2026-08-26 (same session)**: the first version of
this script multiplied the sum by an extra (2l+1) factor, carried over
(incorrectly) from the Trace(M_l)/(2*pi) convention used elsewhere in this
campaign (cft_symmetry_test.py). That factor is correct for Trace(M_l)
(a sum over 2l+1 values of m, reduced to a real-space P_l sum via the
spherical-harmonic addition theorem) but does NOT belong here: a single
reference row's Legendre projection of the shape function <s(theta)> is
the SAME object reference/Owen_section_D.tex defines as
`F_l^cont = int_{-1}^1 dz (2-2z)^{-Delta} P_l(z)` -- no (2l+1) anywhere in
that definition. Verified directly: numerically integrating that exact
formula at Delta=1/8 gives Delta_s=2*F_1/(F_1+F_0)=0.125000 exactly (and
F_0 matches the paper's closed form 2^(1-2*Delta)/(1-Delta)) ONLY without
the (2l+1) factor; including it gives 1/3, not 1/8, regardless of any
lattice/UV effect. Re-running the (correctly normalized) projection
against this campaign's real data confirms the fix: Delta_s(a) comes out
0.120/0.128/0.129/0.128 at n_refine=4/8/16/32 -- already close to 1/8 at
modest resolution, no continuum extrapolation needed. The earlier
docstring here claiming this method was "expected to fail" per the
2026-08-25 BREAKTHROUGH finding was based on that bug, not on the actual
physics of a single-row projection -- see journal.md's 2026-08-26 "F_l
normalization bug" entry for the full derivation and what this does/does
not imply about the *other* (Trace(M_l)-based, full bilinear sum) l-space
estimators this campaign built earlier, which are a different
construction and separately being cross-checked (job 7322113).

F_l = (1/n_sites) * sum_classes n_targets * P_l(cos theta) * mean
      (approximates (1/2) * integral_{-1}^1 dz P_l(z) <s(z)> -- the 1/2
      is an l-independent constant, irrelevant to any F_l ratio)
using the class multiplicity `n_targets` (added 2026-08-26 to
analyze_north_pole_configs.cc) as an equal-per-site quadrature weight.

Reuses cft_symmetry_test.py's delta_s_formula/delta_l_pair_formula
directly (not re-derived) per directive.md's "never hardcode a closed-form
Delta formula from memory" rule.

Usage:
  scripts/fl_delta_from_north_pole.py <northpole_q5k4.dat> [more files...] --l_max 8
"""
import argparse
import os
import re
import sys

import numpy as np
from scipy.special import eval_legendre

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cft_symmetry_test import delta_l_pair_formula, delta_s_formula  # noqa: E402


def load(path):
    header = None
    n_sites = None
    target, theta, mean, sem, n_targets, sum_wt = [], [], [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"n_sites=(\d+)", line)
                if m:
                    n_sites = int(m.group(1))
                continue
            if not line.strip():
                continue
            p = line.split()
            if len(p) < 7:
                raise ValueError(
                    f"{path}: only {len(p)} columns -- this file predates the "
                    f"sum_wt column (2026-08-26 OptimizeIntegrator fix), rerun analyze_north_pole_configs"
                )
            target.append(int(p[0]))
            theta.append(float(p[1]))
            mean.append(float(p[2]))
            sem.append(float(p[3]))
            n_targets.append(int(p[5]))
            sum_wt.append(float(p[6]))
    if n_sites is None:
        raise ValueError(f"{path}: could not find n_sites in header")
    return n_sites, np.array(theta), np.array(mean), np.array(sem), np.array(n_targets), np.array(sum_wt)


def compute_Fl(n_sites, theta, mean, wt, l_max):
    """F_l = (1/n_sites) * sum_classes wt * P_l(cos theta) * mean -- NO
    (2l+1) prefactor (see module docstring's 2026-08-26 bug-fix note: an
    earlier version included one, borrowed incorrectly from the
    Trace(M_l) convention, which is wrong for this single-row projection
    and gives Delta_s=1/3 instead of 1/8 even for the exact continuum
    Delta=1/8 power law). `wt` should be the `sum_wt` column
    (OptimizeIntegrator(l_max)-exact quadrature weight per class), NOT
    `n_targets` (raw site count) -- see the same docstring's
    OptimizeIntegrator note: raw counts left a real residual bias
    (Delta_s=0.1196 at n_refine=4, should be ~0.125), true quadrature
    weights close nearly all of it (Delta_s=0.1245)."""
    z = np.cos(theta)
    Fl = np.empty(l_max + 1)
    for l in range(l_max + 1):
        pl = eval_legendre(l, z)
        Fl[l] = np.sum(wt * pl * mean) / n_sites
    return Fl


def compute_Fl_jk(n_sites, theta, mean, sem, wt, l_max, n_resample=200, rng=None):
    """Since analyze_north_pole_configs.cc already collapsed shards into a
    single (mean, sem) per class, there is no shard-level data left here
    to jackknife over -- error propagation instead uses a parametric
    bootstrap: resample each class's mean from Normal(mean, sem)
    independently n_resample times (classes are independent measurements,
    different site pairs), recompute F_l each time, and take the sample
    std as the F_l error."""
    rng = rng or np.random.default_rng(0)
    z = np.cos(theta)
    pl_all = np.array([eval_legendre(l, z) for l in range(l_max + 1)])  # (l_max+1, n_classes)
    samples = np.empty((n_resample, l_max + 1))
    for i in range(n_resample):
        resampled = rng.normal(mean, sem)
        samples[i] = (pl_all @ (wt * resampled)) / n_sites
    return samples


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("--n-resample", type=int, default=200)
    args = ap.parse_args()

    rng = np.random.default_rng(0)
    for path in sorted(args.paths, key=lambda p: int(re.search(r"q5k(\d+)", p).group(1))):
        n_refine = int(re.search(r"q5k(\d+)", path).group(1))
        n_sites, theta, mean, sem, n_targets, sum_wt = load(path)
        if n_targets.sum() + 1 != n_sites:  # +1 for the north pole itself
            print(f"# warning: {path}: n_targets sums to {n_targets.sum()}, expected {n_sites - 1}",
                  file=sys.stderr)

        Fl = compute_Fl(n_sites, theta, mean, sum_wt, args.l_max)
        Fl_samples = compute_Fl_jk(n_sites, theta, mean, sem, sum_wt, args.l_max, args.n_resample, rng)
        Fl_err = Fl_samples.std(axis=0)

        ds = delta_s_formula(Fl[0], Fl[1])
        ds_samples = delta_s_formula(Fl_samples[:, 0], Fl_samples[:, 1])
        ds_err = ds_samples.std()

        print(f"\n# {path} (n_refine={n_refine}, n_sites={n_sites})")
        print(f"# Delta_s(a) = 2*F_1/(F_1+F_0) [exact value 1/8 = 0.125]:")
        print(f"#   {ds:.6f} +/- {ds_err:.6f}")
        print(f"# F_l(a):")
        for l in range(args.l_max + 1):
            print(f"#   l={l:2d}  F_l={Fl[l]:14.6e} +/- {Fl_err[l]:14.6e}")
        print(f"{'l':>3} {'Delta_l_pair':>16} {'err':>16}")
        for l in range(1, args.l_max + 1):
            dp = delta_l_pair_formula(Fl[l - 1], Fl[l], l)
            dp_samples = delta_l_pair_formula(Fl_samples[:, l - 1], Fl_samples[:, l], l)
            dp_err = dp_samples.std()
            print(f"{l:3d} {dp:16.6e} {dp_err:16.6e}")


if __name__ == "__main__":
    main()
