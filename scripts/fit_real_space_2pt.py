#!/usr/bin/env python3
"""Fit the real-space two-point function (bin/real_space_2pt_test output)
directly to A*(2-2*cos(theta))**(-Delta) (single operator) or a two-term sum
(two operators), entirely independent of the spherical-harmonic/Legendre-
recursion machinery (cft_symmetry_test.py, two_operator_fit.py) -- built
2026-08-25 per user direction to cross-check the F_l-based Delta extraction
against a direct real-space fit.

**Updated 2026-08-26** for real_space_2pt_test.cc's rewrite that dropped
z=cos(theta) histogram binning entirely in favor of exact icosahedral
equivalence classes (see that file's header comment) -- one row per class
now, format `ref_site target_site theta mean err n orbit_size` (no more
z_lo/z_hi/bin columns). Older .dat files from before the rewrite are NOT
compatible with this loader.

Usage:
  scripts/fit_real_space_2pt.py <real_space_2pt_dat_file> [--theta-min 0.3]
    [--theta-max 2.8] [--two-op]
"""
import argparse
import numpy as np
from scipy.optimize import least_squares


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append((float(p[2]), float(p[3]), float(p[4]), int(p[5])))
    theta, mean, err, n = map(np.array, zip(*rows))
    return theta, mean, err, n


def model1(theta, A, Delta):
    return A * (2 - 2 * np.cos(theta)) ** (-Delta)


def model2(theta, A, D1, B, D2):
    return A * (2 - 2 * np.cos(theta)) ** (-D1) + B * (2 - 2 * np.cos(theta)) ** (-D2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("--theta-min", type=float, default=0.2)
    ap.add_argument("--theta-max", type=float, default=2.9)
    ap.add_argument("--two-op", action="store_true")
    args = ap.parse_args()

    theta, mean, err, n = load(args.path)
    mask = (theta >= args.theta_min) & (theta <= args.theta_max)
    th, y, ye = theta[mask], mean[mask], err[mask]
    print(f"# {args.path}: using {mask.sum()}/{len(theta)} classes, theta in [{args.theta_min},{args.theta_max}]")

    def resid1(p):
        A, D = p
        return (model1(th, A, D) - y) / ye

    sol1 = least_squares(resid1, (y[len(y) // 2], 0.15), bounds=([0, 0.001], [np.inf, 3.0]))
    chi2_1 = np.sum(sol1.fun ** 2)
    dof1 = len(th) - 2
    print(f"single-operator fit: A={sol1.x[0]:.5e} Delta={sol1.x[1]:.5f} chi2/dof={chi2_1/dof1:.3f}")
    try:
        cov1 = np.linalg.inv(sol1.jac.T @ sol1.jac) * (chi2_1 / dof1)
        print(f"  errors: A_err={np.sqrt(cov1[0,0]):.2e} Delta_err={np.sqrt(cov1[1,1]):.5f}")
    except np.linalg.LinAlgError:
        pass

    if args.two_op:
        def resid2(p):
            A, D1, B, D2 = p
            return (model2(th, A, D1, B, D2) - y) / ye

        best = None
        for g in [(y[0] * 0.5, 0.125, y[0] * 0.5, 0.6), (y[0] * 0.8, 0.1, y[0] * 0.2, 1.0),
                  (y[0] * 0.5, 0.3, y[0] * 0.5, 0.05)]:
            sol = least_squares(resid2, g, bounds=([0, 0.001, -np.inf, 0.001], [np.inf, 3.0, np.inf, 3.0]))
            chi2 = np.sum(sol.fun ** 2)
            if best is None or chi2 < best[0]:
                best = (chi2, sol)
        chi2_2, sol2 = best
        dof2 = len(th) - 4
        A, D1, B, D2 = sol2.x
        print(f"two-operator fit: A={A:.5e} D1={D1:.5f} B={B:.5e} D2={D2:.5f} chi2/dof={chi2_2/dof2:.3f}")
        try:
            cov2 = np.linalg.inv(sol2.jac.T @ sol2.jac) * (chi2_2 / dof2)
            errs = np.sqrt(np.diag(cov2))
            print(f"  errors: A={errs[0]:.2e} D1={errs[1]:.5f} B={errs[2]:.2e} D2={errs[3]:.5f}")
        except np.linalg.LinAlgError:
            pass


if __name__ == "__main__":
    main()
