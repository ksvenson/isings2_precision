#!/usr/bin/env python3
"""Extract Delta from windowed F_0/F_1, computed with a *consistent*
windowed integration scheme on both the data and the theory side --
per user question 2026-08-25 ("are you able to get delta from F0 and F1
with the correct integration scheme?").

Why the earlier attempt failed (journal.md's "windowed-Legendre-integral
attempt... FAILS"/"SETTLES the question" entries): that attempt computed
a windowed F_l from data, then tried to invert it via the *unwindowed*
closed-form recursion F_l^cont(Delta) (Owen_section_D's
(l-1+Delta)/(l+1-Delta) recursion, and the shortcut identity
Delta = 2*F_1/(F_1+F_0)) -- both are only valid for a FULL-SPHERE
integral of the pure power law. Using them on a windowed F_l is a
mismatched comparison, not a fundamental failure of F_0/F_1 as
Delta-carrying quantities.

This script instead computes F_0^model(Delta), F_1^model(Delta) by
DIRECT NUMERICAL INTEGRATION of A*(2-2z)^(-Delta)*P_l(z) over the exact
same [z_min, z_max] window used for the data sum -- no recursion
anywhere -- then fits Delta by matching the data ratio 2*F_1/(F_1+F_0)
(amplitude-independent) to the model ratio at the same window. Since the
window is theta in [theta_min, theta_max], i.e. z in
[cos(theta_max), cos(theta_min)], and the real_space_2pt_test bins are
uniform in z, F_l_data = sum_bins mean(bin) * P_l(z_center) * dz exactly
(sin(theta)dtheta = -dz, so no extra Jacobian is needed once z is the
integration variable).

Usage:
  scripts/windowed_f0_f1_delta.py <real_space_2pt_dat_file>
    [--theta-min 0.3] [--theta-max 2.8]
"""
import argparse

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append((float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])))
    z_lo, z_hi, z_center, theta, mean = map(np.array, zip(*rows))
    return z_lo, z_hi, z_center, theta, mean


def p1(z):
    return z


def data_Fl(z_center, mean, theta, theta_min, theta_max, dz, l):
    mask = (theta >= theta_min) & (theta <= theta_max)
    pl = np.ones_like(z_center) if l == 0 else p1(z_center)
    return np.sum(mean[mask] * pl[mask] * dz)


def model_Fl(Delta, z_min, z_max, l):
    pl = (lambda z: 1.0) if l == 0 else p1

    def integrand(z):
        return (2.0 - 2.0 * z) ** (-Delta) * pl(z)

    val, _ = quad(integrand, z_min, z_max, limit=200)
    return val


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("--theta-min", type=float, default=0.3)
    ap.add_argument("--theta-max", type=float, default=2.8)
    args = ap.parse_args()

    z_lo, z_hi, z_center, theta, mean = load(args.path)
    n_bins = len(z_center)
    dz = 2.0 / n_bins

    F0_data = data_Fl(z_center, mean, theta, args.theta_min, args.theta_max, dz, 0)
    F1_data = data_Fl(z_center, mean, theta, args.theta_min, args.theta_max, dz, 1)
    ratio_data = 2.0 * F1_data / (F1_data + F0_data)
    print(f"# {args.path}")
    print(f"windowed data:  F0={F0_data:.6e}  F1={F1_data:.6e}  2F1/(F1+F0)={ratio_data:.5f}")

    z_min = np.cos(args.theta_max)
    z_max = np.cos(args.theta_min)

    def ratio_model(Delta):
        F0m = model_Fl(Delta, z_min, z_max, 0)
        F1m = model_Fl(Delta, z_min, z_max, 1)
        return 2.0 * F1m / (F1m + F0m)

    # ratio_model(Delta) is monotonic in Delta over (0,1); solve for match
    def resid(Delta):
        return ratio_model(Delta) - ratio_data

    lo, hi = 1e-3, 0.99
    r_lo, r_hi = resid(lo), resid(hi)
    print(f"model ratio at Delta={lo}: {ratio_model(lo):.5f}   at Delta={hi}: {ratio_model(hi):.5f}")
    if r_lo * r_hi > 0:
        print("no sign change in [1e-3, 0.99] -- cannot bracket a root; data ratio outside model's achievable range")
        return

    Delta_sol = brentq(resid, lo, hi, xtol=1e-6)
    print(f"\nwindowed F0/F1 solution: Delta = {Delta_sol:.5f}  (exact CFT value: 0.125)")


if __name__ == "__main__":
    main()
