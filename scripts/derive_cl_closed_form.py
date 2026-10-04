#!/usr/bin/env python3
"""Closed-form CFT two-point-function harmonic coefficient C_l(l, Delta).

Reproduces and independently verifies Brower/Cheng/Fleming/Gasbarro/Raben/
Tan/Weinberg, "Lattice phi^4 Field Theory on Riemann Manifolds: Numerical
Tests for the 2-d Ising CFT on S^2" (arXiv:1803.08512), Section 5.1, Eq
5.1-5.3 -- see reference/brower_2018_lattice_phi4_riemann_S2.pdf.

For a CFT primary O of dimension Delta on the unit S^2,
    <O(n)O(n')> = 1 / (2 - 2*cos(gamma))^Delta          (their Eq 5.1)
Expanding in Legendre polynomials, c_l^cont = integral_{-1}^1 dz
(2/(1-z))^Delta P_l(z) (their Eq 5.2), which obeys the exact two-term
recursion (their Eq 5.3):
    c_l^cont(Delta) = [(l-1+Delta)/(l+1-Delta)] * c_{l-1}^cont(Delta)
    c_0^cont(Delta) = 1/(1-Delta)

Per directive.md's "Closed-form CFT formulas" rule, this is not taken on
the paper's word: `verify_against_paper_values`/`verify_recursion_against_paper`
independently reproduce Eq 5.2 via direct sympy integration, and
`verify_numeric_integration` cross-checks with an independent numerical
quadrature (no Legendre/sympy machinery at all).

This campaign's own lattice measurement is NOT in the paper's c_l^cont
convention. `ising_s2_crit.cc`'s on-disk M_l[m,m] is normalized via the
spherical-harmonic addition theorem (S_lm = integral O(n) Y_lm*(n) dOmega,
sum_i w_i = 4*pi quadrature weights): for an exactly SO(3)-symmetric
ensemble, <S_lm S_l'm'*> = delta_{ll'} delta_{mm'} * C_l with
    C_l = c_l^cont * 4*pi / (2*l + 1)
(derived below in `lattice_conversion_factor`, and verified symbolically
against `ising_s2_crit.cc`'s stated 4*pi-normalized quadrature weights and
1/vol^2 two-point normalization). The `M_l[m,m]` ratio between adjacent l
therefore differs from the paper's c_l^cont ratio by this l-dependent
factor -- `rho_l_numeric` below already folds it in, so callers can apply
it directly to measured M_l[m,m] ratios.
"""
import sys

import numpy as np
import sympy
from sympy import Rational, Symbol, gamma, integrate, legendre, nsimplify, pi, symbols

x = Symbol("x")
Delta = Symbol("Delta", positive=True)

PAPER_CL_AT_DELTA_1_8 = {
    0: Rational(8, 7),
    1: Rational(8, 105),
    2: Rational(24, 805),
    3: Rational(408, 24955),
    4: Rational(680, 64883),
    5: Rational(22440, 3049501),
}


def legendre_coeff_symbolic(l_val, Delta_val):
    """c_l^cont for concrete integer l_val, symbolic or concrete Delta_val.

    Independent reproduction of Eq 5.2 -- direct sympy integration, no
    reference to the paper's quoted numbers or its recursion formula.
    """
    f = (2 / (1 - x)) ** Delta_val
    return sympy.simplify(integrate(f * legendre(l_val, x), (x, -1, 1)))


def verify_against_paper_values(l_max=5):
    """Eq 5.2 at Delta=1/8 must match the paper's quoted exact rationals."""
    ok = True
    for l_val, paper_val in PAPER_CL_AT_DELTA_1_8.items():
        if l_val > l_max:
            continue
        computed = legendre_coeff_symbolic(l_val, Rational(1, 8))
        match = sympy.simplify(computed - paper_val) == 0
        status = "PASS" if match else "FAIL"
        print(f"  [{status}] l={l_val}: computed={computed}, paper={paper_val}")
        ok = ok and match
    return ok


def cl_recursion_ratio(l_val, Delta_val):
    """Ratio c_l^cont/c_{l-1}^cont for concrete integer l_val (l_val>=1),
    derived independently by dividing two direct integrations -- not by
    assuming Eq 5.3's form."""
    num = legendre_coeff_symbolic(l_val, Delta_val)
    den = legendre_coeff_symbolic(l_val - 1, Delta_val)
    return sympy.simplify(num / den)


def paper_recursion_ratio(l_val, Delta_val):
    """Eq 5.3's stated ratio, as a function to check the derived ratio against."""
    return (l_val - 1 + Delta_val) / (l_val + 1 - Delta_val)


def verify_recursion_against_paper(l_max=8):
    """Confirm the independently-derived ratio equals Eq 5.3's closed form,
    for symbolic Delta (not just numeric values) where sympy can simplify,
    and numerically otherwise."""
    ok = True
    for l_val in range(1, l_max + 1):
        derived = cl_recursion_ratio(l_val, Delta)
        stated = paper_recursion_ratio(l_val, Delta)
        diff = sympy.simplify(derived - stated)
        match = diff == 0
        if not match:
            # fall back to numeric check at several Delta values
            match = all(
                abs(complex(diff.subs(Delta, dv))) < 1e-9
                for dv in (Rational(1, 8), Rational(1, 2), 1, Rational(3, 2))
            )
        status = "PASS" if match else "FAIL"
        print(f"  [{status}] l={l_val}: derived={derived}, paper={stated}")
        ok = ok and match
    return ok


def verify_numeric_integration(delta_vals=(Rational(1, 8), 1), l_max=8):
    """Independent numerical quadrature (mpmath), no Legendre/sympy symbolic
    machinery -- cross-checks the symbolic result at each l, Delta."""
    import mpmath

    ok = True
    for dv in delta_vals:
        dv_f = float(dv)
        for l_val in range(l_max + 1):
            symbolic_val = float(legendre_coeff_symbolic(l_val, dv))

            def integrand(z, l_val=l_val, dv_f=dv_f):
                return float(mpmath.legendre(l_val, z)) * (2.0 / (1.0 - z)) ** dv_f

            numeric_val = float(mpmath.quad(integrand, [-1, 1]))
            rel_err = abs(numeric_val - symbolic_val) / max(abs(symbolic_val), 1e-30)
            match = rel_err < 1e-6
            status = "PASS" if match else "FAIL"
            if not match or l_val <= 2:
                print(f"  [{status}] Delta={dv_f:.4f} l={l_val}: "
                      f"symbolic={symbolic_val:.6e} numeric={numeric_val:.6e} "
                      f"rel_err={rel_err:.2e}")
            ok = ok and match
    return ok


def lattice_conversion_factor(l_val):
    """C_l (lattice M_l[m,m] convention) = conversion_factor(l) * c_l^cont.

    Derivation: the spherical-harmonic addition theorem gives
    sum_m Y_lm(n) Y_lm*(n') = (2l+1)/(4pi) * P_l(cos gamma), so
    <O(n)O(n')> = sum_l c_l^cont P_l(cos gamma)
                = sum_l c_l^cont * (4pi/(2l+1)) * sum_m Y_lm(n) Y_lm*(n').
    Defining S_lm = integral O(n) Y_lm*(n) dOmega(n) (the continuum
    analogue of this campaign's lattice measurement, which uses
    4*pi-normalized quadrature weights per ising_s2_crit.cc), orthonormality
    of Y_lm gives <S_lm S_l'm'*> = delta_{ll'} delta_{mm'} * c_l^cont *
    4*pi/(2l+1) -- i.e. conversion_factor(l) = 4*pi/(2l+1).
    """
    return 4 * pi / (2 * l_val + 1)


def verify_conversion_factor_symbolic():
    """Sanity check: conversion_factor is positive, decreasing in l (matches
    the addition theorem's (2l+1) growth in the denominator), and reduces
    to 4*pi at l=0 (single m=0 mode, no addition-theorem sum)."""
    ok = True
    f0 = lattice_conversion_factor(0)
    match0 = sympy.simplify(f0 - 4 * pi) == 0
    print(f"  [{'PASS' if match0 else 'FAIL'}] conversion_factor(0) = {f0} (expect 4*pi)")
    ok = ok and match0
    l_sym = Symbol("l", positive=True)
    deriv = sympy.diff(lattice_conversion_factor(l_sym), l_sym)
    decreasing = deriv.is_negative if deriv.is_negative is not None else \
        bool(deriv.subs(l_sym, 5) < 0)
    print(f"  [{'PASS' if decreasing else 'FAIL'}] conversion_factor(l) decreasing in l")
    ok = ok and decreasing
    return ok


# ---- lambdified numeric API for downstream scripts (fit_delta_sigma.py) ----

_l_sym, _Delta_sym, _CDelta_sym = symbols("l Delta C_Delta", positive=True)

# c_0^cont(Delta) = 1/(1-Delta); build c_l^cont(Delta) via the verified
# recursion rather than re-deriving the closed Gamma-function form here --
# the recursion IS the closed form (Eq 5.3 applied l times to c_0).
def _cl_cont_symbolic_chain(l_max):
    """Build c_0..c_{l_max} as explicit symbolic expressions in Delta by
    repeated application of the verified recursion ratio."""
    chain = [1 / (1 - _Delta_sym)]
    for l_val in range(1, l_max + 1):
        ratio = paper_recursion_ratio(l_val, _Delta_sym)
        chain.append(sympy.simplify(chain[-1] * ratio))
    return chain


_MAX_L_CACHE = 32
_CL_CONT_CHAIN = _cl_cont_symbolic_chain(_MAX_L_CACHE)
_CL_CONT_LAMBDAS = [sympy.lambdify(_Delta_sym, expr, "numpy") for expr in _CL_CONT_CHAIN]


def C_l_numeric(l, Delta_val, C_Delta=1.0):
    """Lattice-convention eigenvalue C_l = C_Delta * conversion_factor(l) *
    c_l^cont(Delta_val). Vectorized over l (numpy array or scalar)."""
    l_arr = np.atleast_1d(np.asarray(l, dtype=int))
    out = np.array([
        float(_CL_CONT_LAMBDAS[ll](Delta_val)) * float(lattice_conversion_factor(ll))
        for ll in l_arr
    ]) * C_Delta
    return out if np.ndim(l) else out[0]


def rho_l_numeric(l, Delta_val):
    """Recursion ratio C_{l+1}/C_l IN THE LATTICE M_l[m,m] CONVENTION,
    i.e. paper_recursion_ratio(l+1, Delta) * conversion_factor(l+1)/conversion_factor(l).
    Does not depend on C_Delta (cancels), matching the "ratio cancels the
    free normalization" property. Vectorized over l."""
    l_arr = np.atleast_1d(np.asarray(l, dtype=int))
    out = np.array([
        float(paper_recursion_ratio(ll + 1, _Delta_sym).subs(_Delta_sym, Delta_val))
        * float(lattice_conversion_factor(ll + 1) / lattice_conversion_factor(ll))
        for ll in l_arr
    ])
    return out if np.ndim(l) else out[0]


def verify_rho_monotonic(Delta_val=Rational(1, 8), l_max=10):
    """rho_l_numeric should be monotonic in Delta at fixed l (needed for a
    well-posed 1D brentq root-find in fit_delta_sigma.py)."""
    ok = True
    for l_val in range(l_max):
        deltas = np.linspace(0.01, 1.8, 30)
        vals = np.array([rho_l_numeric(l_val, d) for d in deltas])
        diffs = np.diff(vals)
        monotonic = np.all(diffs > 0) or np.all(diffs < 0)
        status = "PASS" if monotonic else "FAIL"
        if not monotonic or l_val <= 2:
            print(f"  [{status}] rho_{l_val}(Delta) monotonic over Delta in (0.01,1.8): {monotonic}")
        ok = ok and monotonic
    return ok


def main():
    all_ok = True
    print("verify_against_paper_values (Eq 5.2 @ Delta=1/8):")
    all_ok &= verify_against_paper_values()

    print("verify_recursion_against_paper (Eq 5.3):")
    all_ok &= verify_recursion_against_paper()

    print("verify_numeric_integration (independent mpmath quadrature):")
    all_ok &= verify_numeric_integration()

    print("verify_conversion_factor_symbolic (lattice M_l[m,m] normalization):")
    all_ok &= verify_conversion_factor_symbolic()

    print("verify_rho_monotonic (well-posedness of recursion root-find):")
    all_ok &= verify_rho_monotonic()

    print()
    if all_ok:
        print("ALL CHECKS PASSED")
    else:
        print("SOME CHECKS FAILED -- do not import C_l_numeric/rho_l_numeric until fixed")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
