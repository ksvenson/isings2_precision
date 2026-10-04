#!/usr/bin/env python3
"""Legendre-coefficient conformal/spherical-symmetry test, following
reference/Owen_section_D.tex ("Agreement with the Ising CFT", Sec.
CFTtest) -- a sharper, dimension-independent companion to the campaign's
existing M_l[m,m'] eigenvalue-cumulant test (symmetry_test.py) and to
Brower et al. 2018 (arXiv:1803.08512) Sec. 5.1
(reference/brower_2018_lattice_phi4_riemann_S2.pdf), which this method
matches up to normalization convention.

Method
------
Owen_section_D defines lattice Legendre coefficients

    F_l(a) = (2l+1)/(2*pi*4*pi) * sum_ij sqrt(g_i g_j) P_l(r_i . r_j) <s_i s_j>

via the spherical-harmonic addition theorem
sum_m Y*_lm(r_i) Y_lm(r_j) = (2l+1)/(4*pi) * P_l(r_i . r_j). This campaign's
own measurement already computes the (2l+1)x(2l+1) matrix
M_l[m,m'] = <S_lm S*_lm'> with S_lm = sum_i w_i Y*_lm(x_i) s_i and
quadrature weights normalized so sum_i w_i = n_sites, NOT 4*pi as this
docstring previously (incorrectly) claimed -- confirmed 2026-08-25 via
code audit (`ising_s2_crit.cc` calls `QfeLatticeS2::UpdateWeights()`, not
`OptimizeIntegrator`, and `UpdateWeights` normalizes to `sum wt = n_sites`,
see CLAUDE.md). ising_s2_crit.cc's actual per-configuration accumulation
also divides the S_lm*conj(S_lmp) product by n_sites^2 before writing it
to the jackblocks file the C++ side calls this quantity `M_l[m,m']` but,
combined with the sum-to-n_sites weights, it differs from the algebraic
M_l[m,m'] used in this derivation by a uniform (l-independent) factor of
1/(4*pi)^2 -- verified 2026-08-25 that this factor is exactly constant
across l, so it cancels completely in every ratio-based diagnostic this
module computes (delta_l, Delta_s, Delta_l_pair) and does not affect any
conclusion drawn from them; it does mean this module's F_0 is NOT
literally equal to Owen_section_D's `2*<M(a)^2>` in absolute units (off
by (4*pi)^2 ~ 157.9), so do not use F_0 for that specific absolute
cross-check without correcting for it. Also note `S_lm` as actually
computed in C++ uses `Y_lm` (not the conjugate `Y*_lm` this docstring
writes) -- shown 2026-08-25 to leave every diagonal (trace, i.e. F_l) and
magnitude-based (R_l/kappa_r cumulants) quantity unchanged, since tracing
or taking |.| of a globally-conjugated Hermitian matrix is invariant; it
only flips the phase of individual off-diagonal M_l[m,m'] (m!=m') entries
relative to this docstring's stated convention, which no downstream
script currently reads directly. sqrt(g_i g_j) above plays the same role
as w_i w_j regardless of these normalization details.
Tracing M_l over m and applying the same addition theorem gives

    Trace(M_l) = sum_m M_l[m,m] = (2l+1)/(4*pi) * sum_ij w_i w_j P_l(x_i.x_j) <s_i s_j>

so, comparing coefficients directly,

    F_l(a) = Trace(M_l) / (2*pi)                                        (*)

This is computed directly here from the *_ylm_2pt_full_jackblocks_*.dat
files already produced for the existing kappa2_r/R_l symmetry test --
no new C++ measurement or production run is needed. (This also sidesteps
CLAUDE.md's flagged-open-item *_legendre_2pt_*.dat proxy: (*) is an
independently-derived formula for the same quantity, verified against
Owen_section_D's stated normalization, not a re-use of the unconfirmed
on-disk proxy.)

Diagnostics (Owen_section_D, eliminating Delta from the F_l^cont recursion
F_l^cont = [(l-1+Delta)/(l+1-Delta)] * F_{l-1}^cont):

  delta_l(a) = 1 - (F_0-F_1)(F_{l-1}+F_l) / [l*(F_0+F_1)(F_{l-1}-F_l)]
      dimension-independent conformal/spherical-symmetry-breaking measure,
      built only from measured F_l -- should -> 0 in the continuum limit.
      Defined for l >= 2 (needs F_{l-1}, F_l with l-1 >= 1).

  Delta_s(a) = 2*F_1 / (F_1 + F_0)
      direct lattice estimator of the spin scaling dimension (exact value
      1/8), from inverting the l=1 recursion step alone.

All point estimates come from the full (pooled) dataset; errors are
leave-one-block-out jackknife over the same blocks symmetry_test.py uses,
propagated through the (nonlinear) delta_l/Delta_s formulas per jackknife
sample -- consistent with this campaign's existing error-propagation
convention (see symmetry_test.py's docstring).
"""

import argparse
import sys

import numpy as np

from symmetry_test import read_jackblocks_multi


def diag_sums_per_l(l_max, blocks):
    """total_diag_sum[l] = sum over all blocks/pairs of sum[(l,m,m)] (complex,
    but imaginary part must vanish since M_l is Hermitian); per_block_diag[l]
    is the same per individual block, for jackknife leave-one-out.
    """
    n_blocks = len(blocks)
    total_diag = np.zeros(l_max + 1, dtype=complex)
    block_diag = np.zeros((n_blocks, l_max + 1), dtype=complex)
    for k, blk in enumerate(blocks):
        for (l, m, mp), val in blk["sum"].items():
            if m != mp:
                continue
            block_diag[k, l] += val
            total_diag[l] += val
    return total_diag, block_diag


def compute_F_l(l_max, blocks):
    """Returns (F_full[l_max+1], F_jk[n_blocks, l_max+1]) via Eq (*) above,
    DIVIDED by (2l+1) -- bug found and fixed 2026-08-26 (see journal.md's
    "F_l(a) missing (2l+1) division" entry). Trace(M_l)/(2*pi) (the
    quantity this function computed before the fix) is NOT F_l^cont in the
    continuum limit; it is (2l+1)*F_l^cont.

    Derivation: F_l(a) := Trace(M_l)/(2*pi) = (2l+1)/(8*pi^2) *
    sum_ij w_i w_j P_l(x_i.x_j) <s_i s_j> (via the addition theorem). Take
    the continuum limit with w_i normalized so sum_i w_i g(x_i) ->
    int dOmega g(x) (i.e. true solid-angle quadrature weights, matching
    Owen_section_D's F_0(0)=2<M^2> normalization check -- an
    l-independent amplitude mismatch from this codebase's actual
    sum_i w_i = n_sites convention cancels in every ratio here regardless,
    same as documented in this module's docstring). Expand the two-point
    function G(x_i.x_j) = sum_l 2*pi*F_l^cont(Delta) * sum_m
    Y_lm(x_i) Y*_lm(x_j) (inverting Owen_section_D Eq 25's G(z) =
    sum_l (2l+1)/2 * F_l^cont * P_l(z) via the addition theorem). Then

      int dOmega_i dOmega_j [sum_m Y*_lm(x_i)Y_lm(x_j)] G(x_i,x_j)
        = sum_l' 2*pi*F_l'^cont * sum_{m,m'}
              [int dOmega_i Y*_lm(x_i)Y_l'm'(x_i)]
              [int dOmega_j Y_lm(x_j)Y*_l'm'(x_j)]
        = 2*pi*F_l^cont * (2l+1)                       (orthonormality
                                                          collapses l'=l,
                                                          m'=m, leaving a
                                                          bare sum over
                                                          the (2l+1)
                                                          values of m)

    so F_l(a) -> (1/(2*pi)) * 2*pi*(2l+1)*F_l^cont = (2l+1)*F_l^cont, NOT
    F_l^cont. delta_s_formula/delta_l_pair_formula below were derived
    against F_l^cont's recursion (Owen_section_D Eq 60-63) and require
    their F_0/F_1/... inputs to actually BE F_l^cont -- feeding them
    Trace(M_l)/(2*pi) directly multiplies the l=1/l=0 ratio by an
    unaccounted (2*1+1)/(2*0+1)=3, which is why Delta_s(a) computed this
    way plateaus near ~1/3 (0.33-0.4 with lattice corrections) instead of
    1/8: 3*0.125=0.375, matching the observed plateau almost exactly.
    Confirmed numerically: dividing by (2l+1) moved a lean_harmonic_stats
    n_refine=4 test run's Delta_s from 0.3595 to 0.1362 -- landing right
    where the real-space method's n_refine=4 point (0.139) sits, the
    correct sign/magnitude for a finite-lattice-spacing correction, unlike
    the uncorrected value.

    This bug predates this fix and affects every prior Trace(M_l)-based
    Delta_s/Delta_l_pair result from this module (cft_symmetry_test.py,
    build_ylm_matrix_from_configs.py, build_ylm_matrix_from_corr.py, the
    original in-simulator ising_s2_crit.cc l-space measurement) -- none of
    those numbers should be trusted without rerunning through this fix.
    The single-row north-pole method (fl_delta_from_north_pole.py) is
    UNAFFECTED: F_l there is already a direct single Legendre integral
    (literally F_l^cont's lattice discretization, no addition-theorem/
    double-sum step), which is exactly why removing ITS erroneous (2l+1)
    fixed it while this construction needs the opposite correction.
    """
    n_blocks = len(blocks)
    total_n = sum(blk["n"] for blk in blocks)
    block_n = np.array([blk["n"] for blk in blocks])

    total_diag, block_diag = diag_sums_per_l(l_max, blocks)

    max_imag = np.max(np.abs(total_diag.imag))
    max_real = np.max(np.abs(total_diag.real))
    if max_real > 0 and max_imag / max_real > 1e-8:
        print(
            f"warning: Trace(M_l) has non-negligible imaginary part "
            f"(max|Im|/max|Re|={max_imag / max_real:.2e}) -- Hermiticity "
            f"assumption may be violated",
            file=sys.stderr,
        )

    two_l_plus_1 = 2.0 * np.arange(l_max + 1) + 1.0

    F_full = total_diag.real / total_n / (2.0 * np.pi) / two_l_plus_1

    F_jk = np.empty((n_blocks, l_max + 1))
    for k in range(n_blocks):
        loo_n = total_n - block_n[k]
        loo_diag = (total_diag - block_diag[k]).real
        F_jk[k] = loo_diag / loo_n / (2.0 * np.pi) / two_l_plus_1

    return F_full, F_jk


def jackknife_err(full, jk):
    n = len(jk)
    bar = np.mean(jk, axis=0)
    return np.sqrt((n - 1) / n * np.sum((jk - bar) ** 2, axis=0))


def delta_l_formula(F0, F1, Flm1, Fl, l):
    denom = l * (F0 + F1) * (Flm1 - Fl)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 - (F0 - F1) * (Flm1 + Fl) / denom


def delta_s_formula(F0, F1):
    return 2.0 * F1 / (F1 + F0)


def delta_l_pair_formula(Flm1, Fl, l):
    """Delta estimate from inverting the single recursion step at this l
    alone (l>=1): rho=F_l/F_{l-1}, Delta_l = [rho*(l+1)-(l-1)]/(1+rho).
    Generalizes delta_s_formula, which is just this at l=1 (reduces to
    2*F1/(F1+F0) exactly, since l-1=0, l+1=2 there). If the ensemble is a
    genuine single-Delta CFT, every l should give the same Delta_l --
    checking that agreement directly (not just trusting l=1) is the point;
    see journal.md 2026-08-25 user prompt."""
    rho = Fl / Flm1
    return (rho * (l + 1) - (l - 1)) / (1.0 + rho)


def analyze(paths, min_blocks_factor=5):
    l_max, jack_block_size, blocks = read_jackblocks_multi(paths)
    n_blocks = len(blocks)
    F_full, F_jk = compute_F_l(l_max, blocks)

    # Delta_s(a)
    ds_full = delta_s_formula(F_full[0], F_full[1])
    ds_jk = delta_s_formula(F_jk[:, 0], F_jk[:, 1])
    ds_err = jackknife_err(ds_full, ds_jk)

    results = dict(
        l_max=l_max, n_blocks=n_blocks,
        F_full=F_full, F_err=jackknife_err(F_full, F_jk),
        delta_s=ds_full, delta_s_err=ds_err,
        delta_l={}, delta_l_err={},
        delta_l_pair={}, delta_l_pair_err={},
    )

    for l in range(2, l_max + 1):
        dl_full = delta_l_formula(F_full[0], F_full[1], F_full[l - 1], F_full[l], l)
        dl_jk = delta_l_formula(F_jk[:, 0], F_jk[:, 1], F_jk[:, l - 1], F_jk[:, l], l)
        results["delta_l"][l] = dl_full
        results["delta_l_err"][l] = jackknife_err(dl_full, dl_jk)

    # Delta_l_pair(a): the single-step recursion inversion at every l,
    # not just l=1 (delta_s_formula is the l=1 case of this). If this is
    # really one CFT primary of dimension Delta, every l should agree.
    for l in range(1, l_max + 1):
        dp_full = delta_l_pair_formula(F_full[l - 1], F_full[l], l)
        dp_jk = delta_l_pair_formula(F_jk[:, l - 1], F_jk[:, l], l)
        results["delta_l_pair"][l] = dp_full
        results["delta_l_pair_err"][l] = jackknife_err(dp_full, dp_jk)

    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "jackblocks_path", nargs="+",
        help="*_ylm_2pt_full_jackblocks_*.dat file(s); multiple shard files "
        "for the same (n_refine, coupling_rule, mesh_mode) point are pooled",
    )
    args = ap.parse_args()

    r = analyze(args.jackblocks_path)

    print(f"# {args.jackblocks_path}")
    print(f"# l_max={r['l_max']} n_blocks={r['n_blocks']}")
    print(f"# Delta_s(a) = 2*F_1/(F_1+F_0) [exact value 1/8 = 0.125]:")
    print(f"#   {r['delta_s']:.8f} +/- {r['delta_s_err']:.8f}")
    print(f"# F_l(a):")
    for l in range(r["l_max"] + 1):
        print(f"#   l={l:2d}  F_l={r['F_full'][l]:14.6e} +/- {r['F_err'][l]:14.6e}")
    print(f"# delta_l(a) [dimension-independent symmetry-breaking measure, -> 0 in continuum limit]:")
    print(f"{'l':>3} {'delta_l':>16} {'err':>16}")
    for l in sorted(r["delta_l"]):
        print(f"{l:3d} {r['delta_l'][l]:16.6e} {r['delta_l_err'][l]:16.6e}")
    print(f"# Delta_l_pair(a) [single-step recursion inversion at each l -- should all "
          f"agree with each other AND with Delta_s=Delta_1_pair if this is one CFT "
          f"primary; l=6 (and other icosahedral-invariant l for this mesh's q) is "
          f"expected to be contaminated regardless of continuum limit -- see journal.md "
          f"2026-08-25 fem_scalar_test finding, flag/exclude rather than trust]:")
    print(f"{'l':>3} {'Delta_l_pair':>16} {'err':>16}")
    for l in sorted(r["delta_l_pair"]):
        print(f"{l:3d} {r['delta_l_pair'][l]:16.6e} {r['delta_l_pair_err'][l]:16.6e}")


if __name__ == "__main__":
    main()
