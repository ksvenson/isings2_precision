#!/usr/bin/env python3
"""Which harmonic levels CAN show a nonzero kappa2_r on a q=5 mesh.

Motivation
----------
The campaign's standing objective-1 result is "l<=2 clean, l>=3 a
non-shrinking kappa2_r plateau" (PLAN.md step 6, journal.md 2026-09-02).
This script shows the l<=2 half of that split is not evidence of anything:
it is forced by symmetry and cannot come out any other way.

Argument
--------
A q=5 mesh is exactly icosahedrally symmetric -- to machine precision for
`naive` and `equal_area` (CLAUDE.md "Mesh construction modes"; verified in
src/test_ico_symmetry.cc / src/test_gd_symmetry.cc). The lattice action,
the quadrature weights and the MC measure all inherit that symmetry, so
M_l[m,m'] = <S_lm S*_lm'> commutes with the representation of the
icosahedral group I on the (2l+1)-dimensional space spanned by Y_lm.

By Schur's lemma, M_l is therefore a scalar on each irreducible component:
if the (2l+1)-dim SO(3) irrep restricted to I decomposes into k
*inequivalent* irreps, M_l has at most k distinct eigenvalues, with
multiplicities equal to those irreps' dimensions. In particular:

  - l=1 restricts to T1 alone  -> M_1 = c*Identity EXACTLY
  - l=2 restricts to H  alone  -> M_2 = c*Identity EXACTLY

and every spectral cumulant of order >=2 vanishes identically at those two
levels, at ANY mesh resolution, however coarse. l>=3 is the first level
whose restriction contains two or more inequivalent irreps, so it is the
first level at which a nonzero kappa2_r is allowed at all.

Consequence for the gate: the two levels that "pass" PLAN.md step 6 are
exactly the two that are incapable of failing it. The l>=3 plateau is not
a l>=3-specific pathology -- it is the signal appearing at the first l
where it is kinematically permitted to appear.

This does NOT by itself explain why the l>=3 plateau fails to *shrink*
under refinement. Icosahedral symmetry is exact at every n_refine, so it
constrains which l can break SO(3), not how fast the breaking dies away.

Cross-check: this reproduces, analytically, the l=6 fact the campaign
previously found empirically via a deterministic free-scalar FEM run
(src/fem_scalar_test.cc, journal.md 2026-08-25 / PLAN.md "invariant
subspace starting at l=6, multiplicity 1") -- l=6 is the first l>0 whose
restriction contains the trivial rep A.

Committed per directive.md's "Closed-form CFT formulas" rule: the
character table and the decomposition are computed here from the group's
own class structure, not quoted from memory.

Usage:
  .venv_plot/bin/python scripts/icosahedral_harmonic_decomposition.py [--l_max 12]
"""
import argparse

import numpy as np

# Icosahedral rotation group I (order 60), by conjugacy class: the class
# sizes and the rotation angle of each class's elements.
CLASS_NAMES = ["E", "12C5", "12C5^2", "20C3", "15C2"]
CLASS_SIZES = np.array([1, 12, 12, 20, 15])
CLASS_ANGLES = np.deg2rad(np.array([0.0, 72.0, 144.0, 120.0, 180.0]))
GROUP_ORDER = int(CLASS_SIZES.sum())

PHI = (1.0 + np.sqrt(5.0)) / 2.0

# Character table of I, rows in the class order above. Dimensions are the
# E-column entries: A=1, T1=3, T2=3, G=4, H=5 (1+9+9+16+25 = 60 = |I|).
CHARACTERS = {
    "A":  np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
    "T1": np.array([3.0, PHI, 1.0 - PHI, 0.0, -1.0]),
    "T2": np.array([3.0, 1.0 - PHI, PHI, 0.0, -1.0]),
    "G":  np.array([4.0, -1.0, -1.0, 1.0, 0.0]),
    "H":  np.array([5.0, 0.0, 0.0, -1.0, 1.0]),
}


def check_character_table():
    """Verify orthonormality of the stored character table, so a typo in it
    cannot silently propagate into the decomposition below."""
    names = list(CHARACTERS)
    for i, a in enumerate(names):
        for b in names[i:]:
            inner = (CLASS_SIZES * CHARACTERS[a] * CHARACTERS[b]).sum() / GROUP_ORDER
            expected = 1.0 if a == b else 0.0
            assert abs(inner - expected) < 1e-12, f"<{a}|{b}> = {inner}, expected {expected}"
    assert sum(int(CHARACTERS[n][0]) ** 2 for n in names) == GROUP_ORDER


def so3_character(l, theta):
    """chi_l(theta) = sin((l+1/2)theta)/sin(theta/2), the SO(3) character."""
    if abs(theta) < 1e-12:
        return float(2 * l + 1)
    return float(np.sin((l + 0.5) * theta) / np.sin(theta / 2.0))


def decompose(l):
    """Multiplicity of each irrep of I in the (2l+1)-dim SO(3) irrep."""
    chi = np.array([so3_character(l, t) for t in CLASS_ANGLES])
    out = {}
    for name, char in CHARACTERS.items():
        n = (CLASS_SIZES * chi * char).sum() / GROUP_ORDER
        n_int = int(round(n))
        assert abs(n - n_int) < 1e-9, f"l={l} {name}: non-integer multiplicity {n}"
        if n_int:
            out[name] = n_int
    assert sum(CHARACTERS[k][0] * v for k, v in out.items()) == 2 * l + 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--l_max", type=int, default=12)
    args = ap.parse_args()

    check_character_table()
    print("character table of I verified orthonormal (|I| = 60)\n")
    print(f"{'l':>3} {'dim':>4}  {'restriction to I':<26} {'eigenvalue mult.':<22} "
          f"{'kappa2_r':<12} note")
    print("-" * 100)
    for l in range(args.l_max + 1):
        dec = decompose(l)
        label = " + ".join(f"{n if n > 1 else ''}{k}" for k, n in dec.items())
        mults = sorted((int(CHARACTERS[k][0]) for k, n in dec.items() for _ in range(n)),
                       reverse=True)
        # Distinct *inequivalent* irreps; a repeated irrep leaves a mixing
        # block, so it does not force extra degeneracy.
        n_distinct = len(dec)
        forced_zero = n_distinct == 1
        kappa = "== 0 exactly" if forced_zero else "may be nonzero"
        note = ""
        if forced_zero and l > 0:
            note = "single irrep -> M_l proportional to I at any resolution"
        if "A" in dec and l > 0:
            note = "contains trivial rep: icosahedrally invariant (exclude from SO(3) claims)"
        print(f"{l:>3} {2*l+1:>4}  {label:<26} {str(mults):<22} {kappa:<12} {note}")

    first = next(l for l in range(1, args.l_max + 1) if len(decompose(l)) > 1)
    print(f"\nFirst l at which a nonzero kappa2_r is symmetry-allowed: l = {first}")
    print("Levels l=1,2 cannot show SO(3) breaking in this observable at all --")
    print("their 'passing' the PLAN.md step-6 gate carries no information.")


if __name__ == "__main__":
    main()
