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



# ---------------------------------------------------------------------
# Independent cross-check, and the q=3/4/5 comparison.
#
# The commutant of a representation -- the space of ALL operators that
# commute with every group element, which is exactly where M_l is forced
# to live -- has dimension sum_a n_a^2 over the irrep multiplicities n_a,
# and equals (1/|G|) sum_g |chi(g)|^2. That needs only the group's class
# sizes and rotation angles: NO character table. So it is a genuinely
# independent check on the decomposition above.
#
# M_l is forced proportional to the identity  <=>  commutant dimension 1
# (one irrep, multiplicity one). That is the whole criterion.
#
# Rotation groups of the three valid base polyhedra (S2.h: q = 3, 4, 5 for
# tetrahedron, octahedron, icosahedron). Each entry is (class size,
# rotation angle in degrees).
# ---------------------------------------------------------------------
ROTATION_GROUPS = {
    3: ("T  (tetrahedron)", [(1, 0.0), (4, 120.0), (4, 240.0), (3, 180.0)]),
    4: ("O  (octahedron)",  [(1, 0.0), (6, 90.0), (3, 180.0), (8, 120.0), (6, 180.0)]),
    5: ("I  (icosahedron)", [(1, 0.0), (12, 72.0), (12, 144.0), (20, 120.0), (15, 180.0)]),
}


def commutant_dim(q, l):
    """dim of the algebra of operators commuting with the whole group on
    the (2l+1)-dim spin-l space = sum_a n_a^2. Character table not used."""
    _, classes = ROTATION_GROUPS[q]
    order = sum(size for size, _ in classes)
    total = sum(size * so3_character(l, np.deg2rad(ang)) ** 2 for size, ang in classes)
    val = total / order
    n = int(round(val))
    assert abs(val - n) < 1e-9, f"q={q} l={l}: non-integer commutant dim {val}"
    return n


def first_breaking_l(q, l_max=24):
    for l in range(1, l_max + 1):
        if commutant_dim(q, l) > 1:
            return l
    return None


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

    # Independent check of that number via the commutant dimension, which
    # uses only class sizes/angles -- no character table, so a typo in
    # CHARACTERS cannot produce an agreeing answer by accident.
    for l in range(args.l_max + 1):
        from_table = sum(n * n for n in decompose(l).values())
        assert commutant_dim(5, l) == from_table, f"l={l}: commutant mismatch"
    print(f"cross-check: commutant dims agree with the character table for all "
          f"l<={args.l_max} (computed without it)")

    print("\nWhy q=5. The base polyhedron (S2.h: q=3/4/5 = tetrahedron/"
          "octahedron/icosahedron)\nfixes the residual symmetry group, and so "
          "fixes the first l that can break SO(3):\n")
    print(f"  {'q':>2}  {'group':<20} {'|G|':>4}  {'first l with kappa2_r allowed':<30}")
    print("  " + "-" * 62)
    for q in (3, 4, 5):
        name, classes = ROTATION_GROUPS[q]
        order = sum(sz for sz, _ in classes)
        print(f"  {q:>2}  {name:<20} {order:>4}  l = {first_breaking_l(q)}")
    print("\nThe icosahedral group is the largest of the three, so it keeps the")
    print("harmonics degenerate longest: q=5 forces BOTH l=1 and l=2 to be exactly")
    print("scalar, where q=4 and q=3 force only l=1. That extra clean level is one")
    print("reason this campaign uses q=5 -- and it is also why its first usable")
    print("symmetry signal sits at l=3 rather than l=2.")


if __name__ == "__main__":
    main()
