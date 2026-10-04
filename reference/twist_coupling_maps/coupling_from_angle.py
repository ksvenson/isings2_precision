"""Flat-triangle coupling-angle maps imported from the `Twist/` campaign.

Both maps below give the (dimensionless, already beta*J) Ising coupling `K`
that sits on one edge of a FLAT Euclidean triangle, as a function of that
edge's own geometry -- specifically the angle `theta` opposite it (the
other two triangle angles are irrelevant beyond fixing theta via
theta_v + theta_u + theta_w = pi). Source campaign: `../../Twist/`
(`scripts/fit_permutation_symmetric_coupling_map.py`,
`scripts/compute_geometry_boundary_angles.py`,
`campaign_runs/phase6_coupling_map/kappa_sa_fit.json`), copied in
2026-08-21. `kappa_sa_fit.json` is a static snapshot, not synced -- if
Twist's Phase 6 fit changes, re-copy deliberately and note it in this
campaign's journal.md, don't assume it stays current.

IMPORTANT / open item (see PLAN.md): these maps were derived for the flat
2D triangular lattice with three coupling directions per unit cell
(K_A, K_B, K_C). `IsingS2`'s simplicial-S2 engine already assigns each link
its own geometric weight via `QfeLatticeS2::UpdateWeights()` (a cotangent
/dual-area Regge-type weight, `link->wt`, entering the action as
`beta * link->wt * spin_i * spin_j` -- see `../../IsingS2/include/ising.h`
and `S2.h`). Whether the maps below should replace that existing weight,
validate it locally (each S2-mesh triangle is close to flat except at the
finitely many curvature-defect vertices), or are simply not applicable at
this mesh's resolution, is UNRESOLVED. Do not wire either map into a
production S2 run's coupling assignment without resolving this first.
"""
import json
import math
from pathlib import Path
from typing import Optional

_HERE = Path(__file__).resolve().parent


def exact_sinh_rule(theta_rad: float) -> float:
    """K(theta) = 0.5*asinh(cot(theta)), i.e. sinh(2K) = cot(theta).

    theta_rad: the angle (radians) of the flat triangle opposite the edge
    carrying this coupling. Exact result for the flat anisotropic
    triangular-lattice Ising model (classical star-triangle/duality
    relation) -- not an empirical fit. K0 = 0.5*asinh(1/sqrt(3)) at the
    isotropic point theta=pi/3, matching the known isotropic triangular
    Ising critical beta*J.
    """
    return 0.5 * math.asinh(1.0 / math.tan(theta_rad))


def s_a_from_angles(theta_v: float, theta_u: float, theta_w: float):
    """Shape parameters (s_i, a_i) per direction i in {v,u,w}, from angles
    in radians. delta_i = theta_i - pi/3 (sums to zero for a flat
    triangle). For direction i with the other two directions {j,k}:
        s_i = delta_j + delta_k = -delta_i   (exactly, since sum is 0)
        a_i = delta_j - delta_k
    Returns a dict {"v": (s_v,a_v), "u": (s_u,a_u), "w": (s_w,a_w)}.
    """
    d = {
        "v": theta_v - math.pi / 3.0,
        "u": theta_u - math.pi / 3.0,
        "w": theta_w - math.pi / 3.0,
    }
    out = {}
    for i, (j, k) in {"v": ("u", "w"), "u": ("v", "w"), "w": ("v", "u")}.items():
        out[i] = (d[j] + d[k], d[j] - d[k])
    return out


def _load_kappa_sa_fit(degree: int = 2, path: Optional[Path] = None) -> dict:
    path = path or (_HERE / "kappa_sa_fit.json")
    with open(path) as f:
        data = json.load(f)
    return data[str(degree)]


def empirical_kappa(s: float, a: float, degree: int = 2, path: Optional[Path] = None) -> float:
    """kappa(s,a) = sum_{p,q} coef[p,q] * s^p * a^q, monomial basis, the
    ACCEPTED production fit from Twist Phase 6
    (campaign_runs/phase6_coupling_map/kappa_sa_fit.json; degree 2 or 3
    are the real candidates, degree 10 is a discarded overfit check).
    Agrees with `exact_sinh_rule` to ~0.1-0.3% over the fitted domain.

    Note this is kappa = beta_c(K1,K2,K3) * K_i (a gauge-invariant physical
    coupling), not K_i itself in isolation -- see Twist's Phase 6 scripts
    for the full gauge-fixing context before using this as a bare K.
    """
    fit = _load_kappa_sa_fit(degree, path)
    total = 0.0
    for (p, q), c in zip(fit["terms"], fit["coef"]):
        total += c * (s ** p) * (a ** q)
    return total


if __name__ == "__main__":
    # sanity check: exact rule vs accepted degree-2 empirical fit, isotropic point
    theta0 = math.pi / 3.0
    s0, a0 = 0.0, 0.0
    print("exact K(pi/3) =", exact_sinh_rule(theta0))
    print("empirical kappa(0,0) deg2 =", empirical_kappa(s0, a0, degree=2))
