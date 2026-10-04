#!/usr/bin/env python3
"""Convention-free numerical confirmation that M_l is forced scalar for l<=2.

Companion to scripts/icosahedral_harmonic_decomposition.py, which derives
the same result from the icosahedral character table. This script derives
it WITHOUT any character table, without a Wigner D matrix, and without
committing to any spherical-harmonic phase convention -- and it uses the
campaign's OWN group data (grp/elem/o3q5.dat, the same file the simulator
reads), so it also closes the loop that the abstract group in that other
script is the concrete group the code actually uses.

Method: for each icosahedral rotation R, DISCOVER the matrix T(R) that
implements R in the Y_lm basis by least squares -- i.e. solve
Y_lm(R^-1 x) = sum_m' T(R)_{m'm} Y_lm'(x) for T, rather than asserting
what T is. Then check, at each l:

  (a) the least-squares residual is ~0
      -> span{Y_lm : m} really is closed under rotation. This is the one
         non-trivial input to the Schur argument, and here it is measured
         rather than assumed. (If the span were not invariant, no T could
         fit and the residual would be O(1).)
  (b) T(R) is unitary
  (c) T(R1) T(R2) = T(R1 R2) -> the T's are a genuine representation
  (d) dim{M : [T(g), M] = 0 for all g}, the commutant dimension.
      M_l is forced proportional to the identity iff this equals 1.

Any stray phase convention in Y_lm would change T by conjugation with a
diagonal unitary, which changes neither unitarity, nor the representation
property, nor the commutant dimension -- so (d) is convention-independent
by construction, which is the point of doing it this way.

Expected output: commutant dim 1 for l=0,1,2 (kappa2_r identically zero)
and >1 from l=3 up, matching icosahedral_harmonic_decomposition.py's
sum_a n_a^2 exactly.

Usage:
  .venv_plot/bin/python scripts/verify_harmonic_rep_numerically.py
(run from the repo root -- it reads grp/elem/o3q5.dat by relative path)
"""
import numpy as np
from scipy.special import sph_harm_y

# --- icosahedral group from the repo's own data (proper rotations only) ---
rows = np.loadtxt("grp/elem/o3q5.dat")
proper = rows[rows[:, 1] > 0]
assert len(proper) == 60, len(proper)

def quat_to_R(w, x, y, z):
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
        [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)]])

G = [quat_to_R(*r[2:6]) for r in proper]
for R in G:
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-12)
    assert abs(np.linalg.det(R) - 1) < 1e-12
print(f"loaded {len(G)} proper rotations from grp/elem/o3q5.dat; all orthogonal, det=+1")

# closure check: it really is a group
Gstack = np.stack(G)
def find(R):
    d = np.abs(Gstack - R).reshape(len(G), -1).max(axis=1)
    i = int(d.argmin())
    assert d[i] < 1e-9, f"product not found in group (closest dist {d[i]:.2e})"
    return i
mult = np.array([[find(G[a] @ G[b]) for b in range(60)] for a in range(60)])
print("closed under multiplication (verified all 3600 products)")

# --- sample points, build Y and rotated-Y, discover T(R) by least squares ---
rng = np.random.default_rng(7)
npts = 400
v = rng.normal(size=(npts, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)

def Ymat(pts, l):
    th = np.arccos(np.clip(pts[:, 2], -1, 1))
    ph = np.arctan2(pts[:, 1], pts[:, 0])
    return np.stack([sph_harm_y(l, m, th, ph) for m in range(-l, l+1)], axis=1)

print(f"\n{'l':>3} {'max fit resid':>14} {'max |T T^H - I|':>16} {'max rep error':>14} "
      f"{'commutant dim':>14} {'forced scalar?':>15}")
print("-"*82)
for l in range(0, 9):
    Y = Ymat(v, l)
    Ts, resid = [], 0.0
    for R in G:
        Yrot = Ymat(v @ R, l)            # x -> R^{-1}x  is  rows v @ R
        T, *_ = np.linalg.lstsq(Y, Yrot, rcond=None)
        resid = max(resid, np.abs(Y @ T - Yrot).max())
        Ts.append(T)
    unit = max(np.abs(T @ T.conj().T - np.eye(2*l+1)).max() for T in Ts)
    reperr = max(np.abs(Ts[a] @ Ts[b] - Ts[mult[a, b]]).max()
                 for a in range(0, 60, 7) for b in range(0, 60, 11))

    # commutant: stack [T,M]=0 as a linear system on vec(M)
    d = 2*l+1; I = np.eye(d)
    A = np.vstack([np.kron(I, T) - np.kron(T.T, I) for T in Ts])
    sv = np.linalg.svd(A, compute_uv=False)
    cdim = int((sv < 1e-8 * max(sv[0], 1.0)).sum()) + (d*d - len(sv))
    print(f"{l:>3} {resid:>14.2e} {unit:>16.2e} {reperr:>14.2e} {cdim:>14d} "
          f"{'YES' if cdim == 1 else 'no':>15}")
