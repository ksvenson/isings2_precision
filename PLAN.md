# High-precision spherical-symmetry and CFT-data verification for 2D Ising on S^2

Version: v001
Date: 2026-08-21
Scope: scientific plan only; no jobs are launched by this document.

## Scientific objective

Simulate the critical 2D Ising model directly on a simplicial discretization
of the round sphere S^2, and use the spherical harmonic decomposition of the
spin two-point function to:

1. **Verify spherical symmetry** — check that the discretized theory at
   criticality respects the full SO(3) isometry of the sphere to controlled
   precision, as a diagnostic of simplicial discretization artifacts (mesh
   curvature-defect concentration, quadrature-weight optimization order,
   finite mesh resolution).
2. **Extract CFT data** — the scaling dimensions Delta_sigma (spin operator)
   and Delta_epsilon (energy operator), by fitting the harmonic-level
   eigenvalues of the two-point function to the closed-form prediction for a
   CFT primary on the sphere, and comparing to the exact 2D Ising values
   (Delta_sigma = 1/8, Delta_epsilon = 1).

This method (project a lattice operator onto spherical harmonics, build the
resulting `(2l+1)x(2l+1)` correlation matrix per harmonic level `l`, read off
Delta from the eigenvalues) parallels the "fuzzy sphere" harmonic-decomposition
approach used elsewhere in the CFT-on-sphere numerics literature. Find and
cite the specific papers/conventions before Phase 1 closes (open item below)
rather than reimplementing from memory.

## Starting point

No preloaded data. This is a from-scratch campaign; nothing is carried over
from `../IsingS2/` except the reusable engine code (mesh generation, Ising
updates, Ylm/quadrature-weight infrastructure — see `CLAUDE.md`). No old run
output, checkpoint, or fit from `../IsingS2/data/` may be treated as
production input.

## Method

### Operators

- **Spin**: `S_lm = sum_i w_i * Y_lm(x_i) * s_i`, where `s_i in {+1,-1}` is
  the site spin, `x_i` its position on the unit sphere, and `w_i` the
  quadrature weight from `QfeLatticeS2::OptimizeIntegrator(l_max)`
  (`sites[i].wt`), normalized so `sum_i w_i = 4*pi`.
- **Energy**: analogous projection of a bond/plaquette energy density
  observable onto `Y_lm`; exact lattice definition (bond energy assigned to
  which site/face, and how it is symmetrized) is an open item for Phase 1.

### Correlation matrices and the symmetry/CFT tests

For each harmonic level `l`, form the `(2l+1)x(2l+1)` matrix
`M_l[m,m'] = <S_lm S*_lm'>` averaged over independent thermalized
configurations (and, for the spin operator at the critical point, over the
Z2 sector — `<S_lm>` itself vanishes by symmetry and only the sector-averaged
second moment is physical, matching how magnetization-squared observables
are normally handled at criticality).

- **Symmetry test**: on the exact sphere, SO(3) invariance forces `M_l` to be
  proportional to the identity (`M_l[m,m'] = C_l * delta_{m,m'}`), with no
  dependence on the mesh's discrete orientation. Off-diagonal power and
  variation of the diagonal across `m` measure the discretization's SO(3)
  breaking, and should shrink under mesh refinement (increasing `l_max` /
  finer triangulation) if the continuum limit is being approached correctly.
- **CFT-data test**: for a CFT two-point function `<O(n)O(n')> = C_Delta /
  (2*(1-cos(gamma)))^Delta` (gamma = geodesic angle between n, n'), the
  eigenvalue `C_l` has a known closed form as a function of `l` and `Delta`
  from expanding `(1-cos(gamma))^{-Delta}` in Legendre polynomials / the
  spherical harmonic addition theorem. **Derive and symbolically verify this
  closed form in Phase 1** (e.g. with sympy, cross-checked against known
  special cases: Delta -> integer limits, and the free-fermion value
  Delta_epsilon = 1) before using it in any fit — do not trust a
  from-memory formula. Fit `Delta` from the measured `C_l(l)` shape
  (l-dependence only; overall normalization `C_Delta` is a free nuisance
  parameter, not part of the fit target).

### Mesh

Reuse `../IsingS2/`'s simplicial S2 mesh construction
(icosahedral-refinement family, `QfeLatticeS2`). Resolution is set by the
refine level `n_refine` (or equivalently `n_sites`); `l_max` for the harmonic
decomposition must stay well below what the mesh resolution can resolve
(rule of thumb, to be validated empirically in Phase 1: `l_max` a small
fraction of `sqrt(n_sites)`, since that sets the shortest wavelength the
triangulation supports without aliasing). Multiple mesh resolutions are
required to see discretization artifacts shrink and to take a continuum
extrapolation.

### Critical coupling (resolved for the spherical-symmetry test, 2026-08-22)

**Superseded.** This section originally assumed a single tunable global `K`
that would need locating via a Binder-cumulant crossing scan. As of
2026-08-22, `ising_s2_crit.cc`'s `--coupling_rule exact_sinh` assigns every
link's `K` directly from local mesh geometry
(`K = 0.5*asinh(cot(theta))`, averaged over the link's two adjacent faces)
with **no free global scale factor** — there is nothing left to scan.

Per user direction, the spherical-symmetry-test phase (below) takes this
geometry-fixed coupling as given: criticality is asserted by construction
(the same local self-dual/star-triangle logic behind `Twist/`'s flat-lattice
result), not verified by a K_c crossing search. Whether the locally
self-dual assignment produces a genuinely critical (scale-invariant) system
once assembled on a curved, defected mesh is exactly what the
spherical-symmetry test itself checks — indirectly, via whether `R_l`/the
diagonal chi-square shrink and `U4` stabilizes under mesh refinement, not
via a separate tuning step. `--coupling_rule duality` (the pre-existing
dual-link-length formula) remains available as a free secondary comparison
run at the same resolutions — informative, but non-blocking for the
symmetry-test conclusion.

### Local coupling assignment (imported from `Twist/`)

`../Twist/` (a sibling campaign mapping flat-triangle geometry to critical
Ising couplings) supplies two candidate maps from a flat triangle's
per-edge opposite angle `theta` to that edge's coupling `K`, now vendored
into `reference/twist_coupling_maps/coupling_from_angle.py` (snapshot
copied 2026-08-21, not live-synced to `Twist/`):

- **Exact rule**: `K(theta) = 0.5*asinh(cot(theta))`
  (`sinh(2K)=cot(theta)`), the classical flat-triangular-lattice
  star-triangle/duality result — not fit to any data.
- **Empirical rule**: `kappa(s,a)`, Twist's accepted Phase 6 monomial fit
  (`kappa_sa_fit.json`, degree 2/3), agreeing with the exact rule to
  ~0.1-0.3% over Twist's fitted domain.

**Decision (2026-08-21, user direction)**: run the first precision
production pass using the **exact sinh rule** to assign local per-edge
couplings on the S2 mesh (each mesh triangle's own opposite angle feeding
`K(theta)`, replacing or benchmarked against `IsingS2`'s existing
`link->wt` cotangent/dual-area weight — resolving that comparison is part
of this pass, see open items). Only if that checks out (spherical symmetry
+ correct CFT scaling dimensions to target precision) do we separately
verify the empirical `kappa(s,a)` rule as a second, independent pass. The
empirical rule is explicitly second in line, not run in parallel with the
exact-rule check.

## Open items (resolve before Phase 1 closes)

- Pin down the specific literature source/convention for the eigenvalue
  closed form and confirm this campaign's operator normalization matches it.
- Decide the exact lattice energy-density operator definition.
- Decide the statistics policy: number of independent configurations per
  mesh resolution/l_max needed for the eigenvalue precision this campaign
  targets, and how errors are estimated (jackknife/bootstrap over
  configurations, consistent with autocorrelation time).
- SCC build of `../IsingS2/`'s C++ engine (see `CLAUDE.md` — Makefile is
  currently Mac-only).
- Confirm whether `S_lm` measurement needs a new small C++ measurement pass
  bolted onto `isingS2.cc`, or whether raw spin configurations should be
  checkpointed and all harmonic-transform analysis done offline in Python
  (Eigen already computes/caches `Y_lm` on the C++ side via `UpdateYlm`, so
  there is a real choice here about where projection happens).
- **(2026-08-25, central open question)** `Delta_l_pair(a)` (the recursion
  inversion generalized to every adjacent `(l-1,l)` pair, not just l=0,1 --
  `scripts/cft_symmetry_test.py`) shows l=1 gives ~0.33 while every l=2..8
  independently agrees on a completely different value, ~0.58-0.63 --
  stable across mesh resolution (32 vs 128), mesh mode (naive/equal_area/
  equal_rp), and coupling rule (exact_sinh/duality). Neither cluster is
  near the exact value 1/8. Rules out mesh non-uniformity, resolution, and
  the l=6 icosahedral-invariant contamination (see next item) as the
  explanation. See `journal.md` 2026-08-25 for detail and open hypotheses
  — this supersedes the earlier "Delta_s(a) plateaus at ~0.33" framing as
  the campaign's sharpest open anomaly. **RESOLVED (2026-08-25, same
  session): `equal_rp` does NOT fix this either, at production
  statistics.** `push8_2026-08-25_eqrp_exactsinh` (ladder
  `n_refine={4,6,8,12,16,24,32}`, `l_max=8`, `--coupling_rule exact_sinh`,
  32 shards/point, 224/224 tasks completed) gives the identical split at
  its best-resolved points: l=1 -> 0.337-0.338, l=2..8 -> 0.58-0.63,
  bit-for-bit the same shape as `naive`/`equal_area`. Mesh non-uniformity
  (resolution, area, and now joint circumradius+perimeter+area) is fully
  ruled out as the explanation — see `journal.md`'s "equal_rp does NOT fix"
  entry for the full table. **This is now the central open question for
  the campaign** — the three untested hypotheses ((a) admixture of a
  second, non-CFT-primary operator, (b) `F_0`'s zero-mode having special
  normalization, (c) a coupling/criticality distortion specific to
  long-wavelength modes) are the next thing to investigate, not another
  mesh-smoother variant.
- **(2026-08-25, same session) Light Binder-cumulant (U4) scan
  (`scripts/binder_scan.py`) argues the system IS genuinely scale
  invariant, disfavoring hypothesis (c) above.** On
  `push7_2026-08-25_dual_naive_eqarea`'s naive-mesh ladder
  (`n_refine=2..128`, `exact_sinh`), U4 rises steeply out to
  `n_refine~16` then plateaus (0.8486 at 32 -> 0.8522 at 128, a small,
  slowly-decaying drift, not a runaway) at a value in the right ballpark
  for the known 2D Ising Binder universal amplitude (~0.856-0.861). A
  bulk order parameter that looks this well-behaved, paired with a sharp,
  resolution/mesh-mode/coupling-rule-independent anomaly only in the
  harmonic-decomposition observable (`Delta_l_pair(a)`), points at a bug
  or convention error in the `S_lm`/`M_l[m,m']`/`F_l` analysis pipeline
  rather than a real failure of criticality. **Next step: deep audit of
  the harmonic-decomposition code path** (`ising_s2_crit.cc`'s `S_lm`
  measurement, `scripts/cft_symmetry_test.py`'s `F_l`/`Delta_l_pair`
  derivation) — see `journal.md` for the audit findings.
  **Audit done (2026-08-25, same session): no code bug found.** The
  `F_l = Trace(M_l)/(2*pi)` formula and the `Delta_l_pair` recursion
  inversion were both re-derived by hand against
  `reference/Owen_section_D.tex` and match exactly; `S_lm` array indexing
  is correct. Two real discrepancies were found and fixed in
  documentation (`ising_s2_crit.cc` calls `UpdateWeights()`, not the
  `OptimizeIntegrator(l_max)` harmonic-exact quadrature solver
  `CLAUDE.md` claimed; and the code's `S_lm` uses `Y_lm` not `Y*_lm`) but
  both are provably **l-independent** (verified algebraically) and so
  cancel exactly in every ratio-based diagnostic (`Delta_l_pair`,
  `Delta_s`, `delta_l`) — confirmed not the cause. The coincident-point
  ("self-contact", `i=j`) lattice artifact is real and l-dependent but
  three-plus orders of magnitude too small (`~2e-5` vs. measured `F_l`
  `~1e-4`-`~4e-3`) to explain the split. **Combined with the Binder scan
  above, this now favors hypothesis (a)** (a genuine admixture of a
  second, lattice-artifact operator into the interacting Ising two-point
  function — something the deterministic free-scalar `fem_scalar_test.cc`
  structurally cannot probe even though it shares the exact same
  mesh/weight/harmonic-projection code and comes back clean) **over (b)
  or (c)**. Next step: fit `F_l(a)`/`<s_i s_j>` to a two-term model (CFT
  primary + a second operator) instead of a single `Delta`, not another
  code audit or mesh-uniformity pass. See `journal.md`'s full audit entry
  for the item-by-item derivation.
- **(2026-08-25, same session) The per-l `Delta_l_pair(a)` values are
  individually STABLE (precisely, cleanly converging) — the anomaly is
  that there are (at least) two different stable numbers, not that any
  one of them fails to converge.** Continuum extrapolation
  (`a_inf + b/n_refine`, weighted fit, push7's full 13-point naive
  ladder) gives l=1 -> 0.3322(6) with a real, 10+ sigma O(1/n_refine)=O(a)
  correction (still visibly draining away even at n_refine=128); l=2 ->
  0.5790(13) and l=3 -> 0.6055(15), both already flat (slope consistent
  with zero) by n_refine=8; l=8 -> ~0.617(5) restricted to its
  best-resolved points. **l=1's O(a) correction is the wrong power for
  an ordinary irrelevant lattice artifact** (those go as a^2 or faster)
  — the right signature for a relevant/marginal operator contaminating
  specifically the longest-wavelength channel. Even the "l>=2 cluster"
  shows a small but real drift with l (0.579 -> 0.606 -> ~0.617), not one
  exactly-shared constant. See `journal.md`'s "individually STABLE"
  entry for the full fit table. Next step unchanged: a genuine
  two-operator joint fit of `F_l(a)` across all l, using these
  individually-precise per-l numbers as the fit target.
- **(2026-08-25, same session) Done: two-operator fit
  (`scripts/two_operator_fit.py`), `F_l = A*F_l^cont(Delta_1) +
  B*F_l^cont(Delta_2)`.** Decisively preferred over a single power law
  (chi2/dof improves from ~7500 to <1 at n_refine>=16). At n_refine=128:
  `Delta_1=0.0724(35)`, `Delta_2=0.6364(14)`, chi2/dof=0.42. `Delta_2` is
  stable ~0.63-0.65 from n_refine=16 on. `Delta_1` is much closer to 1/8
  than any single-l estimate (vs. l=1's clean-but-wrong 0.332
  extrapolation) but a naive `a_inf+b/n_refine` extrapolation lands at
  `0.060(2)`, still notably below 0.125 and with an imperfect
  extrapolation fit (chi2/dof=2.95) — **not yet a clean confirmation of
  1/8**. Two open caveats before trusting `Delta_1` precisely: (1) fit
  used a diagonal chi-square, not the full jackknife covariance across l
  (every `F_l` at a point comes from the same jackknife blocks, so they're
  correlated); (2) the continuum-extrapolation functional form
  (`a_inf+b/n_refine`) hasn't been validated for this fitted, nonlinear
  parameter. Next steps: build the full covariance into the fit, try
  `+c/n_refine^2`, and/or fit `A,B,Delta_1,Delta_2` jointly across the
  whole ladder simultaneously rather than point-by-point then
  extrapolate. See `journal.md`'s "two-operator fit" entry for the full
  table and reasoning.
- **(2026-08-25, same session) REJECTED: the second operator is NOT
  epsilon (Delta=1).** Fixing `Delta_1=1/8` (exact sigma) and letting
  `Delta_2` float still converges to `~0.656-0.660`, not 1, and
  degrades chi2/dof from <1 to 5-23 just from fixing `Delta_1`. Fixing
  BOTH to their exact CFT values (`Delta_1=1/8`, `Delta_2~1`, epsilon)
  gives chi2/dof ~16,000-20,000 — as bad as no second operator at all.
  Also: `F_l^cont(Delta)` is analytically divergent at `Delta=1` for
  every l (a real 2D-CFT feature — a dimension-1 operator's chordal
  power law is marginally log-divergent when integrated over the full
  sphere), so epsilon structurally cannot enter this decomposition as a
  simple second term regardless of amplitude. Since the c=1/2 minimal
  model has no other primaries besides `1`, `sigma`, `epsilon`, the
  freely-fit `Delta_2~0.63-0.66` is most likely a **non-universal
  lattice artifact**, not a second physical CFT operator. See
  `journal.md`'s "REJECTED" entry for the two tests and reasoning.
- **(2026-08-25, same session) BREAKTHROUGH — RESOLVED: the whole
  `Delta_l_pair`/two-operator anomaly is an artifact of the l-space
  (Legendre/recursion) analysis method, not the physics.** Per user
  direction, built `src/real_space_2pt_test.cc` (new standalone driver,
  bypasses `S_lm`/spherical harmonics entirely — direct angular-bin
  `<s(theta)>` measurement) and `scripts/fit_real_space_2pt.py`. Fitting
  the real-space correlator directly to `A*(2-2cos(theta))^{-Delta}`,
  excluding only the smallest angular separations (`theta_min>=0.3`),
  gives **Delta = 0.128, chi2/dof = 0.02-0.03** — close to the exact 1/8,
  robust across 8 window choices and confirmed at two independent
  resolutions (n_refine=16 and 32). This is night-and-day different from
  every l-space result. Resumming the harmonic-MC F_l data into a
  truncated real-space curve and comparing pointwise to the direct
  measurement shows a genuine, non-constant shape distortion (not just
  the known l-independent normalization factor) — likely because every
  Legendre coefficient F_l integrates over the *whole* sphere including
  the short-distance (small-theta) region where ordinary lattice
  discretization artifacts live, which a real-space fit can simply
  exclude with a window but no F_l-based estimator can. **Practical
  conclusion: extract Delta_sigma via direct real-space fits going
  forward, not via `Delta_l_pair`/the two-operator F_l fit** — those
  should be considered unreliable Delta estimators (though `symmetry_test.py`'s
  `R_l`/`kappa2_r` off-diagonal-vs-diagonal spherical-symmetry checks are
  a different question and unaffected). Next steps: continuum-extrapolate
  `Delta(a)` from real-space fits across a resolution ladder; test the
  short-distance-contamination mechanism directly via a windowed Legendre
  transform of the fine-grained real-space data. See `journal.md`'s
  "BREAKTHROUGH" entry for full detail.
  **Update (2026-08-25, later same session)**: the windowed-Legendre-
  transform idea was tried three ways (hard cutoff, smooth taper, exact
  pair-restriction via new `src/save_configs.cc`/
  `scripts/analyze_saved_configs.py`) and conclusively **fails for a
  structural mathematical reason** (windowed subtraction cannot obey the
  full-sphere Legendre recursion) — see `journal.md`'s "windowed-Legendre-
  integral attempt... FAILS" and "SETTLES the question" entries. Not worth
  revisiting. The continuum-extrapolation half of the next step was then
  started: user chose full multi-shard jackknife error bars (over the
  cheaper single-run-per-point option), so built
  `cluster/sge/real_space_shard_task.sh`/`submit_real_space_ladder.sh`
  (each independent-seed shard is its own jackknife unit, delete-one-shard
  resampling, since `real_space_2pt_test` has no in-C++ jackknife
  accumulator) and `scripts/fit_real_space_ladder.py` (jackknifes
  `Delta(n_refine)` per ladder point, then fits `Delta_inf +
  c*n_refine^-p` across the ladder). **Job submitted**: SGE job-array
  7311120, ladder `n_refine={8,16,24,32,48,64,96,128}`, naive mesh,
  `exact_sinh`, 24 shards/point, 192 tasks, `h_rt=2:00:00`. See
  `journal.md`'s "submitted first production-scale real-space Delta(a)
  ladder" entry for full design/verification detail, and its "SGE
  behavior" note about `qstat` unreliably reporting an array as finished
  while most tasks are still genuinely running — corroborate with `qstat
  -j <jobid>` before trusting completion.

  **Update (2026-08-25, later same session)**: per user idea, tested
  icosahedral-orbit averaging as a variance-reduction/contamination-check
  technique (`src/test_orbit_avg.cc` smoke test, then built into
  `real_space_2pt_test.cc` as `--orbit_avg`, per user direction "build it
  in next" — new `SpatialHash` class for O(1) position lookup at
  production `n_refine`). Confirmed unbiased (~4-5x free error reduction
  for a single fixed pair) and, in a same-seed direct comparison against
  the non-orbit z-binned method, found the two agree everywhere except the
  single shortest-distance bin (6 sigma disagreement, already excluded by
  the existing `theta_min=0.3` fit window) — real evidence of z-binning's
  cross-orbit contamination, but not one that should change the ladder's
  `Delta_inf` result given the existing exclusion window. Submitted a
  second, orbit-averaged production ladder for direct comparison:
  `cluster/sge/submit_real_space_ladder_orbit.sh`, job-array 7311202, same
  ladder/coupling/stats as job 7311120 (`SEED_BASE=810000`, disjoint).
  **Both jobs submitted, neither yet analyzed as of this update — next
  session**: confirm both complete (`qstat -j 7311120`/`7311202`), then run
  `scripts/fit_real_space_ladder.py` against both
  `campaign_runs/real_space_ladder_2026-08-25/` and
  `campaign_runs/real_space_ladder_orbit_2026-08-25/`, compare `Delta_inf`
  between the two methods, and read off the result vs. the exact 1/8. Full
  detail: `journal.md`'s "icosahedral-orbit-averaging test" and "built
  --orbit_avg into real_space_2pt_test.cc" entries; `CLAUDE.md`'s new
  "Real-space Delta_sigma extraction" section has the consolidated
  reference for all of this method's files/flags.
- **(2026-08-25)** `q=5` meshes have exact icosahedral (not full SO(3))
  symmetry, which has its own invariant spherical-harmonic subspace
  starting at l=6 (multiplicity 1, matching `OptimizeIntegrator`'s own
  `n_l` formula) — confirmed via a deterministic free-scalar FEM test
  (`src/fem_scalar_test.cc`) that l=6's off-diagonal `R_l` never shrinks
  under refinement or mesh relaxation, unlike every neighboring l which
  converges cleanly. l=6 (and, at higher `l_max`, l=10/12/15/...) should be
  excluded from any "does SO(3)-breaking shrink under refinement" claim.
  Not yet applied to `symmetry_test.py`/`run_cft_ladder.py`'s default
  output — open item. Does not explain the `Delta_l_pair` split above (l=6
  is not an outlier in that diagnostic).

## Phases

**Phase 0: scaffolding.** Done — build (`cluster/build.sh`), the
`--coupling_rule`-gated driver, the full `M_l[m,m']` measurement, the
block-jackknife accumulator, and `scripts/symmetry_test.py` all exist and
are smoke-tested (mechanism only, no production statistics). See
`journal.md` 2026-08-22 entries.

**Phase 1: spherical-symmetry production plan (designed 2026-08-22;
executing since 2026-08-23 — see status note after step 7).** Scope:
campaign objective (1), **verify spherical symmetry**, only. Objective (2),
CFT-data extraction (Delta_sigma/Delta_epsilon), is explicitly **out of
scope** for this phase — it needs the energy-operator definition and the
`C_l(l,Delta)` closed-form derivation, neither of which the spin-only
symmetry test requires. Both stay open items (see above) but do not block
Phase 1 from closing.

1. **Coupling**: `--coupling_rule exact_sinh` for the production
   conclusion; `duality` run in parallel at the same resolutions as a free,
   non-blocking comparison. See "Critical coupling" above. **Status**: only
   `exact_sinh` has been run at production scale so far
   (`production_2026-08-23`, `push5_2026-08-24_lmax12`); `duality` has only
   been smoke-tested (`campaign_runs/smoke_test_2026-08-22/`), never run as
   a parallel production comparison — still open.

2. **Mesh resolution ladder**: `q=5`, `n_refine in {4, 8, 16, 32}`
   (`n_sites = 10*n_refine^2 + 2` → 162 / 642 / 2562 / 10242 sites, formula
   confirmed against the `n_refine=1,2` examples in `journal.md`). Four
   points, enough to see a refinement trend without committing to a large
   grid before the method is validated. **Status**: production ladders were
   extended past this plan to `{4,8,16,32,64,128}` (`production_2026-08-23`,
   `l_max=8`) and `{8,16,32,64,128}` (`push5_2026-08-24_lmax12`, `l_max=12`,
   `n_refine=4` dropped — see step 3).

3. **`l_max` per resolution — via an explicit aliasing check, not the
   unvalidated `0.3*sqrt(n_sites)` rule of thumb above**: at each mesh
   resolution, run short/low-stats diagnostics pushing `l_max` upward until
   `R_l`/chi2 at the highest `l` blow up unphysically (an aliasing
   signature) or `QfeLatticeS2::OptimizeIntegrator` stops converging.
   Record the safe `l_max` per resolution. Production then uses a single
   shared `l_max` — the smallest safe value across the whole ladder — so
   the same harmonic levels are compared at every resolution. **Status**:
   `production_2026-08-23` used `l_max=8` as a **user override, not derived
   from this check** (flagged in `journal.md` 2026-08-23 — its
   `n_refine=4` point is a plausible aliasing risk by extrapolation, never
   separately confirmed). The check itself was later done, but only for
   `l_max=12` (`push4_2026-08-22_lmax12`, run 2026-08-23): safe threshold
   found to be `n_refine>=8` (`n_refine=2,4` contaminated for `l>=5`,
   `kappa2_r` up to O(10)). `push5_2026-08-24_lmax12` applied that finding
   by dropping `n_refine=4` from its ladder. No aliasing check has been
   done at `l_max=8` specifically.

4. **`jack_block_size` — via a post-hoc block-merge scan, not driver
   reruns**: stored per-block partial sums in
   `*_ylm_2pt_full_jackblocks_*.dat` are additive, so a single
   `jack_block_size=1` run lets an offline script merge every `k` native
   blocks into a superblock and recompute the jackknife error as a function
   of `k`. Add a `--merge-blocks` sweep mode to `scripts/symmetry_test.py`;
   pick the smallest `k` where the jackknife error stops growing (the
   tau_int signature) as the production block size. No need to relaunch the
   driver at multiple `--jack_block_size` values. **Status**: the formal
   tau_int-based merge scan was never implemented; in practice
   `jack_block_size` has instead been picked directly at driver-invocation
   time to land the native block count in the few-hundred range (see
   `CLAUDE.md` "Cluster conventions") — `100` for `production_2026-08-23`
   (`n_traj=20000` → 200 blocks), `4000` for `push5_2026-08-24_lmax12`
   (`n_traj=1,000,000` → 250 blocks, chosen to land exactly on step 5's
   gate at `l_max=12`). This satisfies step 5's numeric gate without ever
   running the tau_int check this step called for.

   **Closed by argument, not by running the scan (2026-08-25, user
   direction)**: `tau_int ~ O(1)` in units of cluster sweeps is expected
   for Wolff cluster updates near criticality — the algorithm exists
   specifically to eliminate the critical slowing down that makes
   `tau_int` diverge for local (Metropolis-only) updates. Every push in
   this campaign includes `n_wolff=5` cluster updates per trajectory
   (`ising_s2_crit.cc`), so native block count and effective-independent
   sample count should closely track each other here — the formal
   `--merge-blocks` sweep was never built/run, but the theoretical
   justification is now recorded (`journal.md` 2026-08-25 "push7
   resubmission fix" entry) as the reason existing and future jackknife
   error bars can be read at face value rather than as an optimistic
   lower bound.

5. **Statistics target**: production gate is `n_blocks >= 10*(2*l_max+1)`
   at the largest `l` measured (double `scripts/symmetry_test.py`'s
   existing `--min-blocks-factor` floor of 5 — that floor keeps the
   jackknife covariance from being nonsense, this gate asks for enough
   margin to draw an actual conclusion). Size `n_traj` after step 4 fixes
   `jack_block_size`. **Status**: nominally met by both production pushes
   (200 blocks >= 170 needed at `l_max=8`; 250 blocks >= 250 needed at
   `l_max=12`, exact) under the block-count arithmetic above — but see step
   4's caveat that this counts native blocks, not tau_int-verified
   independent samples.

6. **Gate / success criteria**: the primary measure is the **reduced
   spectral cumulants** `kappa2_r, kappa3_r, kappa4_r` of each `M_l`'s own
   eigenvalue distribution (`scripts/symmetry_test.py`, added 2026-08-22),
   built from `Tr(M_l^k)` for `k=1..4` via the standard moment-to-cumulant
   formulas and normalized by `kappa1^k`. These are basis-independent —
   unlike `R_l`/diagonal-chi2, which depend on the mesh's arbitrary
   quantization axis and so can't distinguish genuine SO(3) breaking from
   an accidental (mis)alignment of that axis with a mesh symmetry
   direction. Under exact SO(3) symmetry `M_l` is proportional to the
   identity, so every eigenvalue is equal and `kappa2_r=kappa3_r=kappa4_r=0`
   exactly. Spherical symmetry is considered verified for a resolution
   ladder if, at each fixed `l`, `kappa2_r/kappa3_r/kappa4_r` do not grow
   (ideally shrink toward zero) as `n_refine` increases, land within the
   step-5 statistics target at the finest resolution tested, and `U4`
   stabilizes across resolutions (a scale-invariance sanity check standing
   in for the crossing-point check that isn't available without a free
   coupling parameter — see "Critical coupling" above). `R_l`/diagonal-chi2
   are kept as a secondary diagnostic — they test whether the *specific*
   mesh-fixed axis looks anomalous, a genuinely different question from
   whether `M_l` is symmetric at all. **Status (2026-08-24, evaluated)**:
   run against all completed production data (`production_2026-08-23`
   l_max=8, full ladder; `push5`+`push6` l_max=12, full ladder) — see
   `journal.md` 2026-08-24 "ran PLAN.md step-6 gate" entry and
   `campaign_runs/analysis_lmax12_2026-08-24/`. **Gate not passed**: l=1,2
   are consistent with zero / no trend (passes) but l>=3 shows a real,
   multi-sigma, non-vanishing `kappa2_r` plateau (order 1e-5 to 1e-4) that
   does not shrink with refinement out to `n_refine=128` at either l_max —
   at l_max=8 it actually rises monotonically from `n_refine=16` to `128`.
   This is now the central open question for Phase 1, not a statistics
   gap. `U4` stabilization (the gate's other half) and `kappa3_r/kappa4_r`
   were not separately analyzed yet.

7. **Execution shape**: one SGE array task per `(n_refine, coupling_rule,
   seed)` tuple under `cluster/sge/`, following the provenance-first task
   design already used by `Twist`/`Twist_xR` (`/projectnb/qfe/misra/CLAUDE.md`)
   — task identity fixed by the tuple, not by worker count or completion
   order. **Status**: built and used twice. `cluster/sge/production_shard_task.sh`
   is the shared, fully-parameterized per-task driver invocation;
   `cluster/sge/submit_production.sh` and `cluster/sge/submit_push5.sh` are
   per-push submit scripts, each pinning that push's ladder/`l_max`/stats/
   `jack_block_size`/`h_rt` and writing a `manifest.json` into
   `campaign_runs/<job_tag>/` before submitting (task index = shard within
   ladder point, not `coupling_rule` — every production task so far has
   used a single fixed `coupling_rule` per push, not the per-task axis
   originally described here).

**Status summary (2026-08-24)**: two production pushes have run end to end
— `production_2026-08-23` (`l_max=8`, ladder to `n_refine=128`, complete,
96/96 tasks, no errors) and `push5_2026-08-24_lmax12` (`l_max=12`, ladder
`{8,16,32,64,128}`, 1,000,000 traj/shard, in progress as job-array 7301188).
Despite that, **Phase 1 has not actually closed**: no run has been
evaluated against step 6's gate, the step-4 tau_int block-merge check was
never built (block size picked by a size heuristic instead), `l_max=8` was
never aliasing-checked, and the `duality` comparison (step 1) has never
been run past a smoke test. Next real step once `push5` output exists:
run `symmetry_test.py`'s gate logic across the full `l=1..12` range on both
production datasets, not just add more raw statistics.

**Status update (2026-08-25)**: `push7_2026-08-25_dual_naive_eqarea`
submitted — the first production-scale head-to-head of `naive` vs.
`equal_area` `--mesh_mode` (the two surviving candidates from the
2026-08-24 investigation below), ladder `{2,3,4,6,8,12,16,24,32,48,64,96,
128}`, `l_max=8`, ~10M pooled trajectories split across 24
`(n_refine,mesh_mode)` combinations. See `journal.md`'s "push7" entries
(2026-08-25, several) for the full design/seed/stats rationale and two
real bugs found and fixed along the way (a login-node CPU-contention stall
in the mesh-cache prewarm step, and a filename-padding mismatch in the
prewarm-verification check) — as of the last check this session, the
prewarm stage was running/re-verifying and the two production arrays had
not yet been confirmed submitted; confirm that first next session. A
separate, unrelated `gd_stress_1024_2026-08-25` job (job-array `7305751`)
is also running concurrently — a pure mesh-geometry performance
characterization of `EqualizeFaceAreas` out to `n_refine=1024`, not a
production Ising run; see journal for details.

**Step 6 gate evaluated 2026-08-24 (same day, later): FAILS for l>=3** —
see the status note on step 6 above and `journal.md`'s "ran PLAN.md step-6
gate" entry. `kappa2_r` for `l>=3` plateaus at a nonzero, multi-sigma value
instead of shrinking with refinement on the naive mesh. This motivated a
same-day investigation into whether naive-mesh triangle-area
non-uniformity (std/mean ~0.13, not shrinking with resolution) is the
systematic behind the plateau — see `CLAUDE.md`'s new "Mesh construction
modes" section and `journal.md`'s full 2026-08-24 mesh-mode entries for
the detailed story. Summary: three `--mesh_mode` alternatives to `naive`
were built (`equal_face_area`, `equal_area_analytic`; a third,
`equal_dual_area`, was tried, found to diverge, and dropped per user
direction). Both surviving alternatives reduce area non-uniformity at
some resolutions but **neither is currently production-ready**:
`equal_face_area`'s fixed iteration/step schedule diverges at
`n_refine>=32`, and `equal_area_analytic` — after two real bugs were found
and fixed to get it to retain exact icosahedral symmetry (verified to
machine precision via `src/test_ico_symmetry.cc`) — currently has *worse*
area uniformity than naive due to an overly-conservative symmetry-safety
fallback. **No mesh-mode alternative has yet been run at production scale
against the step-6 gate** — the original question (does a more uniform
mesh fix the l>=3 plateau) is still open pending either fix.

**ROOT CAUSE FOUND 2026-08-25 (superseding the paragraph below — kept for
history): `Delta_s(a)` plateau is a mesh-uniformity problem, not a
coupling-formula bug.** `scripts/cft_symmetry_test.py`/
`scripts/run_cft_ladder.py` (new, reuse existing jackblocks — no new
production run needed) implement the dimension-independent symmetry test
and direct `Delta_s(a) = 2*F_1/(F_1+F_0)` scaling-dimension estimator from
`reference/Owen_section_D.tex`. Found `Delta_s(a)` plateaus at ~0.33-0.34
(exact value 0.125) on push7's full ladder and independently on
`production_2026-08-23`, mesh-mode-independent (naive/equal_area agree).
Tested all coupling rules this campaign has (`duality`, `exact_sinh`, and
a new `owen_dual` implementing arXiv:2407.00459 Eq 15/37 literally) — all
four give the same plateau, ruling out the coupling formula. Read the
actual source paper (arXiv:2407.00459, `reference/owen_2407.00459.pdf`):
their coupling is only exactly critical when every triangle has **equal
circumradius AND equal perimeter** (Eq 33-37), not equal area. Their own
"basic"/unmodified mesh (= this campaign's `naive`) shows the identical
documented failure, and they explicitly tested an equal-*area* mesh
smoother (same style as this campaign's `equal_area`) and report it also
fails, for the same reason. Fixing this requires their Appendix B.1/B.2
circumradius+perimeter mesh smoother (icosahedral-orbit-constrained
Newton's method minimizing `E_R+E_P`, Eq 39) — a new, substantial
implementation this campaign has never built (only area-based smoothers
exist so far). Full chain of reasoning, dead ends, and exact paper
quotes: `journal.md` 2026-08-25 "ROOT CAUSE FOUND" entry. **This is now
the concrete, well-evidenced next step for the campaign** — not another
area-based mesh variant, and not further coupling-rule tuning.

**Objective (2) (Delta_sigma extraction) CLOSED 2026-08-25/26: clean
confirmation of Delta_sigma = 1/8.** The real-space method's two
production ladders (z-binned job 7311120, orbit-averaged job 7311202,
both 192/192 tasks clean) continuum-extrapolate to `Delta_inf =
0.12497(20)` (-0.14 sigma) and `0.12515(23)` (+0.65 sigma) respectively
against the exact 1/8, **once the extrapolation ansatz is corrected to a
linear `Delta(a) = Delta_inf + c/n_refine` leading correction** (chi2/dof
0.82/0.97) rather than the previously-assumed quadratic one (chi2/dof
3.06/1.60, ~4.4 sigma off — a bad-ansatz artifact, not a real
discrepancy). Two independent binning methods agree at every ladder
point. See `journal.md`'s "CLEAN CONFIRMATION" entry (2026-08-25/26) for
the full table and `campaign_runs/real_space_ladder_2026-08-25/plots/
delta_vs_invL_power1.png` for the plot. **Objective (1) (spherical-
symmetry verification, `R_l`/`kappa2_r` gate) remains open** — see the
step-6 gate-fail note above, unaffected by this result.

**Objective (1) update, 2026-08-28: naive-vs-equal_area real `kappa2_r`
comparison, n_refine<=48 (jobs 7341057/7341058, still running past this
point).** Real off-diagonal `kappa2_r` (not the diagonal proxy) confirms
`equal_area` has a systematically lower l>=3 plateau than `naive`
(roughly 2.5-4x smaller, 4-7 sigma significant, from n_refine=12 to 48) —
this **reverses** the push7-session conclusion above that both mesh modes
"converge to the same plateau by n_refine>=12-16" (that was based on
lower-statistics spot checks). However, **neither mesh mode's own plateau
is shrinking with resolution** over this range — both wobble within a
factor of ~2 but show no continuum-limit-consistent downward trend from
n_refine=12 to 48. So mesh non-uniformity (area, at least) is not a red
herring — it measurably helps — but it does not fully explain the l>=3
symmetry breaking on its own. `equal_rp` (circumradius+perimeter) or a
genuine residual lattice-artifact explanation remain live candidates.
Full table: `journal.md`'s 2026-08-28 entry. Next: rerun at n_refine up
to 128 (and eventually 512, once job 7342648's chain runs) once those
jobs finish, to see whether either mode's plateau eventually shrinks or
both saturate at a nonzero floor.

**Objective (1) update, 2026-08-29: naive-vs-eqarea extended to n_refine
up to 112 (confirmed clean, eqarea 30-40x lower than naive at l=3, gap
widening not narrowing with resolution — neither plateau shrinking
toward zero); third mesh mode (`equal_rp`) data now on disk but not yet
compared.** `equal_rp`'s 100-shard production (job 7358337, the
`lean_mesh_compare_100shards_2026-08-27` ladder's third arm) finished
clean. `scripts/compare_three_mesh_kappa2r.py` + `cluster/sge/
compare_three_mesh_task.sh` (new) were built to finally run the 3-way
`kappa2_r` comparison this section has been waiting on since `equal_rp`
was introduced 2026-08-25 — submitted as **job 7363395**, queued as of
this entry, **no numbers yet**. This is the concrete next check for
whether `equal_rp` (built specifically because arXiv:2407.00459 predicts
equal-area alone is insufficient) shrinks the l>=3 plateau further than
`equal_area` already does, or whether it plateaus at the same nonzero
floor. Check `qacct -j 7363395` before citing any equal_rp kappa2_r
number.

Separately, the `lean_ladder_512_2026-08-27` production (naive+eqarea,
the eventual n_refine up to 512 point for this same gate) had stalled on
a stuck `-sync y` recovery script (see `journal.md`'s 2026-08-29 entry)
— fixed and resubmitted as jobs 7363384/7363385, also queued, not
started.
