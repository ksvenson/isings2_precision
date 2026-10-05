# Journal

Running log of what has actually been done in this campaign, including dead
ends and reversed decisions. Newest entries at the bottom.

## 2026-08-21 — campaign init

Scaffolded `CLAUDE.md`/`PLAN.md`/`directive.md`/`journal.md` and an empty
`src/`/`scripts/`/`cluster/sge/`/`configs/`/`reference/`/`campaign_runs/`
directory skeleton. Initialized as its own git repo (sibling to
`../IsingS2/`, not nested inside it).

Objective (from user): verify spherical symmetry and extract CFT scaling
dimensions for the critical 2D Ising model simulated on S^2, via spherical
harmonic decomposition of the two-point function — build `S_lm = sum_i
sqrt(g_i) Y_lm(x_i) s_i` projections per thermalized configuration, form the
`(2l+1)x(2l+1)` matrix `<S_lm S*_lm'>` at each harmonic level `l`, and check
both (a) it reproduces the SO(3)-diagonal structure expected on the sphere
(symmetry test) and (b) its eigenvalues match the CFT two-point function
prediction for some `Delta` (CFT-data test).

Surveyed `../IsingS2/` for reusable engine code:
- `src/isingS2.cc` — file-driven mesh + Ising MC driver; measurement loop is
  currently minimal (only prints the action), no `S_lm` recording yet.
- `include/S2.h` (`QfeLatticeS2`) — already has `OptimizeIntegrator(l_max)`
  (solves for per-site quadrature weights `sites[s].wt` that integrate
  spherical harmonics exactly up to `l_max`) and `UpdateYlm`/`GetYlm`
  (cached `Y_lm` values per site). This is exactly the machinery the `S_lm`
  projection needs — confirms build-on-top-of rather than reimplement.
- `include/ising.h` (`QfeIsing`) — spin field, `Metropolis`/`WolffUpdate`,
  field checkpointing.
- `src/Makefile` is Homebrew/macOS-flavored (`/opt/homebrew/...` include
  paths) — will not build on SCC as-is; noted as a Phase-1 prerequisite in
  `PLAN.md`.

Did not derive or hardcode the closed-form harmonic-eigenvalue-vs-Delta
formula yet (deliberately — see `directive.md`: it must be symbolically
derived/verified before use, not pulled from memory). Did not start any
build, mesh generation, or simulation. PLAN.md Phase 1 open items (energy
operator definition, literature reference for the eigenvalue formula,
statistics policy, SCC build, C++-vs-Python split for the harmonic
projection) are unresolved — next session should start there.

## 2026-08-21 (cont): imported Twist's coupling-angle maps

Per user direction, loaded both of `../Twist/`'s flat-triangle
coupling-to-angle maps into `reference/twist_coupling_maps/`:
- `coupling_from_angle.py` — exact rule `K(theta)=0.5*asinh(cot(theta))`
  and empirical `kappa(s,a)` monomial-fit loader, with `s_a_from_angles()`
  matching Twist's `delta_i = theta_i - pi/3`, `s_i = -delta_i`,
  `a_i = delta_j - delta_k` convention. Confirmed via subagent
  investigation (grepped `Twist/scripts/fit_permutation_symmetric_coupling_map.py`,
  `compute_geometry_boundary_angles.py`, journal 2026-08-11 entries).
  Sanity check: `exact_sinh_rule(pi/3)=0.27465` vs accepted degree-2
  empirical `kappa(0,0)=0.27570` — ~0.4% apart, consistent with Twist's
  reported ~0.1-0.3% agreement over its full fitted domain.
- `kappa_sa_fit.json` — static copy of Twist's accepted Phase 6 production
  fit (`campaign_runs/phase6_coupling_map/kappa_sa_fit.json` as of
  2026-08-21). Not live-synced; re-copy deliberately if Twist's fit changes.

Flagged, not yet resolved: both maps are derived for the flat 2D
triangular lattice (three coupling directions per unit cell); `IsingS2`'s
engine already assigns each S2-mesh link its own geometric weight via
`QfeLatticeS2::UpdateWeights()` (`link->wt`, entering the action as
`beta*link->wt*spin_i*spin_j`, see `include/ising.h`/`S2.h`) — a
cotangent/dual-area Regge-type weight, not literally this flat-triangle
formula. How/whether the imported maps interact with or replace that
existing weight on the (mostly-flat, curvature concentrated at vertices)
S2 mesh triangles is unresolved.

**User direction**: run the first precision production pass using the
**exact sinh rule** (not the empirical rule) to assign local per-edge S2
couplings. Only if spherical symmetry and CFT scaling dimensions check out
against the exact rule does a second, separate pass verify the empirical
`kappa(s,a)` rule too. See `PLAN.md` "Local coupling assignment" section.
Next session: resolve the exact-rule-vs-`link->wt` interaction, then set
up the exact-rule precision run.

## 2026-08-22 — SCC build achieved; discovered `ising_s2_crit.cc` already does most of Phase 1's measurement

**Build.** `src/isingS2.cc` (the driver `CLAUDE.md` originally pointed at) is
stale against current headers — it uses the base `QfeLattice` type instead
of `QfeLatticeS2`, so it's missing `ReadOrbits`/`UpdateWeights`/`.r`/etc.
that only exist on the S2 subclass; not just "needs a build flag fix."
`src/ising_s2_crit.cc` (the K_c-search driver PLAN.md already pointed to
for K_c-per-mesh-resolution) is the one that's actually current and
correct — it builds and runs cleanly against `QfeLatticeS2`.

Fixed one small pre-existing engine bug in `../IsingS2` to get it building:
`ising_s2_crit.cc:145` called `QfeRng::ReadRng()`, a method that no longer
exists on `QfeRng` (`include/rng.h` — RNG-state checkpoint/resume looks like
a dropped feature). Sibling drivers `2ising_mesh_crit.cc`/
`3ising_mesh_crit.cc` hit the same issue and just comment the call out;
mirrored that exact pattern (commented, did not delete/reimplement) rather
than inventing RNG serialization. One-line change, in `../IsingS2/src/`,
not vendored here.

Working build recipe (now `cluster/build.sh`): vendored Eigen from
`../IsingS2/include/Eigen` needs no module; boost (`S2.h` needs
`boost/math/special_functions/spherical_harmonic` for `Y_lm`) comes from
`module load boost/1.83.0` (`$SCC_BOOST_INCLUDE`); `-DGRP_DIR` must point
at `../IsingS2/grp`. Compiles clean with system `g++ 8.5.0`, no `module
load gcc` needed. Binary lands at `bin/ising_s2_crit` (this campaign's own
copy, `../IsingS2/` untouched except the one-line fix above).

Smoke-tested end to end: `q=5, n_refine=1` (12-site mesh), 50 therm + 200
traj, cold start → ran, wrote `*_bulk_*.dat` (action/mag/mag^2.../U4/susc),
plus **`*_legendre_2pt_*.dat`** and **`*_ylm_2pt_*.dat`**. Note: the
`data_dir/run_id/` output directory must exist beforehand — the driver does
not `mkdir -p` it (no `MakePath` available currently); an SGE task script
will need to create it before invoking the binary.

**Major finding — the measurement side of Phase 1 is mostly already
written**, not a from-scratch task:

- `ylm_2pt` (`ising_s2_crit.cc` ~L253-426): computes
  `ylm(s,y_i) = lattice.CalcYlm(s, l, m)` at every site (via `QfeLatticeS2`,
  for `m >= 0` only — `y_i` indexes `(l,m)` pairs with `0 <= m <= l`, count
  `n_ylm = (l_max+1)(l_max+2)/2`), then per configuration builds
  `S_lm = sum_s wt_s * Ylm(s) * spin_s` (site quadrature weight `wt`, same
  object PLAN.md's `S_lm` formula calls for) and measures
  `|S_lm|^2 / vol^2`. This *is* the diagonal (`m=m'`) of PLAN.md's `M_l`
  matrix. **It does not currently measure off-diagonal `m != m'` terms**, so
  it cannot do the symmetry test (`M_l` proportional to identity) as-is —
  only ever writes/reads diagonal entries. Extending it to the full
  `(2l+1)x(2l+1)` matrix (or at minimum the `m,m'` cross terms) is now a
  much smaller, well-scoped addition rather than new infrastructure.
- `legendre_2pt` folds the `m` values at each `l` into one number
  (`legendre_2pt_sum * coeff`, weighting `m=0` by 1 and `m>0` by 2) — this
  looks like an m-averaged/summed `C_l` proxy, not literally `C_l` from the
  addition theorem; needs reading `coeff`'s definition before trusting it
  as the CFT-test observable (not yet done this session).
- **Local coupling is already non-default and looks like a duality-based
  formula, not `link->wt` from `UpdateWeights()`**: right before the Ylm
  setup (`ising_s2_crit.cc` ~L204-250), every link's `wt` is overwritten via
  `K = 0.5*asinh(1.0 / sinh(2.0*L))`, where `L = atanh(sqrt(...))` is built
  from a cross-ratio of face-angle cosines and edge lengths around each
  link (a star-triangle/Regge-type duality relation). This is **not** the
  same formula as the imported exact rule
  (`K(theta) = 0.5*asinh(cot(theta))`, `reference/twist_coupling_maps/`) —
  different functional form, and it's already active in the driver, not
  optional. Reconciling/comparing this pre-existing formula against the
  Twist exact-rule import is now the concrete content of the still-open
  "exact-rule-vs-`link->wt`" item — not a vs-`UpdateWeights()` comparison
  as PLAN.md currently frames it. Have not yet worked out whether the two
  formulas are equivalent under a change of variables or are genuinely
  different physics; do that before the coupling-assignment open item is
  called resolved.

Did not touch statistics policy, energy-operator definition, or the
CFT-eigenvalue derivation this session — build-only, per user direction.
Next session: pick up the coupling-assignment reconciliation (now better
scoped, see above), or extend `ylm_2pt` to off-diagonal `m,m'` terms.

## 2026-08-22 (cont): added exact-sinh coupling rule + full M_l[m,m'] measurement, smoke-tested both

Per user direction, extended `../IsingS2/src/ising_s2_crit.cc` (same driver
as above, one-line-fix already applied) with two additions, both in that
file only:

**1. `--coupling_rule={duality,exact_sinh}` CLI flag** (default `duality`,
preserves prior behavior exactly). `exact_sinh` applies Twist's flat-triangle
star-triangle rule `K(theta)=0.5*asinh(cot(theta))` per link: for each of
the link's two adjacent faces, get that face's angle opposite the link via
law of cosines from the three chord lengths already computed in the
existing duality block (`len_l0/l1/l2` for face `f0`, `len_l0/l3/l4` for
`f1`), apply the exact rule to each face's angle, and average the two
resulting `K`s. Averaging-over-both-faces is a judgment call, not a derived
result — flagged explicitly as such in the code comment; the S2 mesh has no
single-triangle-per-edge convention the way Twist's flat lattice does
(translational symmetry there makes the two faces bordering an edge
equivalent, so no averaging is needed; not true here). This is the concrete
open item from the previous entry, not fully closed — it's now a testable
implementation, not just an analysis question.

Validated by comparing `K_duality` vs `K_exact_sinh` per-link (temporary
debug printf, reverted): on `q=5, n_refine=1` (12 sites, the bare
icosahedron — every face is an exactly congruent equilateral triangle, no
curvature-defect asymmetry) **the two formulas give bit-identical
`K=0.2746530722` on every link**, matching the known isotropic value
`0.5*asinh(1/sqrt(3))`. This is a real mathematical coincidence at maximal
symmetry, not a bug — confirmed by also running `q=5, n_refine=2` (42
sites, curvature defects present, faces no longer all congruent), where the
two rules genuinely diverge: mean action `-0.509(25)` (duality) vs
`-0.467(22)` (exact_sinh) at matched seed/statistics. So the flag is wired
correctly and the two formulas are only degenerate in the fully-symmetric
limit, as expected.

**2. Full `(2l+1)x(2l+1)` `M_l[m,m'] = <S_lm S*_lm'>` measurement**
(new output file `*_ylm_2pt_full_*.dat`, additive — old `*_ylm_2pt_*.dat`
diagonal-only file and format untouched). Reconstructs `S_l,-m` from the
existing `S_l,m` (`m>=0`) via the same `(-1)^m * conj()` relation `GetYlm`
already uses for `Y_l,-m` (`S2.h`), so no new per-site Ylm evaluation is
needed. Stores only the upper triangle `m<=m'` per `l` (Hermitian, so the
rest is redundant) via a precomputed `(l,m,m')` list built once outside the
trajectory loop. Read/write checkpoint support mirrors the existing
`ylm_2pt` pattern (reconstructs `QfeMeasReal` internal `sum`/`sum2`/`n` from
written mean/error/n, same trick `QfeMeasReal::ReadMeasurement` uses) but is
hand-rolled in `ising_s2_crit.cc` rather than added generically to
`QfeMeasComplex` in `statistics.h`, to keep the shared header untouched.

Smoke-tested (`q=5, n_refine=2`, `l_max=4`, 200 traj, both coupling rules):
diagonal entries in the new file match the old diagonal-only file exactly
(e.g. `l=1,m=1`: `3.8671768722688461e-03` in both), and off-diagonal terms
(e.g. `l=1,m=-1,m'=0`) come out nonzero but visibly smaller than the
diagonal — physically sane for a coarse (42-site) mesh at low statistics
(`n=40` samples), not yet a meaningful symmetry-test result. Smoke outputs
kept at `campaign_runs/smoke_test_2026-08-22/{duality,exact_sinh}/` for
reference only — informal statistics, not production data.

**Not done / still open**: this smoke test does not constitute the actual
symmetry test or CFT-data fit — statistics policy, error estimation
(jackknife/bootstrap over configs, autocorrelation-aware), K_c location per
mesh resolution, and the CFT eigenvalue closed-form derivation are all
still open per PLAN.md. The exact-sinh-vs-duality reconciliation is also
still open in the sense that "average the two adjacent faces" has not been
justified from first principles — it's a reasonable default, not a decided
convention.

## 2026-08-22 (cont): CLAUDE.md refreshed to match build/driver reality

`CLAUDE.md` still described the build as "not yet done" and pointed at
`isingS2.cc` as the driver, both stale after this session's build/measurement
work above. Rewrote it: build section now documents the working
`cluster/build.sh` recipe and the `../IsingS2/` RNG fix; driver section
points at `ising_s2_crit.cc` with its actual CLI flags, the
`--coupling_rule` split, the output-dir-must-preexist gotcha, and the four
output file types (flagging `*_legendre_2pt_*.dat` as an unverified `C_l`
proxy, matching the open item above). No code changes.

Next: designing the actual spherical-symmetry test statistics (see below).

## 2026-08-22 (cont): designed and implemented block-jackknife accumulator for the symmetry test

**Design discussion.** The symmetry test needs error bars on quantities that
combine multiple `M_l[m,m']` entries per `l` — an off-diagonal power ratio
`R_l = Σ_{m≠m'}|M_l[m,m']|² / Σ_m|M_l[m,m]|²`, and a chi-square of the
diagonal entries against their mean using the cross-`m` jackknife
covariance (not a plain coefficient of variation, which would ignore
correlations between different `m`). Both need per-sample resampling
(jackknife/bootstrap), but the driver only ever accumulated per-entry
mean/error/n online — no way to reconstruct cross-entry correlations from
that. Considered two fixes: (a) dump raw per-config `S_lm` to disk and do
all resampling offline in Python, or (b) accumulate block-jackknife partial
sums online in C++, small enough to stay lightweight but still
block-synchronized across all `(l,m,m')` pairs so cross-entry correlations
survive. **User direction: option (b).**

User also pointed out that with Wolff cluster updates dominating the
update mix (`n_wolff=5` vs `n_metropolis=4` per trajectory, already the
driver default), `tau_int` for these observables is plausibly O(1)
trajectory, not the O(10-100) pure-Metropolis would need — so blocks can be
small. Did not hardcode a specific size on that basis, though: at
`block_size=1`, file size approaches the raw-per-config-dump option that
was explicitly passed over (e.g. `l_max=6` has 252 `(l,m,m')` pairs; 2000
blocks x 252 pairs is already ~30MB per run). Made `jack_block_size` a CLI
flag (default 1) instead, to be set from an actual measurement (e.g.
comparing jackknife error at a few block sizes on a real production run) —
same "don't hardcode from memory, derive/verify" spirit as `directive.md`'s
closed-form-formula rule, applied to a statistics-policy number instead of
a physics formula.

**Implementation**, in `../IsingS2/src/ising_s2_crit.cc` (same driver as
above; not vendored here, per `directive.md`):

- New `--jack_block_size N` flag (default 1, must be >=1).
- `JackBlock{n, sum[full_pairs.size()]}` accumulator: each measured sample
  (post `n_therm`, post `n_skip`, same loop iteration that already updates
  `ylm_2pt_full`) adds its `Slm(l,m)*conj(Slm(l,m'))/vol^2` value into the
  current block's per-pair running sum; a block is closed and pushed once
  it hits `jack_block_size` samples, and a trailing partial block (if any)
  is flushed after the trajectory loop rather than dropped — the offline
  jackknife is expected to weight blocks by their own `n`, so a smaller
  final block needs no special-casing.
- New checkpoint file `*_ylm_2pt_full_jackblocks_*.dat`: header
  `# n_pairs=.. jack_block_size=.. n_blocks=..`, then one row per
  `(block, l, m, mp)` with the block's summed real/imag parts and its `n`.
  Read back at start (if present) and treated as already-completed blocks
  that new blocks get appended after — mirrors the existing
  `ylm_2pt_full`/`ylm_2pt` resume pattern (rerun same binary against same
  `data_dir`/`run_id`/`seed`, accumulate more, rewrite combined file at
  exit). **Caveat, not enforced by the code**: resuming with a different
  `--jack_block_size` or `--l_max` than the original run silently produces
  a mix of block sizes / stale pair ordering — keep both fixed for a given
  `data_dir`/`run_id`/`seed`.

**Smoke-tested** (`q=5, n_refine=2, l_max=3`, 50 pairs, `jack_block_size=4`,
300 traj / n_skip=5 = 60 samples): first run produced exactly the expected
15 blocks (`60/4`), file line count `50*15+1=751` matched exactly. Re-ran
the identical command against the same `data_dir`/seed to check resume:
correctly read the 15 prior blocks and appended 15 more, landing at 30
blocks / 1501 lines, same 50-pair header — resume accumulation confirmed
working, not just the fresh-run path. Deleted the scratch output after
(mechanism-only check, no physics conclusion, nothing worth keeping as
reference).

**Not done yet**: the offline analysis script (`scripts/symmetry_test.py`)
that turns these blocks into `R_l`/diagonal-chi-square jackknife estimates
is not written. Also still open: actually measuring `tau_int` on a real
production run to pick a justified `jack_block_size` (currently just the
default of 1), and everything else already listed as open in `PLAN.md`
(energy operator, `C_l(l,Delta)` literature/derivation, K_c location,
coupling-rule reconciliation).

## 2026-08-22 (cont): wrote `scripts/symmetry_test.py`, smoke-tested

Implemented the offline analysis promised above. Reads a
`*_ylm_2pt_full_jackblocks_*.dat` file, and per harmonic level `l`:

- Reconstructs the full Hermitian `(2l+1)x(2l+1)` `M_l[m,m']` matrix by
  leave-one-block-out jackknife (upper triangle from the file, lower via
  conjugate symmetry).
- **`R_l`** = relative off-diagonal power (`Σ_{m≠m'}|M_l|² / Σ_m|M_l[m,m]|²`),
  point estimate from the full dataset, error from the standard jackknife
  variance formula over leave-one-out samples.
- **Diagonal chi-square**: builds the `(2l+1)x(2l+1)` jackknife covariance
  matrix of the diagonal entries across `m`, then does a GLS fit of the
  diagonal to a constant (`c_hat`) and reports `chi2`/`dof` — this is the
  cross-`m`-correlation-aware version discussed earlier, not a plain
  coefficient-of-variation. Skips the chi2 (prints `--`) when
  `n_blocks < 5*(2l+1)` for a given `l`, since the jackknife covariance
  isn't trustworthy with too few samples relative to matrix dimension —
  `--min-blocks-factor` is a CLI flag, not hardcoded silently.

**Smoke-tested** end to end (`q=5, n_refine=2, l_max=3`, `jack_block_size=1`,
1000 traj / n_skip=2 = 500 samples/blocks, well above the `5*(2l+1)`
threshold through `l=3`): script parsed the file, ran without error, and
produced physically sane numbers — `R_l` of order 1-2% with jackknife
errors of comparable size, `chi2/dof` <= ~0.9 for `l=1..3` (consistent with
no significant SO(3) breaking detected at this low statistics on a small
42-site mesh — not a claim of a validated result, just confirms the
pipeline produces sane, not garbage, numbers). `l=0` (trivial `dof=0`,
single diagonal entry) handled correctly. Deleted the scratch run output
after (mechanism check only).

**Not done yet**: this has not been run on real production statistics, at
multiple `l_max`/mesh resolutions, or used to draw any actual symmetry
conclusion — that needs K_c located first (still open) so runs are at the
right coupling. Also still open, unchanged from above: `tau_int`
measurement to justify a production `jack_block_size`, energy operator,
`C_l(l,Delta)` derivation, coupling-rule reconciliation.

## 2026-08-22 (cont): designed the Phase 1 spherical-symmetry production plan; PLAN.md updated

**K_c framing resolved.** The previous entry's "needs K_c located first"
framing is now recognized as stale: `--coupling_rule exact_sinh` assigns
every link's `K` from local mesh geometry with no free global scale factor,
so there is nothing to locate via a Binder-cumulant crossing scan. Per user
direction, production runs take the exact_sinh coupling as given —
criticality is asserted by construction (the same local self-dual logic
behind `Twist/`'s flat-lattice result), and what the symmetry test itself
checks is whether that assumption survives assembly on a curved, defected
mesh (via `R_l`/diagonal-chi2 trends under refinement, not a separate
tuning step). This resolves PLAN.md's "Critical coupling" open item for
Phase 1 (not for objective 2 — CFT-data extraction may still need to think
about this differently, revisit if/when that phase starts). `duality`
stays as a free, non-blocking comparison run.

**Scope split made explicit**: campaign objective (1), spherical-symmetry
verification, does not need the energy-operator definition or the
`C_l(l,Delta)` closed-form (both still open) — only objective (2),
CFT-data extraction, does. PLAN.md's Phase 1 is now scoped to objective (1)
only, so it can close without waiting on either.

**Production plan written into `PLAN.md`** (new Phase 1 section):
mesh ladder `q=5, n_refine in {4,8,16,32}` (162/642/2562/10242 sites, via
`n_sites=10*n_refine^2+2` — confirmed against the `n_refine=1,2` values
already measured); `l_max` to be set via an explicit per-resolution
aliasing check rather than the unvalidated `0.3*sqrt(n_sites)` heuristic;
`jack_block_size` to be set via a post-hoc block-merge scan added to
`scripts/symmetry_test.py` (native per-block sums are additive, so this
needs one `jack_block_size=1` run per resolution, not a rerun sweep);
statistics gate `n_blocks >= 10*(2*l_max+1)`, double the script's existing
`--min-blocks-factor=5` sanity floor; success criteria are refinement
trends in `R_l`/diagonal-chi2 plus `U4` stabilization across the ladder,
standing in for a Binder-crossing check that has no coupling axis to run
along. Execution shape: one SGE array task per `(n_refine, coupling_rule,
seed)`, provenance-first per the workspace-wide `Twist`/`Twist_xR`
convention — `cluster/sge/` job script not yet written.

**Not done this session**: no code changes (no `--merge-blocks` flag added
to the script, no aliasing check run, no SGE script written), no jobs
launched. Written plan only, per user direction. Next session: pick up
step 3 (aliasing check) or step 4 (block-merge scan implementation) from
PLAN.md's new Phase 1 list — either can start independently of the other.

## 2026-08-22 (cont): added basis-independent spectral-cumulant symmetry
test (traces of M_l) to `scripts/symmetry_test.py`, smoke-tested

**Motivation.** `R_l`/diagonal-chi2 are computed in the specific `m` basis
fixed by the mesh's Y_lm quantization axis — an arbitrary choice from the
icosahedral mesh construction. A large `R_l` could mean genuine SO(3)
breaking, or just that axis being poorly aligned with an underlying
(near-)symmetry; a small `R_l` could likewise be an accidental alignment
rather than real symmetry. `Tr(M_l^k)` is invariant under any unitary
change of quantization axis (a symmetric function of `M_l`'s eigenvalues),
so it can't be fooled either way. Under exact SO(3) symmetry `M_l` is
proportional to the identity — every eigenvalue equal — so *every* spectral
cumulant of order >=2 vanishes exactly, not just the variance. Per user
direction, implemented this now (pure offline post-processing of data
already on disk, no new C++, no reruns) rather than just writing it up.

**Implementation**, all in `scripts/symmetry_test.py`:
- `matrix_trace_powers(mat, k_max=4)`: `[Tr(M), Tr(M^2), Tr(M^3), Tr(M^4)]`
  via repeated `mat @ mat` (dim is at most a few dozen — trivial cost, no
  diagonalization needed).
- `moments_to_cumulants(p1,p2,p3,p4)`: standard raw-moment-to-cumulant
  formulas (`kappa1=p1`, `kappa2=mu2`, `kappa3=mu3`,
  `kappa4=mu4-3*mu2^2`).
- `reduced_spectral_cumulants(mat)`: combines the above, normalizes by
  `kappa1^k`, on `p_k = Tr(M^k)/dim` (equal-weight moments over the `dim`
  eigenvalues).
- Wired into `analyze()` alongside the existing `R_l`/chi2 computation:
  full-data point estimate plus one value per leave-one-block-out jackknife
  sample (same loop already computing `r_jk`/`d_jk`), fed through the
  existing `jackknife_error()` helper — no new statistical machinery.
  `kappa3_r` requires `dim>=3`, `kappa4_r` requires `dim>=4` (fewer
  eigenvalues can't support a spectral moment distinct from lower ones);
  reported as `--` otherwise, same pattern as the existing chi2 skip.
- Table output extended with `kappa2_r`, `kappa3_r`, `kappa4_r` and their
  jackknife errors; docstring rewritten to explain the
  basis-independent-vs-dependent distinction between the two families.
  `R_l`/chi2 are kept, not replaced — they test a genuinely different
  question (is *this specific* axis anomalous).

**Smoke-tested**: `q=5, n_refine=2, l_max=3, exact_sinh`, 400 traj /
n_skip=2 = 200 samples/blocks (`jack_block_size=1`). Ran end to end with no
errors: `l=0` (`dim=1`, trivial) reports all cumulants as exactly zero
without crashing; `l=1` (`dim=3`) populates `kappa2_r`/`kappa3_r` but
correctly reports `kappa4_r` as `--` (below the `dim>=4` floor); `l=2,3`
populate all three, finite, non-`nan`, same order of magnitude as `R_l`
(e.g. `l=3`: `R_l=0.0495(19)` vs `kappa2_r=0.072(23)`) — physically
sane, both families picking up the same low-statistics/coarse-mesh
symmetry breaking through different, basis-independent-vs-dependent
lenses. Scratch run output deleted after (mechanism check only, per this
campaign's convention).

**PLAN.md updated**: Phase 1 step 6 (gate/success criteria) now names the
reduced spectral cumulants as the primary, basis-independent
symmetry-breaking measure, with `R_l`/chi2 demoted to a secondary
axis-alignment diagnostic.

**Not done**: this has not been run at production statistics or across the
mesh-resolution ladder — same open items as before (aliasing check for
`l_max`, block-merge scan for `jack_block_size`, K_c/coupling framing
already resolved, energy operator and `C_l(l,Delta)` derivation still
out of scope for this phase).

## 2026-08-22 (cont): quick informal refinement diagnostic — no clean
symmetry-restoration trend visible yet at this statistics level

Per user request ("we should see the restoration of symmetry in the large
lattice limit"), ran a login-node diagnostic (not SGE, not the Phase 1
production statistics gate) across part of the mesh ladder:
`q=5, n_refine in {2,4,8,16}` (42/162/642/2562 sites), `exact_sinh`,
`l_max=3` fixed (well below the aliasing threshold at all four
resolutions), `n_therm=1000, n_traj=4000, n_skip=2` (2000 samples/blocks,
`jack_block_size=1`), seed 9001 throughout. Output kept at
`campaign_runs/diagnostic_restoration_2026-08-22/` (not deleted — this is
a real, if low-statistics, data point, not a mechanism-only smoke test).
Wall time trivial even at the largest size (~5s for 2562 sites, 4000 traj).

**Result (kappa2_r, the most robust of the three spectral cumulants at
this statistics level)**:

| l | n_refine=2 | n_refine=4 | n_refine=8 | n_refine=16 |
|---|---|---|---|---|
| 1 | 1.6(1.7)e-3 | 1.9(2.1)e-3 | 4.5(2.9)e-3 | 2.9(2.5)e-3 |
| 2 | 3.3(1.8)e-3 | 3.5(2.0)e-3 | 5.9(2.6)e-3 | 4.0(2.0)e-3 |
| 3 | 9.9(2.6)e-3 | 5.2(1.9)e-3 | 3.7(1.7)e-3 | 4.7(1.8)e-3 |

`l=1,2`: statistically flat across the whole ladder, no significant
refinement trend — all four points overlap within ~1 jackknife sigma.
`l=3`: a real-looking drop from `n_refine=2` to `4` (~1.4 sigma, marginal),
then a plateau around `4e-3` for `n_refine=4,8,16` — could be genuine
residual symmetry breaking surviving refinement, or could be a
statistics-limited plateau; cannot distinguish with 2000 samples and no
block-size tuning.

**Honest conclusion**: this diagnostic does **not** demonstrate symmetry
restoration in the large-lattice limit — it's consistent with either a
slowly-shrinking or a flat/statistics-limited signal. Seeing the actual
large-lattice limit needs the real PLAN.md Phase 1 production statistics
gate (`n_blocks >= 10*(2l+1)`, sized after the block-merge/tau_int
determination), which this quick pass deliberately skipped, plus larger
`n_refine` (32+) to extend the ladder further before concluding anything
about a continuum trend. Not treating this as a negative result either —
just inconclusive at this precision.

**Not done**: no aliasing check, no jack_block_size tuning, no
`n_refine=32` point, no production-statistics run. This was a quick,
deliberately cheap first look, not the Phase 1 production pass.

## 2026-08-22 (cont): package made fully self-contained — reverses the
"do not fork `../IsingS2/`" convention

Per explicit, repeated user direction ("ALL files needed for the project
should be in the project directory", "also the simulator and grp files",
"you did not move all necessary scripts... everything from the code needed
to make the discretizations and the simulator c script and the analysis"),
vendored the complete build/run dependency closure into this package
rather than continuing to read `../IsingS2/` at build or run time:

- `grp/` (61M after dropping `__pycache__`) — icosahedral group-theory
  lattice/elem/site_g data, `rsync`'d from `../IsingS2/grp`.
- `include/` (14M) — `ising.h`, `S2.h`, `lattice.h`, `rng.h`, `grp_o3.h`,
  `statistics.h`, `timer.h`, `util.h`, plus vendored header-only `Eigen`/
  `unsupported` (the same ones `../IsingS2/include/` already vendored from
  upstream Eigen — one hop further removed now, not a new vendoring
  decision).
- `src/ising_s2_crit.cc` — the measurement driver, including the
  already-applied `QfeRng::ReadRng()` fix (see 2026-08-21 entry).
- `cluster/build.sh` rewritten to compile entirely from `$ROOT/include` +
  `$ROOT/src` + `-DGRP_DIR="$ROOT/grp"`, no more `ENGINE_ROOT=../IsingS2`.
  Analysis (`scripts/symmetry_test.py`, new `scripts/plot_restoration_ladder.py`)
  and the plotting venv (`.venv_plot/`) were already local, unaffected.

**Verified**: rebuilt clean, ran `q=5 n_refine=2` at trivial statistics
before and after the switch (fixed seed) — bit-identical `action`/`mag`/
etc. output, confirming the vendored `grp/`+headers reproduce the same
physics as reading `../IsingS2/` directly.

**Consequence, recorded in `directive.md`**: this reverses the original
anti-fork convention, so the vendored copy now needs manual resync — any
future `../IsingS2/`-side fix relevant to this campaign must be re-copied
into `include/`/`src/`/`grp/` here by hand and logged here, since there is
no longer a build-time link forcing consistency. `CLAUDE.md` and
`directive.md` both updated to describe the new self-contained layout and
this obligation.

## 2026-08-22 (cont): higher-statistics login-node push across the
refinement ladder — same plateau, now resolved at 5x the statistics

Per user request to "push" the informal restoration diagnostic further,
reran the `n_refine in {2,4,8,16,32}` ladder (42/162/642/2562/10242 sites)
on the login node (not SGE) with `n_therm=2000, n_traj=20000, n_skip=2`
(10000 samples/blocks, 5x the prior session's 2000), `q=5, exact_sinh,
l_max=3`, seed 4242. New reusable driver: `scripts/run_ladder.sh`. Wall
time: 2m7s for the whole ladder (dominated by `n_refine=32`, 81s alone) —
still trivial for a login-node run. Output at
`campaign_runs/push_2026-08-22/`, plot at
`campaign_runs/push_2026-08-22/plots/symmetry_signals_vs_invL.png`
(`scripts/plot_restoration_ladder.py`, now a saved script rather than an
ad hoc command).

**Result (`kappa2_r`, mean(err))**:

| l | n_refine=2 | n_refine=4 | n_refine=8 | n_refine=16 | n_refine=32 |
|---|---|---|---|---|---|
| 1 | 7.09(5.22)e-4 | 7.75(5.49)e-4 | 4.68(4.44)e-4 | 4.69(4.38)e-4 | 3.45(3.80)e-4 |
| 2 | 5.29(3.25)e-4 | 9.05(4.43)e-4 | 4.77(3.19)e-4 | 7.10(3.91)e-4 | 9.42(4.52)e-4 |
| 3 | 5.50(0.87)e-3 | 2.16(0.55)e-3 | 7.77(3.32)e-4 | 1.11(0.39)e-3 | 1.18(0.40)e-3 |

Same qualitative picture as the 2000-sample pass, now much better
resolved: `l=1,2` flat within errors across the whole ladder (no
refinement trend, consistent with a pure statistical floor around a few
`e-4`). `l=3` drops sharply from `n_refine=2` to `8` (5.5e-3 -> 7.8e-4,
now several-sigma, not the earlier 1.4-sigma-marginal call) then
**plateaus** at `n_refine=8,16,32` around `(0.8-1.2)e-3` rather than
continuing to shrink toward zero — with 10000-sample errors of ~4e-4, this
plateau looks like it may be a real residual (systematic, not
statistical) effect at `l=3`, not a statistics-limited floor.

**Honest conclusion, updated**: still not a demonstration of full
large-lattice symmetry restoration — `l=3`'s plateau above the `l=1,2`
noise floor argues *against* pure statistics-limited residual breaking and
*for* some as-yet-unidentified systematic (candidates, not investigated:
the `exact_sinh` per-link two-face-averaging judgment call flagged in
`CLAUDE.md`/PLAN.md as not a derived result; `l_max=3` aliasing not yet
checked per-resolution; icosahedral mesh's exact residual symmetry group
imprinting preferentially on `l=3` sensitivity). Still not the Phase 1
production-statistics gate (`n_blocks >= 10*(2*l_max+1)`), and no aliasing
check was done here either — but the `l=3` plateau is now a concrete,
reasonably-resolved anomaly worth chasing before declaring Phase 1's
symmetry test passed, rather than the "can't tell yet" of the earlier
pass.

**Not done**: aliasing check, `jack_block_size` tuning, `n_refine=64+`,
diagnosing why `l=3` specifically plateaus above `l=1,2`, production
statistics run.

## 2026-08-23 — `push4_2026-08-22_lmax12`: found the aliasing threshold for
`l_max=12` directly, and hit/fixed a `jack_block_size=1` file-size problem

**Run** (not logged when it was launched, reconstructed from output on
disk): `n_refine in {2,4,8,16,32}`, `q=5, exact_sinh, l_max=12`,
`n_traj=100000, n_therm=?, n_skip=2` (50000 raw samples/point),
`jack_block_size=1` (the tool default, not yet overridden at this point).

**Problem hit**: with `jack_block_size=1` and `l_max=12` (`n_pairs=1547`
upper-triangle entries per sample), each `*_ylm_2pt_full_jackblocks_*.dat`
file came out ~4.8GB (one per `n_refine`, still on disk — not deleted, see
below) and pure-Python parsing in `scripts/symmetry_test.py`/
`plot_restoration_ladder.py` took >2 minutes and didn't reliably finish,
blocking analysis. **Fix**: rebinned each file post-hoc down to 100
jackknife blocks (native per-block sums are additive, so this is a legit
merge, not a re-run) into `rebinned_100blk/q5k*_ylm_2pt_full_jackblocks_rebinned100.dat`
(9.6MB each, ~500x smaller, parses in ~1s). **Saved as a standing rule**:
always pass `--jack_block_size 100` (not the default 1) to `ising_s2_crit`
for this campaign going forward — see
`~/.claude/projects/-projectnb-qfe-misra/memory/feedback_ising_s2_jack_block_size.md`.
This was applied retroactively to the 2026-08-23 production array below
(caught before any task ran, no wasted compute).

**Aliasing finding (the useful result once parseable)**: ran
`scripts/symmetry_test.py`'s `analyze()` on the rebinned files for all
`l=1..12` (not just the `l<=3` the existing plot script hardcodes) across
the full `n_refine in {2,4,8,16,32}` ladder. Result —
`campaign_runs/push4_2026-08-22_lmax12/plots/symmetry_signals_all_l_vs_invL.png`:

- `n_refine=2` (42 sites), `l>=5`: `kappa2_r` is **O(0.1) to O(20)** —
  wildly unphysical (a genuinely SO(3)-symmetric field has `kappa2_r->0`;
  O(1) means the eigenvalue spread of `M_l` is comparable to its mean).
  This is the aliasing signature PLAN.md's Phase 1 step 3 anticipated:
  `l_max=12` needs 91 spherical harmonics resolved and a 42-site mesh
  cannot support that.
- `n_refine=4` (162 sites), `l>=6`: still order `1e-2` to `1e-1` — better
  but still clearly aliased, not just noisy.
- `n_refine=8` (642 sites) and up: **all** `l=1..12` collapse onto a common
  `~1e-4` to `~6e-4` band with no visible `l`-dependent blowup — this is
  the first resolution where `l_max=12` looks safe.
- `l=1..4` are clean (flat, no blowup) at every resolution tested,
  including `n_refine=2`.

**Conclusion — first real answer to PLAN.md Phase 1 step 3's aliasing
check**: for `l_max=12`, the safe threshold is **`n_refine>=8`
(642 sites)**; `n_refine=2,4` are aliasing-contaminated for `l>=5`ish and
should not be used at this `l_max`. This directly informs (but doesn't
replace — only one seed's worth of statistics, no jack_block_size-tuned
error bars) PLAN.md's "single shared `l_max` = smallest safe value across
the whole ladder" rule.

**Consequence for the same-day production array** (`production_2026-08-23`,
`l_max=8` per user override, ladder `n_refine in {4,8,16,32,64,128}`):
`l_max=8` needs fewer harmonics (45 vs 91) than this diagnostic's
`l_max=12`, so the aliasing threshold should be lower, but **has not been
separately checked** — the `n_refine=4` point in that production run is a
plausible aliasing risk by extrapolation from this table (162 sites, same
mesh as this diagnostic's contaminated `l_max=12` case, just at a lower,
less demanding `l_max=8`). Flagged, not yet resolved — check
`production_2026-08-23`'s `n_refine=4` results for the same order-of-
magnitude blowup signature once it completes before trusting them.

**Not done**: `n_refine=64,128` were not part of this diagnostic (only
{2,4,8,16,32}); no jack_block_size scan to confirm 100 blocks is
sufficient at `l_max=12`; the raw (un-rebinned) 4.8GB files per `n_refine`
are still on disk under `campaign_runs/push4_2026-08-22_lmax12/q5k*/` —
harmless but worth deleting once the rebinned versions are confirmed
sufficient for all needed analysis.

## 2026-08-24 — disk cleanup, and `push5_2026-08-24_lmax12` submitted:
l_max=8->12 production, 1M traj/shard, 32 shards, ladder drops n_refine=4

**Cleanup** (per user direction): deleted the five raw (un-rebinned)
`*_ylm_2pt_full_jackblocks_*.dat` files under
`campaign_runs/push4_2026-08-22_lmax12/q5k*/` (~24GB total) — the
`rebinned_100blk/` versions already produced the aliasing-threshold plot
and are kept. Also deleted two stray core dumps at the repo root
(`core.1434606`, `core.1624295`, ~19MB combined, gitignored, unrelated to
any current analysis).

**push5 design** (`cluster/sge/submit_push5.sh`, new — mirrors
`submit_production.sh`'s shape, `production_shard_task.sh` reused
unchanged since it's already fully parameterized via env vars):

- `l_max=12` (was 8 in `production_2026-08-23`). Per user request to push
  `l_max` up now that the queue is open.
- **Ladder drops `n_refine=4`**: `push4_2026-08-22_lmax12`'s aliasing
  diagnostic (previous entry) found `l_max=12` is contaminated below
  `n_refine=8` (kappa2_r blows up at `n_refine=2,4` for `l>=5`) — no point
  spending compute on a ladder point already known to be unusable at this
  `l_max`. Ladder = `{8,16,32,64,128}`.
- `n_traj=1,000,000/shard` (was 20000), `n_shards=32` (was 16) — per user
  request for high statistics and maximum array parallelism.
- `jack_block_size=4000` (was 100). Native block count on disk =
  `n_traj/jack_block_size`; at the old value of 100 this would give 10000
  blocks/file at 1M traj — the same multi-GB-file problem `push4` hit at
  `jack_block_size=1`, just shifted. 4000 keeps native blocks at exactly
  250, which also happens to exactly match PLAN.md Phase 1's statistics
  gate `n_blocks >= 10*(2*l_max+1) = 250` at `l_max=12` — chosen
  deliberately to land on that gate rather than arbitrarily.
- `h_rt=48:00:00` (was 12:00:00) — calibrated on the login node
  (`n_refine=128, l_max=12, n_therm=100, n_traj=2000, n_skip=2` = 4100
  sweeps, wall time 222s -> ~0.054 s/sweep). This is essentially the same
  per-sweep cost as `l_max=8`'s production run (~0.052-0.057 s/sweep
  measured via `qacct` on `production_2026-08-23`'s `n_refine=128` tasks)
  — the Monte Carlo update cost (`O(n_sites)`, 163842 sites at
  `n_refine=128`) dominates over the harmonic-projection cost even at
  `l_max=12`'s 169 harmonics, so raising `l_max` alone didn't change the
  per-sweep budget. Extrapolated worst case: `n_refine=128` shard =
  `(2000 + 1e6*2)` sweeps `* 0.054 s/sweep` ≈ **30 hours**; `h_rt=48:00:00`
  leaves margin above that estimate. Smaller ladder points finish in
  minutes and are not the binding constraint. `seed_base=500000`, disjoint
  from `production_2026-08-23`'s `90000-90095` range.

**Submitted**: `qsub` (detached, unpinned queue, `-P qfe`) — job-array
`7301188`, tasks `1-160` (5 ladder points x 32 shards), queued as of
2026-08-24 18:25. Output under `campaign_runs/push5_2026-08-24_lmax12/`,
logs under `campaign_runs/push5_2026-08-24_lmax12_logs/`, manifest at
`campaign_runs/push5_2026-08-24_lmax12/manifest.json`. Track with `qstat
-u misra | grep s2prec_push5` / `qacct -j 7301188` once tasks complete.

**Not done**: no analysis of push5 output yet (nothing has finished);
`push4`'s remaining scratch dirs (`smoke_test_2026-08-22`,
`sanity_run_2026-08-22`, small — 74K/31M) and the still-untracked git
state were left alone this session, out of scope for what was approved.
Collaborator-handoff git commit explicitly deferred, per user direction
2026-08-24 — do next session or when requested.

## 2026-08-24 (cont): brought `CLAUDE.md`/`PLAN.md`/`directive.md` up to
date with actual execution state — all three had drifted since 2026-08-22

Per user direction to update all the docs, not just this journal. Found
and fixed:

- `CLAUDE.md`: "What this is" still said `src/`/`scripts/`/`configs/` were
  empty scaffolding — stale since the 2026-08-22 self-containment
  vendoring; now says what's actually there and names both production
  pushes. "Cluster conventions" said `cluster/sge/` was empty — it holds
  `production_shard_task.sh` (shared per-task driver, fully parameterized),
  `submit_production.sh`/`submit_push5.sh` (per-push submit scripts, one
  per push rather than edited in place so each push's config stays
  reconstructible), and `calib_test.qsub` (scaffolded, superseded in
  practice by ad hoc login-node timing runs). Added the
  `jack_block_size`-must-be-resized-per-push rule explicitly here too
  (previously only implicit in the memory feedback note).
- `PLAN.md` Phase 1: was still written as "designed 2026-08-22, not yet
  executed" despite two production pushes having actually run. Added a
  **Status** line under each of the 7 steps recording what actually
  happened vs. what was planned — notably: `l_max=8` in
  `production_2026-08-23` was a user override never aliasing-checked (only
  `l_max=12` got the check, via `push4`, informing `push5`'s dropped
  `n_refine=4`); the step-4 tau_int block-merge scan was never built,
  `jack_block_size` was picked by a size heuristic instead; `duality` (step
  1) has only ever been smoke-tested, never run as the planned parallel
  production comparison; and step 6's gate has **never been evaluated**
  against either production dataset — both pushes ran to completion (or
  are running) without anyone checking `kappa2_r/3_r/4_r` against the
  actual pass/fail criterion. Added an explicit status summary at the end
  making this last point unambiguous, since it's the most consequential
  gap: statistics is not the same as an answer.
- `directive.md`: added the `jack_block_size`-resizing rule as a standing
  directive (was previously only in `journal.md`/memory, not in the file
  whose whole job is cross-task standing rules).

**Not done**: did not implement the still-missing step-4 tau_int scan or
run step 6's gate against `production_2026-08-23` — this pass was
documentation-only, no new analysis.

## 2026-08-24 (cont): ran PLAN.md step-6 gate against both production
pushes — first real answer, and it's a fail for l>=3

Added `analyze_blocks`/`read_jackblocks_multi`/`analyze_multi` to
`scripts/symmetry_test.py` (backward-compatible refactor of `analyze()`,
which now just wraps a single-file read) so that a ladder point's 16 or 32
independent per-shard jackblocks files can be pooled into one combined
jackknife ensemble by simple block concatenation (each shard's blocks are
already independent, so concatenation is itself a valid larger
leave-one-block-out set). `main()` now accepts multiple paths. Ran this
across every completed production ladder point:

- `production_2026-08-23` (l_max=8): all 6 points, `n_refine in
  {4,8,16,32,64,128}`, 16 shards each pooled -> n_blocks=1600,
  total_samples=160000 per point (uniform stats across the ladder).
- `push5_2026-08-24_lmax12` + `push6_2026-08-24_lmax12_1Mpooled` (l_max=12):
  `n_refine in {8,16}` from push5 (32 shards, n_blocks=4000,
  total_samples=16,000,000/point) and `{32,64,128}` from push6 (32 shards,
  n_blocks=416, total_samples=500,000/point) — **not uniform stats across
  this ladder** (push5's points have 32x more samples than push6's), a
  side effect of push6 being a recovery run for push5's OOM-killed points;
  flag this when comparing l_max=12 points to each other, though it
  doesn't affect the l_max=8 comparison below, which is the cleaner one.

All output under `campaign_runs/analysis_lmax12_2026-08-24/` (name is
stale, holds both l_max=8 and l_max=12 results — `l8_q5k*` prefix
distinguishes them).

**Gate result (PLAN.md step 6: kappa2_r/3_r/4_r should not grow, ideally
shrink toward zero, as n_refine increases, at each fixed l)**:

Using the l_max=8 ladder (uniform stats, cleanest comparison), `kappa2_r`
at l=3 across `n_refine in {4,8,16,32,64,128}`:
`1.61e-3, 2.17e-4, 4.49e-5, 6.87e-5, 9.36e-5, 1.05e-4` (jackknife errors
~1-3e-5 from `n_refine=16` on). This drops sharply `n_refine=4->16`
(consistent with the known aliasing artifact at coarse mesh, `push4`'s
prior finding) but then **rises monotonically `n_refine=16->128`**, more
than doubling and every step multiple-sigma resolved from zero
(`kappa2_r/err` ~2-5 sigma per point, and the trend itself, not just each
point, is the signal). This is the **opposite** of the gate's "shrink
toward zero" — l=3 spherical-symmetry breaking is not going away as the
mesh refines, at least not up to `n_refine=128` (163842 sites).

`l=8` at l_max=8 (the highest harmonic tested there): `5.47e-3` at
`n_refine=4` (aliased, as expected), then `1.67e-4, 1.42e-4, 1.50e-4,
1.24e-4, 1.38e-4` for `n_refine=8..128` — a flat plateau around `1.2-1.7e-4`,
values overlapping within ~1-2 sigma of each other rather than trending
down. Same qualitative story as l=3: no visible shrinkage once past the
aliased small-mesh regime, just noise around a nonzero plateau.

`l=1,2` stay near a `~1e-5`-`7e-5` noise floor at every resolution, mostly
consistent with zero within jackknife error and with no resolvable trend —
this part of the gate **passes**.

The l_max=12 dataset (despite the mismatched stats across its two source
pushes) shows the same qualitative shape: l=1,2 flat/consistent with zero;
l=3,8,12 sit at a `kappa2_r` plateau of order `1e-5`-`1e-4`,
multiple-sigma resolved from zero at every resolution from `n_refine=8`
through `128`, again not shrinking.

**Conclusion — PLAN.md step 6 gate is not passed as written**: spherical
symmetry looks good (breaking consistent with the statistical noise floor)
for `l=1,2`, but for `l>=3` there is a real, resolvable, non-vanishing
`kappa2_r` plateau that persists out to the largest mesh tested
(`n_refine=128`, 163842 sites) at both `l_max=8` and `l_max=12`. Either (a)
the continuum limit needs meshes finer than anything run so far and the
plateau is a slow crossover this ladder hasn't reached the tail of, or (b)
there is a genuine systematic — candidates already flagged in the
2026-08-22 entry are still live and unreduced: the `exact_sinh` per-link
two-face-averaging judgment call, the icosahedral mesh's residual
point-group symmetry (which naturally couples most strongly to low-l
harmonics like l=3, not to arbitrary l), or an operator-normalization
issue in the `l>=3` harmonics specifically. Not resolved this session —
this is now the central open question blocking Phase 1 closure, not just
"statistics aren't in yet."

**Not done**: `U4` stabilization (the other half of step 6's gate) was not
checked this session — only the `kappa` cumulants. `kappa3_r`/`kappa4_r`
were computed but not separately analyzed here (same `*_symmetry.txt`
files have them, worth a look). No plot made (all analysis via the text
tables in `campaign_runs/analysis_lmax12_2026-08-24/`). The icosahedral
point-group-coupling-to-l=3 hypothesis is not tested against data (e.g. by
comparing `q=5` against a differently-generated mesh family) — flagged as
a candidate, not investigated.

## 2026-08-24 (cont): built and debugged `--mesh_mode`, three candidate
fixes for the l>=3 spherical-symmetry-breaking plateau found above

Motivated directly by the previous entry's finding: `kappa2_r` for l>=3
plateaus at a nonzero, multi-sigma-resolved value instead of shrinking with
mesh refinement, on the naive mesh (plain flat-subdivide-then-radially-
project). Candidate systematic: the naive mesh's triangle-face-area (and
dual/Voronoi-area) non-uniformity (std/mean ~0.12-0.13, roughly resolution-
independent — confirmed empirically this session, see table below) doesn't
shrink with refinement either, so it's a live suspect for something that
couples to low-but->0 harmonics without vanishing in the continuum limit.
Three mesh-construction alternatives were added to `--mesh_mode`
(`src/ising_s2_crit.cc`, `include/S2.h`'s `QfeLatticeS2`):

**`equal_face_area`** (iterative relaxation, `QfeLatticeS2::EqualizeFaceAreas`):
nudges each site toward/away from adjacent-face centroids based on that
face's area imbalance, re-projecting onto the sphere each iteration; a
fixed, purely local, per-site update rule identical for every site, so it
provably preserves the starting mesh's exact icosahedral symmetry (no
group/orbit bookkeeping needs to change). Smoke-tested across
`n_refine in {8,16,32,64}` at `l_max=2`, default 300 iters/step=0.3:
std/mean of face area goes `0.130->0.023`, `0.132->0.014` (n_refine 8, 16 —
real improvement), but `0.133->0.244`, `0.133->0.235` (n_refine 32, 64 —
**gets worse**, likely the fixed 300-iter/step-0.3 schedule not converging
at these larger site counts; not fixed this session, flagged as an open
item below).

**`equal_dual_area`** (iterative relaxation,
`QfeLatticeS2::EqualizeDualAreas`, analogous but targets each site's own
barycentric-proxy dual (Voronoi) area instead of face area): implemented
and smoke-tested, but **diverges** at the same default schedule — post-
relax std/mean (0.40) came out *worse* than pre-relax (0.12) at
`n_refine=4`. **Removed from `--mesh_mode` 2026-08-24 per user direction**:
equalizing face areas already pulls dual areas toward uniform too on this
icosahedrally-symmetric mesh, so a dedicated dual-area-targeting mode
wasn't judged worth separately debugging right now. `EqualizeDualAreas`
itself is left in `include/S2.h` (unused, not wired into the driver) in
case it's revisited.

**`equal_area_analytic`** (per-point closed-ish form, new
`QfeLatticeS2::EqualAreaGnomonicPoint`, generalizing Snyder's 1992
equal-area polyhedral projection): applied at mesh-construction time
(3rd constructor arg `equal_area_analytic`, previously a declared-but-
never-implemented stub inherited from `../IsingS2/` — confirmed via grep
that `EqualAreaGnomonicPoint` had no implementation anywhere, not even in
the origin repo), not as a post-construction relaxation, so it has no
iteration-count/convergence question the other two modes have. Two real
bugs were found and fixed before this worked; both are worth recording in
full since they're easy to reintroduce if this code is touched again:

1. **Edge/vertex mismatch across faces** (found via mesh-construction
   assertion failure, `s_next <= n_sites`, at `n_refine=4`): the first
   version treated `face_r[0]` as a distinguished apex, placing the
   opposite edge (`x+y=k`) by area-ratio bisection relative to that apex,
   and the two near edges (`x=0`, `y=0`) by a different rule (radial
   equal-area fraction from the apex). Two faces sharing an edge generally
   have *different* apexes for that shared edge, so they computed
   different points for the same site — breaks the `coord_map` dedup that
   stitches faces into one mesh. **Fix**: only apply the analytic placement
   to strictly interior grid points (`0<x, 0<y, x+y<k`); boundary points
   keep the naive linear-interpolation-then-normalize placement, which
   depends only on the two shared edge endpoints and position `t` (not on
   the third vertex), so it's inherently face-order-independent. Trade-off:
   only interior-to-face triangles get the equal-area correction — the
   boundary layer keeps naive-mesh distortion, a shrinking fraction
   (`O(1/n_refine)`) of the mesh as resolution increases (consistent with
   the monotonic improvement seen in the table below).

   A second, related issue surfaced once the mesh itself built successfully:
   `QfeLatticeS2::ReadSymmetryData`'s group-element matching/caching
   (`site_g`, `grp/site_g/o3q<q>k<n>.dat`) and its final `UpdateOrbits()`
   call assume every site is an exact group-orbit image of the *naive*
   per-face formula (`CalcOrbitPos`) — hit a `found_g` assertion, and even
   bypassing that, `UpdateOrbits()` would have silently overwritten the
   analytic placement back to the naive one. **Fix**: `ReadSymmetryData`
   takes a new `trust_naive_positions` parameter (constructor passes
   `!equal_area_analytic`); when false, it skips the group-matching/caching
   and `UpdateOrbits()` entirely and returns early, leaving `site_g`
   unpopulated. Confirmed harmless: grepped the whole driver + `S2.h`/
   `lattice.h` — nothing outside `ReadSymmetryData`/`UpdateOrbits` itself
   reads `site_g` or `G`, and `UpdateOrbits` is only reachable via
   `--orbit_path` (unused by default production runs). **Correction (see
   the icosahedral-symmetry check below): the claim made here at the time
   — that the analytic placement is icosahedrally symmetric in principle
   — turned out to be wrong.** It's only symmetric under swapping
   `face_r[1]<->face_r[2]`, not under swapping the apex `face_r[0]` with
   either other vertex, and the base icosahedron's face-vertex labeling
   doesn't consistently preserve "which vertex is the apex" under the
   group action.

2. **Missing Jacobian factor — the real math bug, found via a regression,
   not a crash**: after fix #1, the mesh built and ran, but empirically
   `equal_area_analytic`'s face-area std/mean got *worse* with refinement
   (`0.54->0.72->0.90` for `n_refine in {8,16,32}`), the opposite of the
   goal. Root cause: the per-point construction parametrizes each face by
   `h=(x+y)/k` (a *radius-type* coordinate — linearly spaced grid index
   from the apex `face_r[0]`) and `r` (an azimuth-like coordinate found by
   bisecting for the edge point whose `SphericalTriangleArea` ratio, via
   L'Huilier's theorem, matches `r`). The bug: the radial law used
   `(1-cos(theta_v)) = h*(1-cos(theta_p))` — linear in `h` — which is only
   correct if `h` itself is an area-fraction coordinate. It isn't: `h` is
   linearly spaced in the flat `(x,y)` grid, which makes it a
   *radius*-type coordinate, exactly like ordinary 2D polar coordinates
   `(r,phi)` needing an `r dr dphi` measure (not `dr dphi`) for a
   uniform-in-`r` grid to preserve area under the Cartesian map. Re-deriving
   the Jacobian confirms the fix is a single squared term:
   `(1-cos(theta_v)) = h^2*(1-cos(theta_p))`. Re-tested after the fix and
   confirmed the construction now behaves correctly and monotonically:

   | n_refine | naive std/mean | equal_area_analytic std/mean |
   |---|---|---|
   | 8  | 0.130 | 0.080 |
   | 16 | 0.132 | 0.054 |
   | 32 | 0.133 | 0.036 |
   | 64 | 0.133 | 0.024 |

   Monotonically shrinking as `n_refine` grows, consistent with the
   boundary-layer-only-correction trade-off from fix #1 (`O(1/n_refine)`
   fraction of the mesh still naive) — this is the expected qualitative
   signature of a genuinely-working equal-area correction, unlike the
   pre-fix version's divergent trend.

Full pipeline smoke-tested end to end after both fixes: `--mesh_mode
equal_area_analytic --n_refine 8 --l_max 4`, 20000 traj, ran to completion,
`scripts/symmetry_test.py` processed the output without error (see
`campaign_runs/smoke_test_2026-08-24_eqarea/`, since deleted — this was a
mechanism smoke test only, no production statistics kept).

**Not done**: no production-scale run of `equal_area_analytic` yet (only
smoke-tested at low stats) — the actual question this was all motivated
by, whether it *fixes* the l>=3 `kappa2_r` plateau from the previous
entry, is still open and needs a real production ladder + the same
step-6-gate analysis. `equal_face_area`'s divergence at `n_refine>=32`
(bad iteration schedule) is unfixed. `--jack_block_size` needs the usual
per-push resizing before any real production run (see
`feedback_ising_s2_jack_block_size` memory /
[[feedback_ising_s2_jack_block_size]]). Boundary-layer-only correction
means a direct A/B against `equal_face_area` (which corrects the whole
mesh, when it converges) at matched `n_refine` would be informative but
wasn't run.

## 2026-08-24 (cont): direct empirical check of icosahedral-symmetry
retention — `equal_area_analytic` FAILS it, ~1.5% level, not safe for
production as-is

Per user request to actually verify (not just argue) that the new mesh
modes retain icosahedral symmetry. Wrote a standalone diagnostic
(`src/test_ico_symmetry.cc`, compiled ad hoc, not part of the production
build/committed — deleted the binary after use, source kept in case this
check needs to be repeated) that, for every site `s` and every icosahedral
group element `g` (loaded the same way `ReadSymmetryData` does, from
`grp/elem/o3q5.dat`), applies `G[g]` to `r[s]` and brute-force searches for
the nearest actual site position, reporting the worst-case mismatch. A
mesh that is exactly closed under the group has mismatch at floating-point
round-off; any larger residual is real symmetry breaking.

- **`naive` mesh, n_refine=4**: max mismatch `6.8e-16` — exactly symmetric,
  as expected (confirms the diagnostic itself works and matches the
  existing `ReadSymmetryData` group-matching machinery's assumption for
  this mode).
- **`equal_area_analytic`, n_refine=4**: max mismatch `1.485e-2` — **not**
  symmetric, a real ~1.5% effect, not round-off.
- **`equal_area_analytic`, n_refine=8**: max mismatch `1.495e-2` — same
  size, doesn't shrink with resolution (consistent with a fixed
  fractional-position construction defect, not a floating-point issue).

**Root cause, and why the previous entry's symmetry claim was wrong**:
`EqualAreaGnomonicPoint` treats `face_r[0]` as a distinguished apex. It
*is* symmetric under swapping `face_r[1]<->face_r[2]` (both the area-ratio
bisection and the radial law only depend on the pair through interchangeable
roles), but it is **not** symmetric under swapping the apex `face_r[0]`
with either other vertex — a completely different formula would be used
for "the same" grid position if a different vertex played apex. The
combinatorial mesh (which face gets which of its 3 vertices labeled
`face_r[0]`, inherited from `ReadBaseLattice`'s file order and the
recursive-subdivision loop) is not guaranteed to assign the apex role
consistently under the icosahedral group's action — and empirically, it
doesn't. The naive formula never had this problem because
`xi0*face_r[0]+xi1*face_r[1]+xi2*face_r[2]` is symmetric under all `3!`
permutations of the vertices, so any labeling choice gives the same
answer.

**Consequence**: `equal_area_analytic` is **not safe to use for
production** in its current form — a mesh that isn't closed under the
icosahedral group would contaminate exactly the `kappa2_r`
spherical-symmetry-breaking measurement this mode exists to fix, making
results uninterpretable (can't distinguish "genuine curved-mesh
discretization artifact" from "artifact of this code's own apex-labeling
asymmetry"). The earlier monotonic-improvement-with-refinement table
(std/mean of face area) is still numerically correct as a measurement, but
doesn't by itself certify the construction is fit for purpose here — a
mesh can have very uniform triangle areas and still not respect the point
group. Not fixed this session (open item below) — flagged immediately
rather than left for a future session to rediscover via contaminated
production results.

**Fixed same session, via a different strategy than per-point symmetrization**
(see next entry): rather than making `EqualAreaGnomonicPoint` itself
invariant under all `3!` vertex permutations (would need a real
re-derivation), exploit the mesh's pre-existing orbit/group-element
bookkeeping (`site_orbit`, `site_g`, `G` — the same machinery
`UpdateOrbits()` already uses) to *propagate* one analytic representative
point per orbit to every other site in that orbit as an exact
group-element image. This guarantees closure under the icosahedral group
by construction, regardless of the per-point formula's own symmetry
properties, as long as the propagation only uses `site_g`/`site_orbit`
values that were computed against the (exactly-symmetric) naive
placement — see `QfeLatticeS2::ApplyEqualAreaAnalytic` in `include/S2.h`
and the constructor, which now always builds the naive mesh first, runs
`ReadSymmetryData` unchanged (`site_g`/`site_orbit`/`G` valid), and only
then — if `equal_area_analytic` — overwrites every site's position via
this propagation pass.

Two more bugs surfaced getting this actually correct (both confirmed via
the `test_ico_symmetry` diagnostic — build ad hoc, source kept at
`src/test_ico_symmetry.cc`, not part of the production build, checks every
`(site, group element)` pair for an exact match to some other site):

1. **Orbit weight/vertex-index mismatch**: `orbit_xi[o]` stores
   sorted-descending barycentric fractions, and `CalcOrbitPos` assigns
   `orbit_xi[o](0)` (the *largest*) to `face_r[0]`, `(1)` to `face_r[1]`,
   `(2)` to `face_r[2]`. Recovering integer grid indices `(x,y)` for
   `EqualAreaGnomonicPoint` (which uses the naive-formula convention
   `weight_on_face_r[0]=(k-x-y)/k`, `weight_on_face_r[1]=x/k`,
   `weight_on_face_r[2]=y/k`) requires `x` from `orbit_xi[o](1)` and `y`
   from `orbit_xi[o](2)` — **not** `(0)` and `(1)`, which is what the first
   attempt used. This alone made the mismatch *worse* (`0.015 -> 0.30`),
   a useful confirming regression once found and fixed (back to machine
   precision at `n_refine=4`, but see #2 below for `n_refine>=8`).

2. **Orbits with a nontrivial icosahedral stabilizer**: an orbit whose
   naive representative sits exactly on a symmetry axis/mirror plane has
   fewer than `n_group` (120, full icosahedral group incl. reflections)
   sites sharing it — e.g. at `n_refine=8`, orbit `xi=(0.5,0.25,0.25)`
   (grid point `(x,y)=(2,2)`, on the median from the face apex) has
   multiplicity 60, not 120. Propagating via the *existing* `site_g` coset
   representatives (computed against the naive, stabilizer-respecting
   point) is only valid if the *replacement* analytic point also respects
   that same stabilizer — otherwise the propagated image set doesn't
   close under the full group (confirmed: this exact orbit was the source
   of the residual `1.5e-2` mismatch at `n_refine=8` after fix #1, despite
   `n_refine=4`'s single interior median point `(1,1)` working fine by
   luck of having a trivial stabilizer at that specific resolution).
   **Fix**: `ApplyEqualAreaAnalytic` now computes each orbit's combinatorial
   multiplicity (`orbit_mult`, via `site_orbit`) and falls back to the
   (always-safe) naive `CalcOrbitPos` for any orbit with
   `multiplicity != G.size()`, in addition to the existing boundary
   fallback. Re-verified with `test_ico_symmetry`: exact machine-precision
   closure (`~1e-15`) at `n_refine in {4,8,16}`, up from the `1.5e-2`
   failure.

**Trade-off exposed by fix #2, not yet resolved**: the multiplicity check
is a *conservative* proxy — it flags every symmetry-axis orbit, but I can
show analytically that at least the `x==y` (median-line) case is actually
safe: by the face's own mirror symmetry (equilateral triangle, swapping
`face_r[1]<->face_r[2]` fixes `p(0.5)` = the edge chord midpoint, and the
area-bisection step already proven to land exactly at `t=0.5` there by a
direct reflection argument), so `EqualAreaGnomonicPoint` genuinely does
respect that specific stabilizer and doesn't need the fallback. Because
the multiplicity check can't distinguish "safe" stabilizers (like this
one) from potentially-unsafe ones without checking each case, it currently
reverts the *entire* median line (not just isolated fixed points) back to
naive placement — and empirically this makes area-uniformity **worse
than naive**, and worsening with resolution:

| n_refine | naive std/mean | equal_area_analytic std/mean (post-fix) |
|---|---|---|
| 8  | 0.130 | 0.154 |
| 16 | 0.132 | 0.213 |
| 32 | 0.133 | 0.283 |
| 64 | 0.133 | 0.342 |

So as of this session's end: **icosahedral symmetry is exactly verified
and safe** (the property the user explicitly asked to check and fix), but
the area-uniformity benefit that motivated this mode in the first place is
currently net-negative due to the conservative fallback. Next step, not
done this session: either (a) prove which stabilizer types are actually
safe (the `x==y` reflection case, and possibly others) and narrow the
fallback to only the genuinely-unsafe orbits, recovering the improvement
seen before fix #2 while keeping exact symmetry, or (b) accept the
symmetry-safe/uniformity-regressed trade-off and prioritize
`equal_face_area` instead. `equal_face_area` itself has not been checked
with `test_ico_symmetry` — it's argued (fixed, purely local, per-site
update rule identical for every site) rather than empirically verified,
and should be before being trusted for production either.

## 2026-08-25 — verified the gradient-descent `EqualizeFaceAreas` rewrite
(true per-vertex triangle-area gradient + backtracking line search,
replacing the earlier fixed-unit-step heuristic that diverged at
`n_refine>=32` regardless of step size — code and derivation already
in `include/S2.h`/`QfeLatticeS2::EqualizeFaceAreas` from a prior
session, but never smoke-tested at production resolutions or checked for
icosahedral-symmetry closure; this session did both)

**Divergence at n_refine=32/64 is fixed** — built two ad hoc diagnostics
(not part of `cluster/build.sh`, same compile pattern as
`test_ico_symmetry.cc`):

- `src/test_gd_equal_area.cc` — runs `EqualizeFaceAreas` across
  `n_refine in {8,16,32,64}` and reports face-area std/mean before/after.
  At the old default schedule (300 iters, step=0.3):

  | n_refine | naive std/mean | relaxed (300 iter) |
  |---|---|---|
  | 8  | 0.130 | 0.023 |
  | 16 | 0.132 | 0.012 |
  | 32 | 0.133 | 0.0058 |
  | 64 | 0.133 | 0.030 |

  32 now converges cleanly (vs. the old heuristic's 0.130->0.244
  divergence). 64 looks worse than 32 at a fixed 300-iter budget, but this
  is just under-iteration, not divergence: bumping to 1000 iters gives
  `n_refine=64 -> 0.0030`, better than 32's 300-iter result and still
  improving. Unlike the old heuristic (more iterations didn't help --
  confirmed drift away from the best point reached regardless of step),
  the backtracking line search means more iterations only ever help or
  plateau, never diverge. Runtime is cheap (full 4-point ladder at 1000
  iters: ~9.6s on the login node), so there's no real cost to just running
  more iterations at higher `n_refine` rather than hand-tuning a schedule.
  **Practical takeaway for production**: `--equal_area_iters` should scale
  up with `n_refine` (or just default to a generously large budget like
  1000-2000 — backtracking makes over-iterating harmless, it just stops
  taking steps once converged).

- `src/test_gd_symmetry.cc` — same brute-force `(site, group element)`
  closure check as `test_ico_symmetry.cc`, but applied after
  `EqualizeFaceAreas` rather than to the naive/analytic constructions.
  `O(n_sites^2 * n_group)`, so only run at `n_refine in {4,8}` (642 sites
  already ~3e10 distance evaluations). Confirms the argued-but-unverified
  claim from the 2026-08-24 entry: exact machine-precision closure
  (`max_mismatch ~ 8e-16`, `n_refine=4`; `~9e-16`, `n_refine=8`) —
  `EqualizeFaceAreas`'s fixed, identical, purely-local per-site update rule
  does preserve the naive mesh's icosahedral symmetry in practice, not just
  in principle.

**Net result**: `equal_face_area` is now the strongest of the three
`--mesh_mode` candidates for the l>=3 `kappa2_r` plateau fix — better
area-uniformity than `equal_area_analytic` at matched `n_refine` (e.g.
0.0030 vs 0.024 at `n_refine=64`, with the boundary-layer-only correction
`equal_area_analytic` is stuck with) and, as of this session, no known
divergence or unverified-symmetry caveat left. `equal_dual_area` remains
unfixed/unwired (2026-08-24 decision, unchanged).

**Not done this session**: no production-scale Ising run with
`--mesh_mode equal_face_area` yet — everything above is mesh-geometry-only
(`FlatArea`/closure checks), not the actual `kappa2_r` spherical-symmetry
observable this whole effort is meant to fix. That's the next step: rerun
the aliasing-safe production ladder (see PLAN.md step 3/step 6 gate) with
`--mesh_mode equal_face_area --equal_area_iters <n_refine-scaled budget>`
and check whether the l>=3 plateau actually shrinks now. `--jack_block_size`
needs the usual per-push resizing first (see
[[feedback_ising_s2_jack_block_size]]). Diagnostic binaries
(`test_gd_equal_area`, `test_gd_symmetry`) are ad hoc, not part of
`cluster/build.sh`, same as `test_ico_symmetry` — keep the source but
don't wire into the production build.

## 2026-08-25 (cont): timing/performance of `EqualizeFaceAreas`

Wall-clock cost of the gradient-descent relaxation itself (login node,
single-threaded, `n_iter=1000, step=0.3` — enough to fully converge at
every ladder point per the previous entry), measured with an ad hoc timing
harness (not kept as a source file — trivial `std::chrono` wrap around
`EqualizeFaceAreas`, same construction pattern as `test_gd_equal_area.cc`):

| n_refine | n_faces | time |
|---|---|---|
| 8   | 1,280   | 0.014s |
| 16  | 5,120   | 0.40s |
| 32  | 20,480  | 1.7s |
| 64  | 81,920  | 7.4s |
| 128 | 327,680 | 29.5s |

Scales close to linearly in `n_faces` per iteration, as expected (each
iteration is one O(n_faces) gradient accumulation pass plus a handful of
O(n_sites) backtracking trial evaluations) — roughly 4x runtime per
doubling of `n_refine` since `n_faces` itself grows ~4x. Negligible next to
production Ising run cost (minutes-hours) even at `n_refine=128`, so mesh
relaxation is not a budget concern for the next production push.

Combined with the area-uniformity result from the same session: a ~20-45x
reduction in std/mean at essentially free wall-clock cost, and — because
backtracking makes over-iterating harmless (it just stops taking steps
once converged) — there is no real reason to hand-tune `--equal_area_iters`
per resolution. **Recommendation for the next production push**: default
`--equal_area_iters` to something generously high (1000-2000) across the
whole ladder rather than reusing the old 300-iter default (which
under-converges at `n_refine=64`, see previous entry's table) — not
changed in the driver's compiled-in default this session, just flagged
here so the next push script sets it explicitly.

## 2026-08-25 (cont): simplified `--mesh_mode` to two options, dropped
`--orbit_path` from the driver CLI (per user direction)

Per explicit user direction ("only three options for now: icosahedron,
naive projection, equal areas projection" -- clarified via AskUserQuestion
to mean: no literal third mode, since `naive` already *is* the
icosahedron-based projection; just drop `equal_area_analytic` and rename
`equal_face_area` -> `equal_area`) and a follow-up direction ("do not use
the grp and orbits — needlessly complicated for now"):

- `--mesh_mode` now only accepts `naive` or `equal_area` (was `naive`,
  `equal_face_area`, `equal_area_analytic`). `equal_area` is the
  gradient-descent `EqualizeFaceAreas` mode verified earlier this session
  (machine-precision symmetry closure, best area-uniformity, no known
  divergence).
- `equal_area_iters` default bumped from 300 to 1000 in the driver, per
  this session's earlier timing/convergence finding (300 under-converges
  at `n_refine=64`; 1000 converges the whole ladder at negligible cost).
- Dropped `--orbit_path`/`ReadOrbits` entirely from the driver -- its only
  purpose was overwriting site positions via the group/orbit machinery
  (`site_g`/`G`/`UpdateOrbits`), unused now that `equal_area_analytic` is
  gone.
- **`include/S2.h` is deliberately untouched** (confirmed `git diff
  --stat include/S2.h` is empty) -- per explicit user direction, this was
  a driver-level (`src/ising_s2_crit.cc`) simplification only, not a
  library cleanup. `QfeLatticeS2`'s 3-arg constructor
  (`equal_area_analytic` param), `EqualAreaGnomonicPoint`,
  `ApplyEqualAreaAnalytic`, `ReadSymmetryData`, `UpdateOrbits`,
  `ReadOrbits`/`WriteOrbits`, `G`, `site_g` all still exist in the library,
  just no longer called from this driver (the constructor call now omits
  the 3rd arg, which defaults to `false`). If a future session wants
  `equal_area_analytic` or explicit orbit loading back, the code is still
  there.
- Smoke-tested after rebuild (`bash cluster/build.sh`, clean): `--mesh_mode
  naive` and `--mesh_mode equal_area --equal_area_iters 200` both run to
  completion at `n_refine=4, l_max=2`; `--mesh_mode equal_area_analytic`
  and `--orbit_path` are correctly rejected (unknown mesh_mode / getopt
  unrecognized-option error).

**Caveat inherited, not fixed this session**: `QfeLatticeS2`'s constructor
still unconditionally calls `ReadSymmetryData(q, k)` internally (reads
`grp/elem/o3q*.dat` and `grp/site_g/o3q*k*.dat`, populates `G`/`site_g`)
regardless of driver flags, since that's inside the untouched library
constructor, not gated by the removed CLI options. This has always been
true (naive-mode production runs already paid this cost) and is harmless
(cheap relative to the Ising run itself, and no longer read by anything in
this driver), but means "not using grp and orbits" is true at the driver's
CLI surface, not literally true of every code path the constructor
touches -- flagging in case it matters for a future speed-focused pass.

## 2026-08-25 (cont): mesh-position cache, iters-needed model, and a real
early-stop bug found and fixed

Per user direction: (1) build a persistent cache so `EqualizeFaceAreas`
never has to rerun for a `(q, n_refine, step)` already relaxed, and (2)
develop a model for how many GD iterations a given `n_refine` needs, since
the cache-miss path needs a sane iteration budget rather than either a
fixed guess or an unbounded run.

**Cache**: added `QfeLatticeS2::WritePositions`/`ReadPositions` to
`include/S2.h` (flat binary dump of `r[]`, `n_sites` header as a staleness
check) and a `--mesh_cache_dir` flag to `ising_s2_crit.cc` -- when set and
`--mesh_mode equal_area`, the driver checks
`<mesh_cache_dir>/q<q>k<n_refine>_step<step>.dat` first; on a hit it loads
positions and skips relaxation entirely; on a miss it relaxes as before
and writes the cache entry. Cache is keyed only by `(q, n_refine, step)` --
not scoped to a single campaign run -- so it's meant to be a shared,
growing library under `mesh_cache/` at the repo root, reused across every
future push at a given resolution. **Concurrency note**: the driver itself
does no locking -- if multiple shards for the same `(n_refine, step)` hit
a cold cache concurrently they will race to `fopen(...,"wb")` the same
file. `submit_push7_dual.sh` avoids this by pre-populating the cache
*serially* on the login node before submitting the SGE array; any future
submit script doing many shards per equal_area point must do the same or
risk a corrupted cache file.

**Iters-needed model / early-stop**: `EqualizeFaceAreas` previously only
stopped when backtracking found literally no improving step in 30 tries
(the "hard stop") -- found not to fire within thousands of iterations at
`n_refine>=32` (it keeps finding minuscule improving steps indefinitely).
Added a practical relative-improvement early stop: track `E`'s relative
decrease per iteration, stop once it's been below `rel_tol` (default
`1e-5`) for `patience` (default `20`) consecutive iterations.
`EqualizeFaceAreas` now takes optional `iters_used`/`rel_tol`/`patience`
params (all default-backward-compatible).

**Bug found while validating this against `n_refine=128`**: the naive
version of this early stop fired after exactly 20 iterations at
`n_refine=128`, leaving the mesh essentially unrelaxed
(`relaxed_std_mean=0.133400` vs naive `0.133403` -- no real improvement).
Root cause (confirmed via a `GD_DEBUG=1` env-gated per-iteration trace
added to `EqualizeFaceAreas`): the backtracking step size ramps up 1.2x
per accepted iteration starting from a fixed absolute `step=0.3`
regardless of mesh scale. At `n_refine=128`, `E` starts far smaller
(~8.6e-6, since individual face areas shrink as `n_faces` grows) than at
coarser resolutions, so the *warm-up* phase where the step is still
ramping up produces relative-decrease values that are themselves still
small and, critically, still monotonically *increasing* every iteration
(1.9e-7 -> 6.0e-6 over iters 0-19) -- the naive "below rel_tol" check
can't distinguish this acceleration phase from genuine stalling and
declared victory right as the real progress was about to start. Fixed by
only counting an iteration toward the "stalled" streak when
`rel_decrease` is *not* still increasing relative to the previous
iteration (`rel_decrease <= prev_rel_decrease` guard, in addition to
`rel_decrease < rel_tol`) -- see the in-code comment on
`EqualizeFaceAreas` for the exact logic. This was a real, first-caught
bug, not a theoretical concern: the buggy binary was already used for a
`bin/test_gd_stress` L=256..2048 stress sweep before the fix was found,
which must be treated as invalid/unrun.

**Data (fixed binary), `campaign_runs/gd_stress_2026-08-25/`,
`sweep_pow2_pow2x3_upto128.log`** -- ladder = powers of 2 and powers of
2x3 up to 128, `n_iter` cap 20000, `step=0.3`, login node:

| n_refine | n_faces | naive std/mean | relaxed std/mean | iters used | relax time |
|---|---|---|---|---|---|
| 2   | 80      | 0.0775 | 0.0775  | 0     | 0.000s |
| 3   | 180     | 0.1062 | 0.0574  | 23    | 0.002s |
| 4   | 320     | 0.1181 | 0.0446  | 28    | 0.002s |
| 6   | 720     | 0.1266 | 0.0304  | 44    | 0.007s |
| 8   | 1,280   | 0.1296 | 0.0230  | 51    | 0.006s |
| 12  | 2,880   | 0.1317 | 0.0154  | 69    | 0.021s |
| 16  | 5,120   | 0.1325 | 0.0116  | 107   | 0.045s |
| 24  | 11,520  | 0.1330 | 0.00773 | 221   | 0.206s |
| 32  | 20,480  | 0.1332 | 0.00580 | 388   | 0.654s |
| 48  | 46,080  | 0.1333 | 0.00387 | 852   | 3.32s |
| 64  | 81,920  | 0.1334 | 0.00290 | 1,532 | 11.1s |
| 96  | 184,320 | 0.1334 | 0.00194 | 3,452 | 58.3s |
| 128 | 327,680 | 0.1334 | 0.00145 | 20,000 (cap; did not early-stop) | 604.8s |

Every point 2-96 triggers the (now-correct) early stop comfortably inside
the 20000 cap; n_refine=128 ran the full cap without stopping (still
slowly improving -- the plateau is real but `patience=20`/`rel_tol=1e-5`
isn't tight enough to catch it at this size in reasonable time). Naive
std/mean saturates at ~0.133 independent of resolution as already known;
relaxed std/mean keeps improving smoothly through the whole range checked
so far, no sign of a plateau by n_refine=96. `iters_used` grows
substantially faster than linear in `n_faces` (roughly `n_faces^1.2`
between adjacent points in the 8-96 range) -- no clean single power law
was fit; the 20000 cap (driver default as of this session) is a safety
ceiling sized to comfortably early-stop every point up to 96, not a
literal prediction formula. **The `n_refine in {256,512,1024,2048}` stress
sweep the user separately asked for has not yet been (re-)run with the
fixed binary** -- the pre-fix run (`sweep_large_256_to_2048.log`) is
invalid and was not redone this session; next session should rerun it if
still wanted.

## 2026-08-25 (cont): push7 -- naive vs equal_area head-to-head, submitted

User direction (after briefly considering and then explicitly dropping a
third "unprojected icosahedron" mesh_mode -- flat, un-normalized
subdivided positions -- which would have needed real design work to
reconcile with the Y_lm/quadrature machinery's unit-vector assumption;
not implemented): run a two-way `naive` vs `equal_area` production
comparison across ladder `n_refine in {2,3,4,6,8,12,16,24,32,48,64,96}`
(powers of 2 and powers of 2x3 up to 96), 10,000,000 total pooled
trajectories split evenly across all 24 `(n_refine, mesh_mode)`
combinations (416,000 each: 32 shards x 13,000 traj/shard), "aggressive
parallelism" (32 shards/combo, both mesh_mode arms submitted as separate
detached SGE arrays, 384 tasks each).

**`l_max=3`, not the usual production `l_max=12`** -- this ladder's
smallest point (`n_refine=2`, 42 sites) cannot safely support `l_max=12`
(169 modes vs 42 sites, guaranteed aliased per push4's established
`n_refine>=8` threshold at that l_max). `l_max=3` (16 modes) reuses the
value the 2026-08-22 informal diagnostic already ran successfully down to
`n_refine=2` -- a conservative reused choice, not a freshly derived
per-resolution aliasing check (PLAN.md step 3 remains open). This push
therefore only tests the mesh-mode question at `l=1..3`, not the full
l-range where the naive-mesh plateau was originally seen at `l>=3`
(borderline included, not the higher l where the plateau was clearest).

**New infrastructure**: `cluster/sge/submit_push7_dual.sh` -- writes a
manifest, then *serially* prewarms the equal_area mesh-position cache
(`mesh_cache/`, one `ising_s2_crit` invocation per ladder point with
`n_traj=1` just to trigger+cache the relaxation, discarded scratch output)
before submitting the equal_area array, to avoid 32 concurrent shards
racing on a cold cache file. `production_shard_task.sh` extended to pass
`--mesh_cache_dir` through when set. Seed ranges: naive
`700000-700383`, equal_area `710000-710383` (disjoint from all prior
pushes' ranges).

**Status**: submitted this session, not yet evaluated -- track with
`qstat -u misra | grep s2prec_push7`. `campaign_runs/push7_2026-08-25_dual_naive_eqarea/manifest.json`
has the exact parameters.

## 2026-08-25 (cont): push7 resubmission fix (login-node contention), ladder
extended to n_refine=128, and PLAN.md step 4 (tau_int) closed by argument

**Original push7 submission never actually ran** -- its serial, login-node
mesh-cache prewarm loop (see previous entry) stalled indefinitely at
`n_refine=64`, so neither production `qsub` was ever reached (`qstat`
showed nothing because nothing was submitted; all 768 pre-created shard
directories stayed empty). Root cause confirmed live: this session's shell
is limited to `nproc`=1 against a login-node load average of ~11-25
(`uptime`, 35 concurrent users) -- reproduced the exact stall by rerunning
the same prewarm invocation with `GD_DEBUG=1`, which never advanced past
the initial gradient pass within 120s, versus the `gd_stress_2026-08-25`
benchmark's uncontended ~11.1s for that same `n_refine=64` point. Not an
`EqualizeFaceAreas` correctness issue.

**Fix**: split the prewarm step out of `submit_push7_dual.sh` into a new
`cluster/sge/prewarm_mesh_cache_task.sh`, submitted as its own `qsub -sync
y` SGE array (one task per ladder point -- independent cache files, no
race) instead of a login-node loop. Runs on dedicated compute-node cores,
avoiding the contention that caused the stall. `submit_push7_dual.sh` now
verifies every expected cache file exists after the prewarm array
completes before submitting the two production arrays.

**Ladder extended per user direction**: added `n_refine=128` (13 points
now: 2,3,4,6,8,12,16,24,32,48,64,96,128; was capped at 96) to match
`l_max=8`'s already-established-safe aliasing range and reach the same top
resolution `production_2026-08-23`/`push5`/`push6` used. Note (see
`campaign_runs/gd_stress_2026-08-25/`): at `n_refine=128` the GD relax
step does not actually early-stop within the 20000-iteration cap (still
slowly improving when capped, ~604.8s uncontended) -- the prewarm SGE
array's `h_rt=00:30:00` covers this with margin, but the reported
area-uniformity at that point is a lower bound, not a fully-converged
value.

**Resubmitted** under the same `push7_2026-08-25_dual_naive_eqarea` tag
(safe -- prior attempt wrote no real output). Prewarm array job-ID
`7305708` (13 tasks) started 01:07:19 EDT; both production arrays
(`naive`, `equal_area`, 416 tasks each = 13 points x 32 shards) submit
automatically once it clears.

**PLAN.md step 4 (tau_int-based block-merge scan) -- closed by argument,
not by running the scan, per user direction this session**: user pointed
out `tau_int ~ O(1)` in units of cluster sweeps is the expected behavior
for Wolff cluster updates near criticality specifically *because* that
algorithm is designed to eliminate critical slowing down (unlike local
Metropolis, where `tau_int` diverges with the correlation length). This
campaign's driver always includes `n_wolff=5` cluster updates per
trajectory (`ising_s2_crit.cc`, all pushes so far), so the standing
open-item framing in `PLAN.md` step 4/5 (200-250 native blocks in past
production runs "counts native blocks, not tau_int-verified independent
samples") is over-cautious for this update scheme -- native block count
and effective-independent-sample count should closely track each other
here. This does not retroactively run the formal `--merge-blocks` sweep
`PLAN.md` step 4 originally called for, but does justify treating existing
and push7's jackknife error bars at face value rather than as a
lower-bound-on-error caveat. Push7 pools 4,160 native blocks/combo
(`l_max=8` gate is 170) -- ~16-20x more blocks than the stats level that
already resolved the `l>=3` `kappa2_r` plateau at multi-sigma significance
in `production_2026-08-23`, so statistics are not the limiting factor for
reading push7's naive-vs-equal_area comparison cleanly. The remaining
caveat is aliasing, not statistics: `n_refine=2,3` sit below the
empirically-established (but only informally confirmed) `l_max=8`
aliasing-safe floor of `n_refine>=4` -- check `R_l`/`kappa2_r` at `l=7,8`
for those two points specifically before including them in the mesh-mode
comparison.

## 2026-08-25 (cont): push7 actually launched; new `gd_stress_1024` optimizer
sweep submitted (both in progress, checking back in a few hours)

**push7 second resubmission -- fixed a real bug in the prewarm-verification
check added earlier this session**: the first fixed submission's cache-file
existence check compared against `q5k${k}_step${EQUAL_AREA_STEP}.dat` using
the raw shell variable (`0.3`), but `ising_s2_crit.cc`'s `WritePositions`
path formats the step with `%.3f` (`0.300`) -- a filename mismatch, not a
missing cache, produced a false `ERROR: prewarm did not produce cache for
n_refine=2` after the (successful) first prewarm array run. Fixed by
padding via `printf '%.3f'` before comparing (`submit_push7_dual.sh`).
Confirmed the underlying prewarm array itself was fine on the first
attempt -- all 13 tasks exited 0, all 13 `mesh_cache/q5k*_step0.300.dat`
files present -- so the resubmission's prewarm stage re-verified cache
hits (no recomputation) and should move straight to submitting the two
production arrays. **Status as of 01:25 EDT**: prewarm re-verification
array (job `7305749`, 13 tasks) running; production arrays not yet
confirmed submitted (script still inside its `qsub -sync y` wait) --
check `campaign_runs/push7_2026-08-25_dual_naive_eqarea_logs/` and `qstat
-u misra | grep s2prec_push7` for `_naive`/`_eqarea` array job IDs next
session. Ladder is the corrected 13-point `{2,3,4,6,8,12,16,24,32,48,64,
96,128}`, `l_max=8`, per the previous entry.

**New: `gd_stress_1024_2026-08-25`, EqualizeFaceAreas performance sweep
extended to `n_refine=1024`** (per user direction: "schedule the runs
necessary to see the performance of the optimizer for L powers of 2 and
powers of 2 times 3 up to 1024. give this lots of wall time and lots of
power"). Supersedes the never-validly-rerun `sweep_large_256_to_2048.log`
(pre-early-stop-fix binary, invalid). New infra:
`cluster/sge/gd_stress_task.sh` (one SGE array task per ladder point,
independent -- no shared state, no serialization needed, unlike the
mesh-cache prewarm) and `cluster/sge/submit_gd_stress_1024.sh`. Rebuilt
`bin/test_gd_stress` fresh from current `include/S2.h` first (confirmed
stale-check: previous binary's mtime exactly matched `S2.h`'s already, so
it was current, but rebuilt anyway to be certain before a multi-hour run).

- Ladder: `2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024` (19
  points, powers of 2 and powers of 2x3, extending
  `gd_stress_2026-08-25/sweep_pow2_pow2x3_upto128.log`'s pattern).
- `n_iter` cap 20000, `step=0.3` -- same as the validated `<=128` sweep.
- `h_rt=72:00:00` and `-pe omp 16` (16 cores reserved purely for the SCC
  per-core memory allocation -- `test_gd_stress` itself is
  single-threaded) per explicit user direction for generous
  wall-time/resources. Sizing rationale: per-iteration cost empirically
  scales ~linearly with `n_faces=20*n_refine^2`, and every point from
  `n_refine=128` up already hits the 20000-iteration cap without
  early-stopping (see previous `gd_stress` entry) -- extrapolating that
  trend, `n_refine=1024` (~21M faces, 64x `n_refine=128`'s face count)
  could take on the order of ~10-11 hours single-threaded if it also runs
  the full cap. This is a rough extrapolation, not a calibrated estimate
  -- the point of this sweep is to actually measure it, including whether
  the linear-per-iteration-cost trend or the "hits the cap" behavior
  itself changes at this much larger scale (memory pressure, cache
  effects, etc. could break the extrapolation in either direction).
- **Submitted**: job-array `7305751`, 19 tasks. As of 01:25 EDT, 5/19
  running (`n_refine<=8` already finished and logged;
  `n_refine=12,16,24,32,64` running), remaining 14 queued behind the
  16-core-per-task reservation -- normal scheduling, not a problem.
  Results land in `campaign_runs/gd_stress_1024_2026-08-25/
  gd_stress_k<n_refine>.log`, one file per point,
  `manifest.json` in the same directory records the exact ladder/params/
  binary checksum.

**Both jobs left running unattended this session** -- user will check back
in a few hours. Next-session pickup: (1) confirm push7's two production
arrays actually got submitted (not just the prewarm) and check their
progress/completion against the ~25min/shard `n_refine=128` estimate from
the earlier entry; (2) read back `gd_stress_1024`'s per-point logs once
enough have finished to see whether the linear-in-`n_faces` per-iteration
cost trend and the "caps out at 20000 iters for n_refine>=128" behavior
both hold all the way to 1024, or break down at some point in between.

## 2026-08-25 (cont): push7 analyzed -- equal_area measurably accelerates
SO(3)-symmetry restoration vs naive, converging to the same plateau

Both push7 production arrays (`naive`, `equal_area`, 13 ladder points x 32
shards each, all 416/416 complete for both modes) finished sometime after
the last entry. Wrote `scripts/compare_push7_modes.py` to pool all 32
shards per `(n_refine, mesh_mode)` combo via `symmetry_test.analyze_multi`
and directly compare `kappa2_r(l)` (basis-independent SO(3)-breaking
cumulant) l=1..8 side by side. Output:
`campaign_runs/push7_2026-08-25_dual_naive_eqarea/plots/
naive_vs_eqarea_overlay{,_loglog}.png`.

**Headline result -- confirms the working hypothesis that motivated
building `equal_area` (2026-08-24 plateau finding)**: for the harmonics
most sensitive to the icosahedral point group (l=4,5,6,8 particularly --
`l=6` is the lowest angular momentum with a nontrivial icosahedral
invariant beyond the trivial l=0), `equal_area`'s `kappa2_r` drops to the
common high-resolution plateau *much faster in n_refine* than `naive`'s
does. In the intermediate range (n_refine 3-8) the gap is large and highly
significant:
- l=6, n_refine=4: naive=5.26e-3 vs eqarea=3.12e-4 (17x, 34 sigma)
- l=4, n_refine=3: naive=1.02e-3 vs eqarea=8.23e-5 (12x, 12.6 sigma)
- l=8, n_refine=4: naive=5.72e-3 vs eqarea=1.04e-3 (5.5x, 36 sigma)
- l=5, n_refine=3: naive=1.82e-3 vs eqarea=6.81e-4 (2.7x, 11.5 sigma)

By n_refine>=12-16 both mesh modes converge onto the *same* plateau
(~5e-5 to 1e-4, consistent within errors for every l tested) -- i.e.
`equal_area` does not change the asymptotic answer, it just gets there
at much lower (cheaper) resolution. l=1,2 (dipole/quadrupole) show no
resolvable naive-vs-eqarea difference at any resolution -- both sit in a
noisy ~1e-5 floor throughout, consistent with these low multipoles being
insensitive to icosahedral discretization defects (the icosahedral group
has no invariant subspace below l=6).

**Two points confirm the aliasing caveat flagged in the previous entry
(n_refine=2,3 below the informally-established l_max=8-safe floor of
n_refine>=4) and should be excluded from any clean mesh-mode read**:
- n_refine=2: naive and eqarea are statistically indistinguishable for
  every l>=3 (ratio~1.00, |sigma|<1) -- both modes are still essentially
  the bare base icosahedron at this refine level (GD has almost no room
  to move vertices), so this point can't discriminate mesh_mode at all,
  aliasing aside.
- n_refine=3, l=7: eqarea is *worse* than naive (naive=1.45e-3 vs
  eqarea=3.75e-3, -18 sigma) -- a sign flip relative to the l=4,5,6,8
  pattern above. l=8 at n_refine=3 is also still huge for both modes
  (0.16-0.22, order-1 symmetry breaking) -- not a real physics point for
  either mesh_mode at this l, just aliasing noise below the safe floor.

**Conclusion for the campaign**: `equal_area` should be the default
`--mesh_mode` going forward -- it reaches acceptable SO(3) restoration at
roughly half the n_refine (and therefore far less compute) that `naive`
needs for l up to 8. The residual common plateau (~5e-5-1e-4) at high
resolution is mesh-mode-independent, so it's not an area-uniformity
artifact -- next open question is whether it's a genuine finite-lattice
effect, a jackknife-bias floor, or something else; not investigated this
session.

## 2026-08-25 (cont): Owen_section_D diagnostic implemented -- Delta_s(a)
plateaus at ~0.33, not 1/8, flat across the whole push7 ladder

User pasted an excerpt (source paper/section not otherwise identified in
this repo; saved verbatim as `reference/Owen_section_D.tex`, referred to
as "Owen_section_D") titled "Agreement with the Ising CFT" and asked to
pick up the campaign's spherical-symmetry-check effort based on it. It
defines lattice Legendre coefficients

    F_l(a) = (2l+1)/(2*pi*4*pi) * sum_ij sqrt(g_i g_j) P_l(r_i.r_j) <s_i s_j>

(the *sum-over-m*, i.e. trace, version of the two-point-function harmonic
decomposition -- same closed-form recursion as Brower 2018 Sec 5.1
(`reference/brower_2018_lattice_phi4_riemann_S2.pdf`, already symbolically
verified in `scripts/derive_cl_closed_form.py`), different overall-
normalization convention) and two diagnostics built from it:

  - `Delta_s(a) = 2*F_1/(F_1+F_0)` -- direct lattice estimator of the spin
    scaling dimension (exact 1/8), from inverting the l=1 recursion step.
  - `delta_l(a) = 1 - (F_0-F_1)(F_{l-1}+F_l) / [l(F_0+F_1)(F_{l-1}-F_l)]`
    -- dimension-independent conformal/spherical-symmetry-breaking measure
    (Delta eliminated from the recursion), should -> 0 in the continuum
    limit.

**Key implementation point**: by the spherical-harmonic addition theorem,
`Trace(M_l) = sum_m M_l[m,m]` is *exactly* `(2*pi) * F_l(a)` up to a single
mesh-resolution-dependent (but l-independent) constant coming from this
driver's `/vol_sq` normalization (`vol = n_sites`, not a geometric area --
checked `include/S2.h`/`include/lattice.h`) -- and that l-independent
constant cancels exactly in both `Delta_s` and `delta_l` since they're
homogeneous ratios of F. So **F_l(a) is computable directly from the
already-collected `*_ylm_2pt_full_jackblocks_*.dat` files** used for the
existing `kappa2_r`/`R_l` test -- no new C++ measurement or production run
needed. New scripts: `scripts/cft_symmetry_test.py` (single-point
Delta_s/delta_l with jackknife errors, reusing `symmetry_test.py`'s
`read_jackblocks_multi`) and `scripts/run_cft_ladder.py` (ladder driver +
plots, mirrors `compare_push7_modes.py`'s pooling/plotting pattern).

**Correctness check before trusting any result** (directive.md's
"never trust a from-memory/hand-derived formula" rule): verified
`F_l = Trace(M_l)/(2*pi)` is internally consistent with the campaign's
*already* symbolically-verified `c_hat`/`derive_cl_closed_form.py`
machinery -- `Trace(M_l) = (2l+1) * c_hat[l]` (c_hat is the GLS per-m
diagonal estimate `symmetry_test.py` already computes), and since
`derive_cl_closed_form.py`'s own `lattice_conversion_factor` shows
`c_hat[l] ~ c_l^cont/(2l+1)`, the `(2l+1)` factors cancel and `F_l` is
l-independently proportional to Brower's `c_l^cont` -- so `Delta_s`/
`delta_l` computed from `F_l` should extract the same physics as the
existing (trusted) `c_hat`-ratio machinery, not something new. Hand-check
at `n_refine=128`: `6*c_hat[1]/(c_hat[0]+3*c_hat[1])` from the existing
`symmetry_test.analyze_multi` output reproduces `cft_symmetry_test.py`'s
`Delta_s` to 4 significant figures (0.3343 vs 0.334206). **Not a script
bug.**

**Result, run on push7's full 13-point ladder x {naive, equal_area}, 32
shards/point (`scripts/run_cft_ladder.py ... --ladder
2,3,4,6,8,12,16,24,32,48,64,96,128 --modes naive,equal_area`)**:
`Delta_s(a)` starts at 0.409 (n_refine=2), decreases monotonically through
n_refine=8 (0.350), then **flattens out around 0.33-0.34 from n_refine=16
through 128 and does not decrease further** (128: naive=0.332275(1230),
equal_area=0.334206(1173) -- statistically indistinguishable between mesh
modes, and both far from the exact value 0.125, ~15-30 sigma given the
~0.001 jackknife errors). `delta_l(a)` for l=4..8 shows the same story:
settles around -0.7 to -0.8 by n_refine~16-24 and stays there out to 128,
not trending toward 0 as the continuum-limit prediction requires. Full
table/plots: `campaign_runs/push7_2026-08-25_dual_naive_eqarea/plots/
cft_owen_secD_delta_s.png` and `..._delta_l.png`.

**This parallels, but is far more dramatic than, the still-unexplained
high-resolution `kappa2_r` plateau** (2026-08-24/25 entries above): both
diagnostics (i) are mesh-mode-independent (rules out area non-uniformity
as the cause, same conclusion as before) and (ii) stop improving well
before `n_refine=128`'s statistics are exhausted -- i.e. this is not a
statistics problem, it's a real systematic. Given `Delta_s(a)` is exactly
the estimator Owen_section_D's own source paper reports converging to
0.03%-level agreement with 1/8 by `L=64` on their lattice, a *flat*
plateau at 2.7x the correct value out to `n_refine=128` here is a strong
signal that something structural differs -- leading (untested) hypothesis:
PLAN.md's standing caveat that `exact_sinh`'s per-link, geometry-only
coupling assignment gives criticality "by construction" but this has never
been independently verified (no K_c scan is possible since there's no free
parameter) -- a small residual mass gap at the assigned coupling would
produce exactly this kind of saturating, non-vanishing `Delta_s(a)`/
`delta_l(a)`, concentrated most strongly at low l (large angular
separation, most IR-sensitive to a mass scale). Not confirmed -- no
alternative coupling/criticality check has been run. **This is now a
central open question for the campaign, arguably more urgent than the
pre-existing `kappa2_r` plateau**, since it bears directly on whether
`exact_sinh`'s "critical by construction" assumption (PLAN.md "Critical
coupling") actually holds. Next steps (not started): (a) check
`Delta_s(a)`/`delta_l(a)` against `production_2026-08-23`'s independent
l_max=8 dataset and the l_max=12 push5/push6 dataset as a cross-check;
(b) consider whether a small explicit mass counterterm / coupling rescan
is warranted despite `exact_sinh` having no nominal free parameter.

**Cross-check against `production_2026-08-23` (independent dataset — same
`exact_sinh`/`naive` combination as push7's naive arm, but a different push
entirely: `l_max=8` not derived from push7, different seeds, 16 shards not
32, ladder `{4,8,16,32,64,128}`)**: `Delta_s(a)` = 0.368/0.347/0.339/0.336/
0.336/0.334 at n_refine=4/8/16/32/64/128 respectively — **same plateau at
~0.334 by n_refine=32-64, same flattening pattern**, fully consistent with
push7's naive-mode numbers above (128: 0.332275(1230) push7 vs 0.334249
(1390) here, well within errors). This rules out the plateau being an
artifact specific to push7's run (different job, different seeds, only
`naive`, no `equal_area` in this dataset) — it is reproduced independently.
Confirms this is a real, robust systematic, not noise or a one-off bug in
a single production push.

**Sharper leading hypothesis (user prompted, 2026-08-25): coupling
assignment mismatch, not (necessarily) a generic "small mass gap"**.
`Owen_section_D.tex` references `Sec.~\ref{sec:WMcontinuum}` (not part of
the pasted excerpt) as where the source paper "establishes the
renormalization-group fixed point" for its couplings, and repeatedly calls
its lattice a "modified discretization of $S^2$" (see the Fig
`cft_break_mod`/`sphere_delta_sigma` captions). This strongly suggests
Owen's coupling assignment involved specific derivation/tuning work this
campaign has not replicated. This campaign's `exact_sinh` rule, by
contrast, is Twist's *flat*-triangular-lattice star-triangle result
(`K=0.5*asinh(cot(theta))`) applied per-link on a *curved* mesh via an
explicit, previously-flagged "judgment call" (two-face averaging, no
single-triangle-per-edge convention on S2 — see PLAN.md "Critical
coupling"/"Local coupling assignment", `directive.md` never covered this).
Criticality was asserted "by construction," never independently verified
against a Binder-crossing scan or any other check. **This plateau is
plausible direct evidence that assumption is wrong** — `exact_sinh` on the
curved mesh may carry a small residual relevant coupling (mass term) that
the flat-lattice derivation cannot see. Do not have access to Owen's
`WMcontinuum` section to compare directly; if the user supplies it, check
whether their coupling construction differs from `exact_sinh`'s two-face
average before assuming a generic "unspecified mass gap" story.

## 2026-08-25 (cont): ROOT CAUSE FOUND -- this is a mesh-uniformity problem,
not a coupling-formula bug. Owen's paper explicitly tested and rejected
equal-area meshing for exactly this reason.

User supplied the actual papers: R.C. Brower & E.K. Owen, "The Ising model
on S2", arXiv:2407.00459 (the full paper -- saved to
`reference/owen_2407.00459.pdf`) and its proceedings summary "Ising on S2 -
The Affine Conjecture", arXiv:2503.05621 (`reference/owen_2503.05621.pdf`;
user Rohan Misra is a co-author). Investigated by chasing why NONE of
`duality`/`exact_sinh`/`owen_dual` (a new rule added this session, see
below) fixed the `Delta_s(a)` plateau. Full chain of reasoning:

1. **Read 2503.05621 Eq 15/38**: the Ising coupling is
   `sinh(2*K_ij) = l*_ij/l_ij`, `l_ij` the direct chord edge length,
   `l*_ij` the "circumcenter dual length perpendicular to the edge". Added
   a new `--coupling_rule owen_dual` implementing this literally (using
   the already-existing `FaceCircumcenter()` for the two adjacent faces'
   circumcenters, `K = 0.5*asinh(|v_f0-v_f1|/len_l0)`) --
   `src/ising_s2_crit.cc`, rebuilt via `cluster/build.sh`. Quick
   login-node diagnostic (`n_refine` 4/8/16/32, 20000 traj, 1 seed,
   `campaign_runs/diagnostic_owen_dual_2026-08-25/`): `Delta_s(a)` =
   0.385/0.368/0.337/0.345 -- **statistically the same plateau as
   `exact_sinh`**, not fixed.
2. **Read 2407.00459 (the full paper) directly**, which clarifies
   `l*_ij` is explicitly NOT the straight 3d embedding-space distance
   between circumcenters (my `owen_dual` implementation) -- "the dual
   lattice edge lengths are not calculated by using the geodesic distance
   between triangle circumcenters in the 3-dimensional embedding space...
   the dual lattice lengths are computed by following intrinsic Regge
   geodesics as straight lines between flat triangular faces" (i.e.
   unfold the two adjacent flat faces about their shared edge into a
   common plane, since a 2d Regge manifold's curvature lives only at
   vertices, not edges). This intrinsic/"kink"-aware construction is Eq
   36-37 in 2407.00459, and matches this campaign's **pre-existing**
   `duality` rule almost exactly (its code comment already says "uses
   dual links on the trianglular lattice with a kink" -- inherited from
   the original `IsingS2/` codebase, predating this campaign, likely
   written by the original paper authors). Ran `duality` at the same
   quick-diagnostic resolutions
   (`campaign_runs/diagnostic_duality_check_2026-08-25/`): `Delta_s(a)` =
   0.354/0.365/0.346/0.349 -- **again the same plateau**. So the
   "correct" formula (already implemented, just never used for
   production per PLAN.md's step-1 status note -- `exact_sinh` was used
   instead per 2026-08-21 user direction) doesn't fix it either.
3. **All four candidate coupling constructions now tested
   (`duality`=Eq.37, `exact_sinh`=Twist's flat-lattice heuristic,
   `owen_dual`=naive/wrong embedding distance) give the same plateau** --
   this rules out "wrong coupling formula" as the (sole) explanation.
4. Cross-checked the quadrature-weight convention next (worried
   `OptimizeIntegrator`'s harmonic-exactness weights might differ from
   the paper's plain area weights) -- but `ising_s2_crit.cc` never calls
   `OptimizeIntegrator` at all, only `UpdateWeights()` (line 297), which
   is the standard circumcenter/cotangent Voronoi-dual-area formula --
   the *same* style of weight 2407.00459 Eq 50-52 uses (their `sqrt(g_i)`
   is a *name*, following the continuum `sqrt(g) d^2x` convention, for a
   per-site weight "proportional to the... circumcenter Delaunay dual
   area", normalized `sum_i sqrt(g_i) = 4*pi` -- not a literal square
   root of an area). `UpdateWeights` normalizes to `sum_i wt_i = n_sites`
   instead of `4*pi`, but since it's still proportional to the physical
   dual area with an l-independent constant, this cancels exactly in
   `Delta_s`/`delta_l` (homogeneous ratios) -- ruled out as the cause,
   not just plausible-but-untested.
5. **Root cause, confirmed directly from 2407.00459 Sec 3.1/3.3/Appendix
   B**: the paper *derives* (Eq 33-37) that `sinh(2K_ij)=l*_ij/l_ij` is
   only exactly critical (zero mass term `m_i` in their Eq 35) when
   **every triangle has equal circumradius AND equal perimeter** -- not
   equal area. Quoting directly: "the triangular faces of the basic
   discretization of the 2-sphere fail to satisfy the constraints of
   equal circumradius and equal perimeter... it is perhaps unsurprising
   that the higher l coefficients... do not converge to zero in the
   continuum limit" (Sec 3.3, describing their own *unmodified* /
   "basic"/`naive`-equivalent mesh -- i.e. this campaign's `naive` mode
   is reproducing exactly their own documented negative result). They
   built a dedicated Newton's-method mesh smoother (Appendix B.1/B.2)
   that explicitly minimizes `E = E_R + E_P` (Eq 39, normalized variance
   in circumradius and perimeter, `E_R = <R^2>/<R>^2 - 1`,
   `E_P=<P^2>/<P>^2-1`), using icosahedral-orbit-constrained barycentric
   degrees of freedom (not this campaign's straightforward vertex-
   position gradient descent). **Only after this specific smoothing does
   their `delta_l` measure go to zero and `Delta_s(a)` converge to
   0.125.** Critically, **they explicitly tested an equal-*area* mesh
   smoother too** (citing Karthik-Narayanan [20], same style as this
   campaign's `equal_area`/`EqualizeFaceAreas`) and report: "this
   construction fails to restore rotational symmetry in the continuum
   limit. The rotational symmetry breaking is especially strong using an
   octahedral base lattice... We therefore conclude that the definition
   of uniformity based on the circumradius and perimeter is indeed
   necessary." **This is a direct, first-party confirmation that this
   campaign's `equal_area` mesh mode was never going to fix the plateau**
   -- it targets the wrong non-uniformity metric, exactly as observed
   (push7 found `naive` and `equal_area` converge to the identical
   plateau).

**Conclusion**: the `Delta_s(a)`/`delta_l(a)` plateau is fully explained
and matches the source paper's own documented (and explicitly rejected)
failure mode. It is not a script bug (already cross-checked, see above
entry) and not a wrong-coupling-formula issue (four rules tested, all
plateau identically). **This campaign has never built the correct mesh
smoother** (Owen's Appendix B.1/B.2 circumradius+perimeter minimization,
via icosahedral-symmetry-orbit-constrained Newton's method) -- only
area-based smoothers (`EqualizeFaceAreas`, `equal_area_analytic`, both
already known/flagged as insufficient by prior sessions for unrelated
reasons, now additionally confirmed wrong-target by the source paper
itself). **This is the concrete, actionable next step for the campaign**:
implement circumradius/perimeter mesh uniformization per Appendix B.1/B.2,
not another area-based variant. Not started this session (substantial new
implementation -- orbit/symmetry-constrained DOF parameterization,
Newton's method with the finite-difference Hessian scheme described in Eq
57-58) -- flagging for explicit user go-ahead before starting given the
scope.

## 2026-08-25 (cont): implemented circumradius+perimeter mesh smoother
(`--mesh_mode equal_rp`) -- user directed to go ahead; area kept in the
joint objective per user direction ("the areas should still be equalized")

Added `QfeLatticeS2::EqualizeCircumPerim` (`include/S2.h`, right after
`EqualizeFaceAreas`) and wired it in as `--mesh_mode equal_rp`
(`src/ising_s2_crit.cc`; mesh-cache path suffix `_eqrp`, matching
`equal_area`'s `_eqarea` convention). Design choices:

- **Joint objective `E = E_R + E_P + E_A`** (`E_X = <X^2>/<X>^2 - 1` for
  circumradius R, perimeter P, and area A), not just `E_R+E_P` as
  arXiv:2407.00459 Eq 39 defines. Per explicit user direction: the
  paper's finding is that equal-area-*alone* fails, not that area
  uniformity is harmful when jointly minimized with R/P -- and this
  campaign's own push7 result (2026-08-25, earlier) already showed area
  equalization measurably speeds up SO(3)-symmetry restoration. Dropping
  it would have thrown away a confirmed win for no reason the paper
  actually gives.
- **Gradient computed via exact local finite differences**
  (`LocalEnergyWithVertexMoved`), not hand-derived closed-form R/P
  gradients. Perturbing a single vertex only changes the R/P/A of its
  (~5-6) incident faces, so the new exact global E can be recomputed in
  O(degree) by subtracting old per-face contributions and adding new ones
  (not an O(1/F)-dropped approximation like `EqualizeFaceAreas`'s
  documented shortcut -- this is the exact partial derivative). Chose
  finite differences deliberately over analytic differentiation of R (a
  more error-prone derivative than area, combining edge-length and
  triangle-area gradients) per directive.md's "never trust an unverified
  from-memory formula" rule -- a numeric gradient can't carry a silent
  sign/algebra bug.
- Same backtracking-line-search architecture as `EqualizeFaceAreas`
  (initial `step`, halve on rejected trial, grow 1.2x on accepted,
  `rel_tol`/`patience` early stop) -- reused wholesale, not reinvented.

**Verification before trusting it**:
- New standalone check `src/test_eqrp_symmetry.cc` (mirrors
  `test_ico_symmetry.cc`'s brute-force closure check: for every site and
  every icosahedral group element, confirm `G*r[s]` lands exactly on
  another site). Result: max mismatch 7.6e-9 at `n_refine=4`, 2.0e-7 at
  `n_refine=16`, 5.6e-6 at `n_refine=32` -- small but **not** machine
  precision like `EqualizeFaceAreas` achieves (~1e-15), because the
  finite-difference gradient (step `h=1e-6`) introduces O(h)-scale noise
  the analytic gradient doesn't have. Judged acceptable (still far below
  any physics-relevant scale) but noted honestly rather than claimed as
  exact.
- Sanity run at `n_refine=4,8,16,32`: circumradius std/mean drops from
  ~5-6% (pre-relax) to ~1-2% (post-relax); perimeter similarly ~6%->~1%;
  area also improves as a side effect of the joint objective (~12-13%->
  ~2-4%). Genuine convergence (checked `GD_DEBUG=1` iteration trace at
  n_refine=8: E plateaus at `7.45654e-4` with relative iteration-over-
  iteration changes down to 1e-13, not an early-stop artifact).

**Physics result, NOT yet resolving the plateau**: ran
`campaign_runs/diagnostic_eqrp_2026-08-25/` at `n_refine=4,8,16` (20000
traj, 1 seed, `--coupling_rule duality --mesh_mode equal_rp`, the
theoretically-correct pairing per the paper): `Delta_s(a)` =
0.361/0.344/0.332 -- **barely different from `naive`/`equal_area`'s
~0.33-0.36 at the same resolutions**, despite the ~5x reduction in R/P
non-uniformity. This is a genuinely open, not-yet-understood result:
either (a) `n_refine<=16` is still too coarse to see the effect (the
earlier naive-vs-`equal_area` comparison, a different non-uniformity
metric, also only diverged meaningfully past `n_refine~8-16` and
converged to a common plateau by `n_refine~12-16` -- `equal_rp` might
need to be pushed further, `n_refine=32,64`, before any trend shows), or
(b) the achieved `E` (~7e-4 at n_refine=8, ~2e-4 estimated at n_refine=16
from the printed std/means) is not yet small enough for this specific
unconstrained per-vertex gradient descent to reach what the paper's
orbit-constrained Newton's method achieves, or (c) something else is
still missing. `n_refine=32,64` diagnostics launched in background
(`campaign_runs/diagnostic_eqrp_2026-08-25/`, log
`/tmp/eqrp_k{32,64}.log` -- not yet moved into the repo) to check whether
the trend actually turns over at larger resolution before drawing a
conclusion either way. **Do not report `equal_rp` as a confirmed fix
until this trend is checked** -- current evidence is inconclusive, not
positive.

**Update**: the backgrounded `n_refine=32,64` `equal_rp` diagnostics (above)
finished: `n_refine=32` mesh relax converged in 5200/8000 iters (107.9s),
`n_refine=64` used the full 8000-iter cap (678.5s, so may not be fully
converged -- post-relax circumradius/perimeter std/mean 0.29%/0.17%, still
an improvement over 32's 0.39%/0.29% but the iteration cap being hit is
worth rechecking with a higher cap before trusting 64 as converged). Not
yet re-run through `cft_symmetry_test.py` at these two new points as of
this entry -- see the two entries below, which supersede this whole
plateau investigation with a sharper diagnostic and an independent
validation of the shared measurement code.

## 2026-08-25 (cont): deterministic free-scalar FEM test (`fem_scalar_test.cc`)
-- validates the shared mesh/quadrature/harmonic-projection code
independent of any Ising-specific assumption, and finds a real (but
separate) icosahedral-symmetry artifact at l=6

Built `src/fem_scalar_test.cc` (new small standalone driver, not touching
`ising_s2_crit.cc`) to directly test whether the campaign's shared
measurement infrastructure -- `QfeLatticeS2::UpdateWeights`'s cotangent-
Laplacian link/site weights, `UpdateYlm`/`GetYlm`, the `M_l[m,m']`
harmonic-projection convention -- could itself be the source of the
`Delta_s(a)` plateau, as an alternative to the "Ising coupling isn't
exactly critical on a curved mesh" hypothesis. Key insight: `UpdateWeights`
already computes exactly the standard cotangent-Laplacian FEM stiffness
weight per link (`half_wt = (sq_edge_1+sq_edge_2-sq_edge)/(8*FlatArea(f))`,
summed over the link's two faces -- the usual law-of-cosines cotangent
formula) and `sites[s].wt` is the matching lumped mass matrix (Voronoi dual
area), so a free massless scalar's lattice propagator can be built with
*zero new geometry code*: assemble the sparse graph Laplacian `K` from
`links[l].wt`, mass `M=diag(sites[s].wt)`, regularize `A=K+reg_eps*M`
(`reg_eps=1e-8`, far below the smallest physical eigenvalue `l(l+1)=2` at
l=1, to remove the constant zero-mode singularity without measurably
shifting any physical eigenvalue), factorize once with `Eigen::
SimplicialLDLT` (vendored sparse Cholesky), then for every `(l,m)` solve
`A x = v_lm` where `v_lm[i] = sites[i].wt * GetYlm(i,l,m)` and form
`M_l[m,m'] = v_lm^dagger . x_lm'` directly -- i.e. the exact discrete
Green's-function projection onto spherical harmonics, computed by linear
algebra alone. **No Monte Carlo, no thermalization, no statistical error
at all** -- this isolates the shared mesh/weight/harmonic code from every
Ising-specific question (coupling formula, criticality, autocorrelation).
Continuum comparison is exact and parameter-free: for the free massless
scalar, `C_l * l(l+1)` must be an l-independent constant for l>=1 (no
`Delta` to fit, no free normalization ambiguity beyond that overall
constant).

**Result (`q=5`, `l_max=8`, `n_refine in {2,4,8,16,32}`, both `naive` and
`equal_area` mesh modes)**: l=1,2,3,4,5,7,8 are clean and get *better* with
refinement in both mesh modes -- `C_l*l(l+1)` spread across those l shrinks
from ~15% at `n_refine=2` to <0.6% at `n_refine=32`, and `R_l` (off-diagonal
SO(3)-breaking power, same formula as `symmetry_test.py`) shrinks
monotonically by 3+ orders of magnitude (e.g. l=3: 7.7e-4 at n_refine=4 to
2.5e-7 at n_refine=32). **This is strong, independent evidence the shared
measurement infrastructure is not the source of the Ising `Delta_s(a)`
plateau** -- if `UpdateWeights`/`UpdateYlm`/the harmonic-projection
convention carried a bug, this deterministic test would show it too, and
for these l it doesn't.

**But l=6 is a permanent, order-1 outlier at every single resolution and
both mesh modes tried**: `R_6` = 1.85/1.85/1.85/1.80/0.29 (naive,
n_refine=2/4/8/16/32) and 0.92 (equal_area, n_refine=16) -- essentially
unchanged in order of magnitude across a 16x refinement in `n_sites`,
while every neighboring l sits at 1e-6..1e-7 by n_refine=16-32. Not noise,
not under-resolution (n_refine=32 has 10242 sites, well past this
resolution's aliasing-safe zone for l_max=8), and not fixed by `equal_area`
relaxation. **Root cause: exact icosahedral group theory, not a
discretization artifact.** The icosahedral rotation group (the exact
symmetry group of a `q=5` mesh -- not full SO(3)) has its own invariant
subspace of spherical harmonics starting at l=6 with multiplicity exactly
1 -- this matches `OptimizeIntegrator`'s own `n_l = (l/2)+(l/3)+(l/q)-l+1`
formula, which evaluates to `n_l=1` at `l=6, q=5` (the code has known about
this exceptional level all along, via the quadrature-exactness
requirement, without it ever being connected to the symmetry-test
diagnostics). A `q=5` mesh can support a nonzero icosahedrally-invariant
component at l=6 that **no mesh refinement or relaxation will ever remove**,
because it is topological to the mesh family's discrete symmetry group,
not a discretization error that shrinks as `a->0`. The same effect recurs
at higher icosahedral-invariant levels (l=10, 12, 15, ... -- rarer, not
checked here since `l_max=8`).

**Consequence for existing/future diagnostics**: this does not explain the
`Delta_s(a)` plateau itself (that estimator only uses `F_0,F_1`, i.e. l=0,1,
both clean in this test) -- the plateau is still open, see the entry below
for a sharper, more damning characterization of it. But it does mean
`delta_l(a)`'s l=4..8 sweep (`push7`'s `cft_owen_secD_delta_l.png`) and
`symmetry_test.py`'s `chi2`/`kappa2_r`/`R_l` at l=6 for any `q=5` run are
silently contaminated by a permanent, non-shrinking artifact unrelated to
either the coupling formula or a physical mass gap. Future ladder/fit
scripts for `q=5` meshes should explicitly exclude l=6 (and flag l=10,12,15
if `l_max` is ever pushed that high) rather than treating every l as an
equally trustworthy data point. Not yet done to `run_cft_ladder.py`/
`compare_push7_modes.py` -- open item.

Binary: `bin/fem_scalar_test` (build: same recipe as `cluster/build.sh`
minus the driver source file, e.g. `g++ -g -O3 -fopenmp -Wall
-Wno-deprecated-declarations -Wno-sign-compare -I include -I
include/unsupported -I "$SCC_BOOST_INCLUDE" -DGRP_DIR="\"$(pwd)/grp\""
src/fem_scalar_test.cc -o bin/fem_scalar_test`, after `module load
boost/1.83.0`). Usage: `bin/fem_scalar_test --q 5 --n_refine <k> --l_max
<lmax> --mesh_mode {naive,equal_area,equal_rp} [--equal_area_iters
--equal_area_step --mesh_cache_dir --reg_eps]`.

## 2026-08-25 (cont): user-prompted generalization of `Delta_s(a)` --
every adjacent-l pair should agree if this is genuinely one CFT primary,
and they do NOT: l=1 gives ~0.33, every l>=2 gives a completely different,
mutually-consistent ~0.58-0.63

User's objection to the existing `Delta_s(a) = 2*F_1/(F_1+F_0)` estimator:
it only ever inverts the recursion at the single step l=0->1 -- there is no
principled reason to trust that pair over any other adjacent pair `(l-1,l)`
if the ensemble really is described by one CFT primary of dimension
`Delta`. Added `delta_l_pair_formula` to `cft_symmetry_test.py`: the same
recursion inversion, `Delta_l_pair(a) = [rho*(l+1)-(l-1)]/(1+rho)` with
`rho=F_l/F_{l-1}`, evaluated at *every* `l=1..l_max`, not just l=1
(`delta_s_formula` is exactly the l=1 special case of this -- reduces
identically since `l-1=0`). If the ensemble is one CFT primary, every
`Delta_l_pair(a)` should agree with each other (and, in the true continuum
limit, with 1/8).

**They do not agree, and the disagreement is large and completely stable**.
Ran on `push7_2026-08-25_dual_naive_eqarea`'s full 32-shard production data
(`--coupling_rule exact_sinh` per that push's `manifest.json`, despite the
directory name -- `naive`/`equal_area` mesh modes) at `n_refine=32` and
`128`:

| l | Delta_l_pair, n_refine=32 (naive) | n_refine=128 (naive) | n_refine=32 (equal_area) | n_refine=128 (equal_area) |
|---|---|---|---|---|
| 1 | 0.3375(12) | 0.3323(12) | 0.3363(12) | 0.3342(12) |
| 2 | 0.5782(26) | 0.5809(26) | 0.5763(26) | 0.5768(25) |
| 3 | 0.6064(30) | 0.6052(30) | 0.6040(29) | 0.6041(30) |
| 4 | 0.6087(33) | 0.6112(33) | 0.6188(32) | 0.6198(33) |
| 5 | 0.6246(36) | 0.6213(37) | 0.6170(36) | 0.6167(35) |
| 6 | 0.6214(39) | 0.6266(40) | 0.6186(40) | 0.6191(39) |
| 7 | 0.6159(41) | 0.6203(43) | 0.6258(42) | 0.6224(43) |
| 8 | 0.6288(44) | 0.6225(46) | 0.6275(44) | 0.6278(46) |

(errors in parens, last digit(s), jackknife). **l=1 sits at ~0.33-0.34,
every l=2..8 clusters tightly at ~0.58-0.63 -- a factor-of-~2 split between
l=1 and everything else, essentially bit-for-bit identical between
n_refine=32 and 128 (i.e. not a shrinking discretization artifact) and
between `naive` and `equal_area` mesh modes (i.e. not fixed by the mesh
non-uniformity axis this campaign has spent the most effort on).** Also
cross-checked against the independent `equal_rp`/`duality` diagnostic
dataset (`diagnostic_eqrp_2026-08-25/`, n_refine=4/8/16/32/64, single seed,
noisier low-stats but same qualitative pattern: l=1 pinned at 0.33-0.36,
l>=2 scattered but consistently ~0.57-0.88, never near l=1's value) --
confirms this is not specific to one coupling rule or one production push.
**Notably l=6 is not an outlier in this diagnostic** (0.621-0.627 range,
squarely inside the l>=2 cluster) -- the icosahedral-invariant contamination
found in the FEM test's `R_l` (off-diagonal, m-dependent mixing) does not
obviously show up in this trace-based, m-summed `Delta_l_pair` estimator,
so the l=1-vs-l>=2 split is a different phenomenon, not explained by the
l=6 finding above.

**This reframes the whole plateau investigation.** The previous framing
("`Delta_s(a)` converges to ~0.33 instead of 0.125, a single wrong but
internally-consistent number") understated the problem: the l>=2 data
consistently point to a *different* apparent dimension (~0.58-0.63) than
l=1 does (~0.33), and neither is 0.125. A single global mass-gap/relevant-
coupling picture (the leading hypothesis so far) would need to explain why
the *shortest*-wavelength pairs (l>=2, more UV, "closer to critical" if a
mass term is an IR effect) plateau at a *larger* apparent Delta than the
longest-wavelength pair (l=0/1), which is the opposite of the naive
expectation that a mass term suppresses long-distance (low-l) correlations
and should bias low-l Delta upward, not downward, relative to high-l.
Possible explanations, none yet tested: (a) the two-point function on this
lattice is not well-described by a single CFT primary's closed form at
all -- e.g. it may be an admixture of the physical sigma operator plus a
lattice-artifact contribution with a different effective dimension that
dominates once l is large enough to resolve it, so no single `Delta` fit
should be expected to work uniformly across l; (b) the l=0,1 pair is
special/less trustworthy for a different reason than l>=2 (e.g. `F_0`'s
magnetization-squared-like zero mode has known special normalization
behavior at criticality -- flagged in PLAN.md's operator definition --
that the higher-l pairs don't share); (c) something about the coupling/
criticality assumption that specifically distorts the longest-wavelength
mode differently from all shorter-wavelength ones. **Not narrowed down --
this is now the sharpest concrete anomaly the campaign has, more
actionable than the single-number "0.33 plateau" framing was**, since it
rules out simple explanations (mesh non-uniformity, resolution, l=6
contamination) that could otherwise have been blamed. Code change:
`scripts/cft_symmetry_test.py` now reports `Delta_l_pair(a)` for every
`l=1..l_max` alongside the existing `Delta_s`/`delta_l` outputs.

## 2026-08-25 (cont, new session): push8 -- first production-scale
`equal_rp` run, generalized the two `equal_area`-only SGE task scripts to
also handle `equal_rp`

Per user direction to generate `equal_rp` meshes for `n_refine in
{4,6,8,12,16,24,32}`, run `--coupling_rule exact_sinh` production Ising on
them, and compute `Delta_l_pair(a)` across many adjacent-`(l-1,l)` pairs.
`--mesh_mode equal_rp` itself was already implemented and smoke-tested
earlier today (see the "implemented circumradius+perimeter mesh smoother"
entry above), but the two SGE task scripts under `cluster/sge/` were still
hardcoded to only special-case `equal_area`:

- `production_shard_task.sh` — `RUN_ID_SUFFIX`/`MESH_ARGS` logic extended
  with an `equal_rp` branch (`_eqrp` suffix, matching the driver's own
  `mesh_mode_suffix`); also now forwards `EQUAL_AREA_STEP` (previously
  only `EQUAL_AREA_ITERS` was passed through, silently leaving
  `--equal_area_step` at the driver's default even when a submit script
  set a different value via env).
- `prewarm_mesh_cache_task.sh` — was unconditionally hardcoded to
  `--mesh_mode equal_area` and the `_eqarea` cache-filename suffix.
  Generalized to take an optional `MESH_MODE` env var (default
  `equal_area`, so `push7`'s existing `-v` list keeps working unchanged)
  and compute the matching suffix, mirroring
  `production_shard_task.sh`'s logic exactly (a mismatch here would mean
  the prewarm writes a cache file the production run never looks for).

New `cluster/sge/submit_push8_eqrp.sh` (copied from `submit_push7_dual.sh`'s
structure): ladder `4:6:8:12:16:24:32`, `l_max=8` (matches
`production_2026-08-23`'s baseline, so `Delta_l_pair(a)` is directly
comparable), `--coupling_rule exact_sinh --mesh_mode equal_rp`, 32
shards/point, `n_therm=2000`/`n_traj=20000`/`n_skip=2`/`n_wolff=5`/
`n_metropolis=4` (same per-shard budget as `production_2026-08-23`),
`jack_block_size=100` (200 native blocks/shard, pooled 6400/point vs. the
`n_blocks>=170` gate), `equal_area_iters=20000`/`equal_area_step=0.3`
(driver defaults, generous cap relative to the `equal_rp` diagnostic's
observed convergence at `n_refine<=32`), `seed_base=800000` (disjoint from
every prior push). Ran the prewarm array first (`qsub -sync y`, job
7310279, 7 tasks) — completed cleanly, all 7 `q5k<k>_eqrp_step0.300.dat`
cache files confirmed present — then submitted the production array (job
7310315, 224 tasks = 7 ladder points x 32 shards) detached. Job 7310315
was still queued (`qw`) as of this entry; a background watcher
(`until ! qstat ... ; do sleep 30; done`) is tracking completion. Next
step once it finishes: run `scripts/cft_symmetry_test.py`/
`run_cft_ladder.py`'s `Delta_l_pair(a)` table across this ladder and
compare against the `naive`/`equal_area` numbers in the "user-prompted
generalization of `Delta_s(a)`" entry above (l=1 ~0.33 vs. l>=2
~0.58-0.63 there) — this is the first time that comparison will include
`equal_rp` at production statistics rather than the earlier single-seed,
`duality`-coupling diagnostic.

**Update, same session: job 7310315 finished (224/224 tasks exit 0, 0
failed) -- `equal_rp` does NOT fix the l=1-vs-l>=2 `Delta_l_pair(a)` split
at production statistics either.**

New `scripts/print_delta_l_pair_ladder.py` (pools all 32 shards' jackblocks
per ladder point via `cft_symmetry_test.analyze`, prints the full `l x
n_refine` table -- `run_cft_ladder.py` only plots `delta_s`/`delta_l`, not
the full `delta_l_pair` table, so this fills that gap; also added
`equal_rp`/`_eqrp` to `run_cft_ladder.py`'s `MODE_SUFFIX` map for future
use with the plotting path). Result on `push8_2026-08-25_eqrp_exactsinh`
(`--coupling_rule exact_sinh --mesh_mode equal_rp`, l_max=8, ladder
`{4,6,8,12,16,24,32}`, 3200 pooled jackknife blocks/point):

| l | k=4 | k=6 | k=8 | k=12 | k=16 | k=24 | k=32 |
|---|---|---|---|---|---|---|---|
| 1 | 0.3687(9) | 0.3542(9) | 0.3479(9) | 0.3439(9) | 0.3415(9) | 0.3379(9) | 0.3374(9) |
| 2 | 0.5829(20) | 0.5777(21) | 0.5801(21) | 0.5774(20) | 0.5778(21) | 0.5766(21) | 0.5757(21) |
| 3 | 0.6327(23) | 0.6162(23) | 0.6075(24) | 0.6041(24) | 0.6018(23) | 0.6103(23) | 0.6066(24) |
| 4 | 0.6709(26) | 0.6372(26) | 0.6241(26) | 0.6142(26) | 0.6153(26) | 0.6103(26) | 0.6122(26) |
| 5 | 0.7188(29) | 0.6622(28) | 0.6417(30) | 0.6304(29) | 0.6183(29) | 0.6168(29) | 0.6193(29) |
| 6 | 0.7853(30) | 0.6835(31) | 0.6465(31) | 0.6317(31) | 0.6279(31) | 0.6223(32) | 0.6221(31) |
| 7 | 0.8494(33) | 0.7109(34) | 0.6760(34) | 0.6330(33) | 0.6325(33) | 0.6282(34) | 0.6181(35) |
| 8 | 0.9325(36) | 0.7470(36) | 0.6807(36) | 0.6503(36) | 0.6378(36) | 0.6195(36) | 0.6285(36) |

(errors in parens, units of 1e-4.) **At the two largest, best-resolved
points (k=24,32) the pattern is bit-for-bit the same split found on
`naive`/`equal_area`**: l=1 -> 0.3379/0.3374, l=2..8 all cluster tightly
at ~0.58-0.63 (0.5757-0.6285) -- indistinguishable in shape from the
`naive`/`equal_area` table in the entry above. The l=6,7,8 blowup at k=4
(0.79/0.85/0.93) is the pre-flagged `l_max=8`-at-`n_refine=4` aliasing
risk (81 modes vs 162 sites) this push's header explicitly called out,
not a new phenomenon -- it recedes monotonically and is gone by k=12-16.

**Conclusion: circumradius+perimeter mesh uniformity (this campaign's
leading hypothesis after the naive/equal_area root-cause finding) is now
ruled out too, at real production statistics with the campaign's primary
`exact_sinh` coupling rule.** Three independent non-uniformity axes have
now been tried and failed to close the l=1-vs-l>=2 split: mesh resolution
(32 vs 128), face-area uniformity (`equal_area`), and joint
circumradius+perimeter+area uniformity (`equal_rp`) -- plus multiple
coupling rules (`duality`, `exact_sinh`, `owen_dual`) from the earlier
root-cause investigation. This reopens the field of candidate
explanations listed in the "user-prompted generalization of `Delta_s(a)`"
entry above ((a) admixture of a second, non-CFT-primary lattice operator,
(b) l=0/1's `F_0` zero-mode having special normalization behavior, (c) a
coupling/criticality distortion specific to long-wavelength modes) --
none of the three have been tested yet. **This is now the central open
question for the campaign, superseding the mesh-uniformity avenue
entirely** -- do not propose another mesh-smoother variant without first
addressing why the split is not a uniformity effect at all (it survives a
factor-of-2 shrink in circumradius/perimeter non-uniformity with
essentially no motion in either the l=1 or the l>=2 cluster's plateau
value).

## 2026-08-25 (cont): light Binder-cumulant (U4) scan -- bulk order
parameter looks genuinely scale-invariant, unlike `Delta_l_pair(a)`,
pointing suspicion at the harmonic-decomposition analysis code rather
than the physics/coupling rule

Per user request, a lightweight cross-check of criticality independent of
the `S_lm`/`M_l[m,m']` harmonic-decomposition machinery entirely: the
Binder cumulant `U4 = 1.5*(1 - <m^4>/(3*<m^2>^2))` from the existing
`*_bulk_*.dat` files (already written by every production run, no new
simulation needed). New `scripts/binder_scan.py`: `QfeMeasReal::
WriteMeasurement`'s on-disk format is `mean, err, n` with `err` the
*naive* (non-autocorrelation-corrected) standard error
(`include/statistics.h`'s `Error()`), so `sum`/`sum2` can be exactly
reconstructed per shard (same inverse formula `ReadMeasurement` itself
uses) and pooled by summation across shards -- exact, not
error-propagated from already-averaged per-shard means. Errors via
leave-one-shard-out jackknife (shard-granularity, not `jack_block_size`
granularity -- deliberately "light", no new file format).

Ran on `push7_2026-08-25_dual_naive_eqarea`'s `naive`-mesh ladder (widest
available resolution range, `n_refine=2..128`, all `exact_sinh`, 32
shards/point, 208,000 pooled traj/point):

| n_refine | U4 | err |
|---|---|---|
| 2 | 0.820627 | 0.000680 |
| 3 | 0.831643 | 0.000613 |
| 4 | 0.834468 | 0.000547 |
| 6 | 0.840047 | 0.000644 |
| 8 | 0.843269 | 0.000460 |
| 12 | 0.846326 | 0.000720 |
| 16 | 0.847193 | 0.000620 |
| 24 | 0.849248 | 0.000506 |
| 32 | 0.848572 | 0.000542 |
| 48 | 0.849661 | 0.000561 |
| 64 | 0.848975 | 0.000646 |
| 96 | 0.851285 | 0.000572 |
| 128 | 0.852182 | 0.000598 |

**Reading**: steep rise `n_refine=2->16`, then a much shallower creep
`n_refine=32->128` (0.8486->0.8522, ~5 sigma over a further 4x
refinement but small in absolute size) -- consistent with a
slowly-decaying finite-size/curvature correction approaching a plateau,
not a runaway toward a trivial value (1 for a fully ordered/frozen state,
2/3 for a Gaussian/non-interacting fixed point). The asymptotic value
(~0.85, still creeping up slightly at the largest resolution tried) is in
the right ballpark for the known 2D Ising Binder-cumulant universal
amplitude (~0.856-0.861 depending on geometry/boundary conditions) --
this bulk order-parameter diagnostic looks like a genuinely
scale-invariant, nearly-critical statistical field theory, i.e. what a
lattice realization of the 2D Ising CFT on S^2 should look like.

**This sharpens the interpretation of the `Delta_l_pair(a)` anomaly.** A
resolution-independent, coupling-rule-independent, mesh-uniformity-
independent l=1-vs-l>=2 split (see the several entries above) in a system
whose bulk order parameter behaves like a genuinely critical theory is
now better explained by a bug or convention error in the harmonic-
decomposition analysis pipeline (`S_lm`/`M_l[m,m']` measurement in
`ising_s2_crit.cc`, or the `F_l`/`Delta_l_pair` derivation in
`scripts/cft_symmetry_test.py`) than by the physical system failing to be
critical -- hypothesis (b) from the entry above (`F_0`'s zero-mode having
special normalization behavior, since l=0/1 is exactly where the split's
"outlier" sits) is now the leading candidate, with (a) (operator
admixture) still open and (c) (coupling-specific long-wavelength
distortion) now disfavored by this Binder result. **Next step: deep code
audit of the harmonic-projection/F_l/Delta_l_pair pipeline**, not another
physics run -- see follow-up entry.

## 2026-08-25 (cont): deep code audit of the F_l/Delta_l_pair pipeline --
found and fixed two genuine documentation bugs, found one harmless
convention mismatch, found NO bug that explains the l=1-vs-l>=2 split

Per user request, given the Binder scan's evidence the bulk theory looks
genuinely critical (previous entry), audited every step from raw MC spin
configurations to the printed `Delta_l_pair(a)` table, checking each
piece against `reference/Owen_section_D.tex`'s stated formulas
line-by-line rather than trusting the campaign's own prior derivation.

**1. `F_l(a) = Trace(M_l)/(2*pi)` (the core formula `cft_symmetry_test.py`
uses) is algebraically exactly right.** Re-derived from the paper's
`F_l(a) = 1/(2*pi) * sum_ij sqrt(g_i g_j) (sum_m Y*_lm(r_i)Y_lm(r_j))
<s_i s_j>` and the campaign's own `S_lm = sum_i w_i Y*_lm(x_i) s_i`,
`M_l[m,m'] = <S_lm S*_lm'>` definitions: `Trace(M_l) = sum_ij w_i w_j
(sum_m Y*_lm(x_i)Y_lm(x_j)) <s_i s_j>`, identical to `2*pi*F_l(a)` term
by term once `sqrt(g_i g_j) = w_i w_j` is used (the campaign's own stated
identification). Also re-derived the paper's `F_0(a) = 2*<M(a)^2>`
identity directly from `Y_00 = 1/sqrt(4*pi)` and confirmed it holds
exactly (not just in the continuum limit) at the level of the formula.

**2. `Delta_l_pair_formula`'s algebra is exactly right.** Solved
`rho = F_l/F_{l-1} = (l-1+Delta)/(l+1-Delta)` for `Delta` by hand:
`Delta = [rho*(l+1)-(l-1)]/(1+rho)`, term-for-term identical to the code.

**3. `S_lm` array indexing (`l*(l+1)/2 + |m|`) is correct** -- matches
the triangular enumeration order the `ylm` cache is filled in
(`ising_s2_crit.cc`'s `y_i/y_l/y_m` loop), verified by hand for l=0,1,2.

**4. Found: `ising_s2_crit.cc` never calls `OptimizeIntegrator(l_max)`**
(the harmonic-exact quadrature weight solver `CLAUDE.md` and this
module's docstring both claimed was in use) -- **it calls
`lattice.UpdateWeights()`** instead, the ordinary cotangent-Laplacian
FEM/lumped-Voronoi-dual-area weight scheme (confirmed via grep: zero
call sites for `OptimizeIntegrator` anywhere outside its own definition
and one `n_l`-formula comment). `UpdateWeights` normalizes
`sum_i sites[i].wt = n_sites` (`vol`), not `4*pi`. **This is a real,
confirmed documentation bug** (fixed in `CLAUDE.md` and this module's
docstring this session) **but does not explain the anomaly**: the
campaign's own `fem_scalar_test.cc` (2026-08-25, earlier entry) calls
the *same* `UpdateWeights` -- not `OptimizeIntegrator` either, confirmed
by re-checking its source this session -- and still comes back clean for
l=1,2,3,4,5,7,8 on a deterministic free-scalar propagator with zero
statistical error. If `UpdateWeights`'s non-harmonic-exact quadrature
caused l-dependent leakage/aliasing of the kind seen in
`Delta_l_pair(a)`, the deterministic scalar test would show it too
(no MC noise to hide behind) -- it doesn't, ruling this out as the cause
rather than confirming it.

**5. Found: `ising_s2_crit.cc`'s `S_lm` accumulation uses `Y_lm`, not the
conjugate `Y*_lm`** this module's docstring (and the paper) state --
`ylm(s,y_i) = lattice.CalcYlm(s,y_l,y_m)` returns
`boost::math::spherical_harmonic(l,m,...)` directly, un-conjugated.
Worked through the algebra: letting `X_lm = sum_i w_i Y_lm(x_i) s_i` be
what the code actually stores (`S_lm_paper = conj(X_lm)` in the paper's
convention), the code's stored pair product is `X_lm * conj(X_lmp)`,
which equals `conj(M_l[m,mp]_paper)` for every `(m,mp)`. Since `conj(z)`
has the same magnitude and, on the diagonal `m=mp`, the same real value
as `z`, **this leaves every trace-based (`F_l`) and magnitude-based
(`R_l`, `kappa2_r`/`kappa3_r`/`kappa4_r`) diagnostic this campaign
actually uses completely unaffected** -- it only flips the phase of
individual off-diagonal `M_l[m,m']` (m!=m') entries relative to the
stated convention, which no current script reads directly (only their
magnitudes/traces). Benign, documented in this module's docstring this
session, not fixed in C++ since nothing depends on the phase.

**6. Found: an overall, l-INDEPENDENT normalization factor of
`1/(4*pi)^2` between this module's `F_l` and the paper's absolute `F_l`
convention**, from combining item 4's `sum wt = n_sites` (not `4*pi`)
with `ising_s2_crit.cc`'s explicit `/vol_sq` (`vol = n_sites`) division
before writing each configuration's `S_lm*conj(S_lmp)` product to the
jackblocks file. Verified algebraically this factor is exactly the same
at every `l` (traces back to a uniform rescaling of every `S_lm`, not an
l-dependent one), so **it cancels completely in `Delta_l_pair`/`Delta_s`/
`delta_l`, which are built entirely from ratios of `F_l`'s** -- confirmed
not the source of the split either. Only consequence: this module's
printed `F_0` is not literally `2*<M(a)^2>` in absolute units (off by
`(4*pi)^2 ~ 157.9`) -- flagged in the docstring, not otherwise
consequential.

**7. Estimated the coincident-point ("self-contact", `i=j`) term's size**
and found it far too small to matter: `P_l(1)=1` for every l, so the
`i=j` terms contribute `(2l+1)/(4*pi) * sum_i w_i^2` to `Trace(M_l)`
before normalization -- grows linearly with `l` (a real, `l`-dependent
lattice artifact, unlike items 4-6), but after the `/n_sites^2`
normalization this is `~(2l+1)/(8*pi^2*n_sites)`, order `2e-5` at
`l=8, n_refine=32` (`n_sites=10242`) -- three to four orders of magnitude
below the measured `F_l` values themselves (`~1e-4` to `~4e-3`) and far
too small to produce an `O(0.1-0.3)` `Delta_l_pair` discrepancy that
holds at full production statistics.

**Conclusion: no bug found in the `F_l`/`Delta_l_pair` analysis code (C++
measurement or Python post-processing) that explains the l=1-vs-l>=2
split.** Every formula checked out algebraically against
`reference/Owen_section_D.tex`; the two real discrepancies found (items 4
and 6) are both provably l-independent and therefore invisible to every
ratio-based diagnostic the campaign actually uses; the one l-dependent
lattice artifact identified (item 7, coincident-point contamination) is
three-plus orders of magnitude too small. Combined with the Binder scan
showing the bulk theory looks genuinely critical, this now weighs
**against** a code bug or a bulk-non-criticality explanation and **toward
hypothesis (a) from the "user-prompted generalization" entry above**: a
genuine admixture in the interacting Ising two-point function itself --
e.g. a lattice-artifact/short-distance operator contribution with a
different effective dimension mixing additively into `<s_i s_j>`, on top
of the physical sigma primary -- something the deterministic
`fem_scalar_test.cc` (a *free*, non-interacting scalar, no criticality
tuning, no lattice-artifact operators to admix) structurally cannot probe
even though it shares the exact same mesh/weight/harmonic-projection
code. Hypothesis (b) (`F_0`'s zero mode having special normalization)
remains formally untested but is now less favored, since items 4-6 show
`F_0`'s measurement machinery has no privileged bug relative to any other
`l`'s. **Recommended next step, not yet done**: directly test hypothesis
(a) by fitting `<s_i s_j>` (or equivalently `F_l(a)`) to a two-term model
(CFT primary + a second power/operator) rather than a single `Delta`, or
by checking whether the l=2..8 cluster's apparent dimension (~0.6) is
consistent with any known lattice-artifact/irrelevant-operator dimension
for this discretization.

## 2026-08-25 (cont): the per-l Delta_l_pair(a) values are individually
STABLE, precisely-converging numbers -- the anomaly is that there are
(at least) two different stable numbers, not noise or non-convergence

Per user pushback ("you should be able to get the deltas, why are they
not stable") -- re-examined push7's full 13-point naive-mesh ladder
(`n_refine=2..128`, `exact_sinh`, 416,000 pooled traj/point) instead of
just eyeballing the n_refine=32/128 snapshot from the entry above.
Weighted least-squares fit of `Delta_l_pair(a)` to `a_inf + b/n_refine`
per l (`n_refine>=8` only, to stay clear of the known aliasing-risk
region at small `n_refine`):

| l | extrap Delta_l_pair(a->0) | slope b (1/n_refine coeff) | chi2/dof |
|---|---|---|---|
| 1 | 0.3322 +/- 0.0006 | +0.141 +/- 0.010 | 0.98 |
| 2 | 0.5790 +/- 0.0013 | +0.008 +/- 0.023 (~0) | 0.94 |
| 3 | 0.6055 +/- 0.0015 | -0.020 +/- 0.026 (~0) | 0.60 |
| 8 | 0.617 +/- 0.005 (fit restricted to n_refine>=32) | +0.35 +/- 0.24 | 0.21 |

**Every one of these fits is clean** (chi2/dof near 1, tight errors) --
this is not a statistics or convergence problem, contrary to how the
"plateau" framing in earlier entries could be read. The real finding:
l=1 and l>=2 extrapolate to two **different, precisely-determined**
numbers, not to a shared badly-converged one. l=2 and l=3 are already
flat (slope consistent with zero) by `n_refine=8` -- fully converged,
stable, and equal to ~0.58-0.61, not 0.125. l=1 is NOT yet flat even at
`n_refine=128` -- it has a real, resolvable **O(1/n_refine) = O(a)**
correction (10+ sigma nonzero slope), extrapolating cleanly to
0.3322(6). l=8 also shows a real slope (noisier, aliasing-adjacent) but
even restricted to its best-resolved points (`n_refine>=32` only)
extrapolates to ~0.617(5) -- close to, but several sigma above, l=2's
0.579(1): **even within the "l>=2 cluster" there is a small but real,
resolvable drift with l** (0.579 -> 0.606 -> ~0.617 from l=2 to 3 to 8),
not one exactly-shared constant, though far tighter than the l=1-vs-l>=2
gap.

**Two things stand out as real physics clues, not artifacts:**
1. **l=1's O(a) (not O(a^2)) lattice correction is the wrong power** for
   an ordinary irrelevant lattice-artifact operator (those correct as
   `a^2` or faster on a properly built discretization) -- O(a) is the
   signature expected from a *relevant or marginal* operator/scale
   contaminating specifically the longest-wavelength (l=1) channel, not
   generic discretization error.
2. l>=2's near-zero slope means those channels are essentially UV-
   insensitive already at coarse `n_refine` -- whatever contaminates
   them saturates quickly and does not shrink with further refinement,
   unlike l=1's contamination which is still visibly draining away.

**Updated framing for the campaign**: `<s_i s_j>` is not well-described
by a single power law `(2-2z)^{-Delta}` at ANY resolution tested here --
not "not yet converged", but converged to two-plus distinct effective
dimensions depending on which harmonic channel probes it. Combined with
the code audit (previous entry, no bug found) and the Binder scan (bulk
theory looks genuinely critical), the most concrete next step is now a
genuine two-operator fit of `F_l(a)` (CFT primary of unknown Delta +
second term, jointly fit across all l using the individually-precise
per-l data above as the target, not a single global Delta) -- not
another mesh/coupling/code investigation.

## 2026-08-25 (cont): two-operator fit of F_l(a) -- decisively rejects a
single power law, gives a stable second dimension ~0.63-0.65, gives a
much-improved but not yet exactly-1/8 first dimension

Per user direction, did the two-operator fit the previous entry
recommended: `F_l = A*F_l^cont(Delta_1) + B*F_l^cont(Delta_2)` (linear
superposition of two copies of Owen_section_D's exact recursion, since
the Legendre decomposition is linear in the correlator), fit via
`scipy.optimize.least_squares` (Levenberg-Marquardt, multiple initial
guesses checked for a consistent global minimum -- all guesses converged
to the same solution at every ladder point tested) against
`push7_2026-08-25_dual_naive_eqarea`'s naive-mesh `F_l(a)` (l=0..8,
`exact_sinh`, diagonal chi-square using each `F_l`'s own jackknife error,
**not yet the full jackknife covariance across l** -- important caveat,
see below).

**At n_refine=128 (finest single point): chi2/dof = 0.42 (9 data points,
4 params, dof=5), Delta_1 = 0.0724(35), Delta_2 = 0.6364(14)** -- versus
a single-power-law fit's chi2/dof = 7569 (catastrophic, decisively
rejected) at the same data. Repeated across the ladder (`n_refine=8..128`
using n_refine>=16 for reliable convergence, low-n_refine points are
severely aliased and give chi2/dof in the hundreds-to-thousands):

| n_refine | Delta_1 | Delta_2 | chi2/dof |
|---|---|---|---|
| 8 | 0.160(18) | 0.698(11) | 16.6 |
| 12 | 0.110(9) | 0.658(4) | 3.7 |
| 16 | 0.0997(34) | 0.6487(15) | 0.48 |
| 24 | 0.0780(34) | 0.6385(14) | 0.45 |
| 32 | 0.0754(50) | 0.6358(20) | 0.90 |
| 48 | 0.0764(44) | 0.6355(18) | 0.69 |
| 64 | 0.0643(24) | 0.6305(9) | 0.20 |
| 96 | 0.0679(43) | 0.6337(17) | 0.68 |
| 128 | 0.0724(35) | 0.6364(14) | 0.42 |

**Delta_2 is rock-stable at ~0.63-0.65 from n_refine=16 on** (small,
resolvable drift, similar in size to the Delta_l_pair cluster's own drift
noted in the entry above). **Delta_1 is much closer to the true 1/8 than
any raw `Delta_l_pair` estimate got** (0.06-0.10 vs. 1/8=0.125, compare
to the single-operator l=1 estimate's clean-but-wrong extrapolation to
0.332) but does not yet land exactly on 1/8: a naive `a_inf + b/n_refine`
extrapolation of `Delta_1(a)` over n_refine>=16 gives `a_inf =
0.0591(22)`, notably *below* 0.125 (30-sigma at face value), with the
extrapolation fit itself imperfect (chi2/dof=2.46, i.e. the 1/n_refine
form doesn't fully capture the curvature -- the same crude ansatz that
underperformed for the single-l1-operator extrapolation earlier in this
file).

**Two caveats, not yet resolved, before trusting the fitted `Delta_1`
value precisely**:
1. This fit used a **diagonal chi-square** (`F_l`'s own jackknife error
   per l, ignoring cross-l covariance) even though every `F_l` at a given
   ladder point is computed from the *same* jackknife blocks and is
   therefore correlated with every other `F_l` at that point -- a proper
   fit needs the full `l x l` jackknife covariance matrix (straightforward
   to build from the same `*_ylm_2pt_full_jackblocks_*.dat` files
   `cft_symmetry_test.analyze` already reads, not yet implemented). This
   could shift both the fitted central values and, especially, the error
   bars in either direction.
2. The `a_inf + b/n_refine` continuum-extrapolation ansatz was already
   flagged as imperfect for the single-operator l=1 case (previous
   entry) and is even less obviously justified for a fitted *nonlinear*
   parameter of a two-term model -- worth trying `a_inf + b/n_refine +
   c/n_refine^2` or extrapolating `A`, `B`, `Delta_1`, `Delta_2` jointly
   with a proper multi-point simultaneous fit instead of point-by-point
   fit-then-extrapolate.

**Bottom line**: the two-operator model is decisively preferred over a
single CFT primary (chi2/dof improves by 3-4 orders of magnitude) and
gives a `Delta_1` far closer to 1/8 than the naive single-l estimators
did, with a stable companion `Delta_2 ~ 0.63-0.65` -- **but does not yet
cleanly hand back exactly 1/8 in the continuum limit**, so this is a real
step forward in modeling, not yet a resolution of the campaign's central
open question. Next steps, not yet done: full jackknife covariance in
the fit (item 1 above), a better continuum-extrapolation functional form
(item 2), and joint (not point-by-point) fitting across the whole ladder.
Ad hoc scripts used: `/tmp/two_op_fit.py`/`/tmp/two_op_ladder.py`
(session scratch, not committed to the repo -- rewrite as a proper
`scripts/two_operator_fit.py` before the next session if this direction
is pursued further, per the "no scratchpad for derived data" directive.md
rule -- these are analysis code, not derived data, but should still live
in the repo to be reproducible).

## 2026-08-25 (cont): BREAKTHROUGH -- direct real-space two-point-function
fit gives a clean Delta ~ 0.128, close to the exact 1/8, entirely bypassing
the F_l/Legendre-recursion machinery. The whole `Delta_l_pair`/two-operator
anomaly is an artifact of the l-space method, not the physics.

Per explicit user direction ("extract the scaling dimension from the two
point functions rather than the recursion relations" -- user's stated
suspicion: "there must be an issue in the analysis"), built a genuinely
independent measurement that bypasses `S_lm`/`M_l[m,m']`/spherical
harmonics entirely.

**New standalone driver `src/real_space_2pt_test.cc`** (mirrors
`fem_scalar_test.cc`'s precedent of a separate diagnostic binary, not a
change to the production driver): reuses `ising_s2_crit.cc`'s
mesh-construction/coupling-assignment code verbatim (kept in sync by
hand), then instead of any harmonic projection, picks `n_ref` random
reference sites, precomputes which of `n_bins` angular bins
(uniform in `z=cos(theta)`, i.e. equal-solid-angle rings) every other site
falls into relative to each reference site, and per configuration
accumulates `sum(s_ref * s_target)` per bin -- a direct, real-space
`<s(theta)>` estimator with **zero spherical-harmonic machinery, no
quadrature weights, no l_max truncation**. `O(n_ref * n_sites)` per
configuration (not the full `O(n_sites^2)` all-pairs cost), so this runs
at full production `n_refine` cheaply: n_refine=32 (10242 sites),
n_ref=300, n_bins=100, 20000 traj/2 skip = 10000 measurements took ~2
minutes.

New `scripts/fit_real_space_2pt.py`: fits the binned `<s(theta)>` output
directly to `A*(2-2*cos(theta))^{-Delta}` via `scipy.optimize.least_squares`
in a chosen `[theta_min, theta_max]` window.

**Result at n_refine=32 (naive mesh, exact_sinh, seed 42, 10000 measured
configs)**: excluding only the smallest angular separations
(`theta_min>=0.3`, i.e. cutting close-neighbor short-distance bins),
**Delta = 0.128(0) with chi2/dof = 0.02-0.03**, essentially flat across
every window tried from `theta_min=0.3` to `0.8` and `theta_max=2.6` to
`3.0` (Delta ranged 0.1275-0.1284 across 8 different window choices, all
excellent fits). **This is dramatically different from every l-space
result this campaign has produced** (Delta_s~0.33, Delta_l_pair l>=2
cluster~0.6, two-operator fit's Delta_1~0.06-0.10) and is close to the
exact CFT value 1/8=0.125. Confirmed independently at n_refine=16 (300
seconds smoke-test-scale stats): `Delta=0.1287(3), chi2/dof=0.020` --
same conclusion at a different resolution. Only including the very
closest-neighbor bins (`theta_min<0.3`) degrades the fit
(chi2/dof jumps to 4, Delta drifts to ~0.134) -- consistent with genuine,
localized short-distance lattice-discretization contamination that a
real-space fit can simply exclude with a window, unlike any l-space
estimator.

**Direct comparison confirms the l-space method itself is distorted, not
just under-precise**: resummed the harmonic-MC-measured `F_l(a)` (push7's
n_refine=32 naive dataset, exact_sinh, same physical parameters as the
real-space run) into a truncated real-space curve
`C_resum(theta) = sum_{l=0}^{8} (2l+1)/2 * F_l * P_l(cos(theta))` and
compared pointwise against the independently-measured real-space
correlator at matched theta. **The two disagree by a factor of ~100-250x,
and the ratio is NOT constant across theta** (varies roughly 0.004 to
0.017 across the 7 angles checked) -- i.e. not merely the previously-
identified, l-independent `1/(4*pi)^2` normalization convention mismatch
(item 6 of the earlier code-audit entry, `~1/158`, in the right ballpark
for the average ratio but not exactly it) but a genuine **shape**
distortion between the l_max=8-truncated harmonic reconstruction and the
true correlator, worst near small theta (largest deviation at the
smallest theta checked, 0.295).

**Mechanism (not yet independently verified, but the natural explanation
given all evidence collected today)**: every Legendre coefficient `F_l`,
even `F_0`/`F_1`, is a *global* integral over the *entire* sphere,
weighted by `P_l(cos theta)` which is O(1) at every theta including the
short-distance region -- there is no way for a finite set of `F_l`'s
(or any recursion/ratio built from them) to *exclude* the short-distance,
lattice-discretization-contaminated region the way a real-space fit can
simply do with a `theta_min` cut. If the actual lattice `<s_i s_j>`
genuinely deviates from a pure `(2-2z)^{-1/8}` power law only near
`theta` comparable to the lattice spacing (an entirely ordinary,
expected lattice artifact for ANY discretized QFT), that contamination
still leaks into every single `F_l`, and a recursion/ratio-based Delta
estimator built from those contaminated `F_l`'s has no way to see or
exclude it -- explaining both why `Delta_l_pair(a)` and the two-operator
`F_l(a)` fit gave wrong, resolution-independent, mutually-consistent-but-
incorrect answers (the contamination itself doesn't have to shrink
quickly with `n_refine` for it to already dominate a global-sum-based
estimator at any resolution tested), and why a windowed real-space fit
recovers the correct answer cleanly by simply avoiding that region.

**This does not (yet) explain why the l>=2 `Delta_l_pair` cluster is so
tightly self-consistent (~0.58-0.63, small spread) if it's "just"
short-distance contamination** -- a plausible follow-up hypothesis is
that the l_max=8-truncated Legendre series of a function with a
genuine near-`z=1` power-law feature converges slowly (Gibbs-phenomenon-
like), so the "contamination" seen at low-to-moderate l may itself have
a fairly universal, mesh/coupling-rule-independent SHAPE (explaining the
earlier campaign-wide observation that the anomaly was identical across
`naive`/`equal_area`/`equal_rp` and across `duality`/`exact_sinh`/
`owen_dual`) even though its physical origin is a lattice-scale (not
continuum) effect -- not yet tested directly.

**Practical conclusion for the campaign**: `Delta_sigma` should be
extracted via direct real-space fits of the two-point function
(`real_space_2pt_test`/`fit_real_space_2pt.py`), not via `F_l`/
`Delta_l_pair`/the two-operator Legendre fit -- those should now be
considered unreliable Delta estimators (though still potentially useful
for the ORIGINAL spherical-symmetry-verification objective, `R_l`/
`kappa2_r`, which is a different, off-diagonal-vs-diagonal question, not
a Delta-extraction question). **Next steps, not yet done**: (1) run
`real_space_2pt_test` across a proper resolution ladder and continuum-
extrapolate `Delta(a)` from the real-space fit, the way the campaign
previously tried (and failed) to do with `Delta_l_pair(a)`; (2) test the
"short-distance contamination leaks into every F_l" mechanism directly by
computing a *windowed* Legendre transform (excluding small-theta bins)
from the fine-grained real-space binned data and checking whether ITS
`Delta_l_pair` recursion behaves cleanly; (3) rerun with `equal_rp` mesh
mode to see whether the short-distance lattice artifact this real-space
method isolates is itself sensitive to the circumradius+perimeter
uniformity work from earlier today, now that there is a clean, windowed
way to probe short-distance behavior specifically instead of a global
`F_l` average.

## 2026-08-25 (cont): tested and REJECTED the sigma+epsilon hypothesis
(Delta_1=1/8, Delta_2=1) for the two-operator fit's second component

Per user suggestion -- the two lowest primaries of the exact 2D Ising
minimal model (c=1/2) are `1` (identity), `sigma` (Delta=1/8), `epsilon`
(Delta=1); no other primaries exist -- tested whether the free
`Delta_2~0.63-0.66` found by the two-operator fit is actually `epsilon`
in disguise, and whether fixing both dimensions to their EXACT known
values fits comparably to the free 4-parameter fit.

**Technical note first**: `F_l^cont(Delta) = int (2-2z)^{-Delta} P_l(z)
dz` genuinely diverges at `Delta=1` for every l, not just l=0 -- checked
via the recursion (`(l-1+Delta)/(l+1-Delta)` has a pole at `Delta=l+1`,
so `F_0` has the `Delta=1` pole and every higher `F_l` inherits it
multiplicatively). This is a real feature of 2D CFT, not a code bug: a
dimension-1 operator's contribution to a chordal-distance power law is
exactly log-divergent when integrated over the full sphere (`d^2x/|x|^2`
is the marginal case in 2D). Used `Delta_2=0.995` as a numerically-stable
proxy just below the pole for the "epsilon" tests below.

**Test 1: fix `Delta_1=1/8` exactly, let `Delta_2` float (3-param fit:
A, B, Delta_2).** Result across n_refine=16..128 (push7 naive):
`Delta_2` converges robustly to **0.656-0.660** at every resolution --
nowhere near 1 -- and chi2/dof degrades from <1 (the earlier fully-free
4-parameter fit) to **5-23** just from fixing `Delta_1` to the exact
sigma value. The data actively resists `Delta_1=1/8`.

**Test 2: fix BOTH `Delta_1=1/8` and `Delta_2=0.995` (epsilon proxy),
fit only the two amplitudes A,B (linear least squares).** Result:
**chi2/dof = 16,000-20,000** across n_refine=16..128 -- as catastrophic
as the original single-power-law fit. The sigma+epsilon combination is
decisively, unambiguously rejected by the data; the fitted `B` (epsilon
amplitude) comes out ~1000x smaller than `A`, unable to do anything
useful for the fit.

**Conclusion: this is NOT sigma+epsilon.** The freely-fit `Delta_2 ~
0.63-0.66` from the earlier two-operator fit does not correspond to
either exact primary of the 2D Ising minimal model, and there are no
other primaries in c=1/2 for it to be. Combined with the divergence-at-
Delta=1 finding above (epsilon structurally cannot contribute as a
simple second Legendre-power-law term to this decomposition in the first
place), this points toward `Delta_2` being a genuine **non-universal
lattice artifact** -- some discretization-specific effective operator or
mixing pattern with no continuum CFT meaning, rather than a second
physical operator dimension. Ad hoc scripts: `/tmp/fixed_d1_fit.py`,
`/tmp/sigma_epsilon_fixed.py` (session scratch, not committed -- the
result is negative/simple enough that committing a dedicated script
wasn't judged worthwhile, unlike the two-operator fit itself).

## 2026-08-25 (cont): windowed-Legendre-integral attempt to exclude the UV
FAILS (hard cut, smooth taper, excise-and-replace all give garbage for
l>=2) -- diagnosed as a confound in the reprojection method, not proof the
mechanism is wrong; led to a proper architecture fix (save_configs.cc)

Per user direction ("modify the integral to exclude the UV"), tried
literally truncating/tapering the domain of `int C(z) P_l(z) dz` (using
the fine, low-noise real-space binned data from the earlier n_refine=32
run) to exclude the near-coincident (`z->1`, small-theta) region. Three
variants tried: hard cutoff, smooth raised-cosine taper, and "excise and
replace" (swap the untrusted region for the already-fitted clean power
law instead of zeroing it, to avoid deleting the genuine large near-
coincident signal). **All three still give wildly unstable/negative F_l
for l>=2** (only l=1 improved, best in the excise-and-replace variant:
0.120, close to 1/8). Diagnosed the likely cause: re-Legendre-transforming
only 100 coarse, statistically-noisy real-space bins via naive
trapezoidal integration against `P_8(z)` (which oscillates 8 times across
`[-1,1]`) is numerically ill-conditioned regardless of windowing -- a
confound in *this specific reprojection method*, not proof that excluding
the UV doesn't work on the actual (much higher-precision) harmonic
measurement. The real `ising_s2_crit.cc` pipeline computes `S_lm` exactly
per configuration from ~10^4 site positions directly; re-deriving it from
a lossy 100-bin summary throws away exactly the precision needed for
l>=2.

**Fix: decouple data collection from analysis entirely**, per user
direction ("save the real-space correlators and do the spherical symmetry
tests by projecting afterwards"). New `src/save_configs.cc` (mirrors
`real_space_2pt_test.cc`'s mesh/coupling setup, drops all in-C++ analysis):
dumps site positions + `UpdateWeights` quadrature weights once
(`<run_id>_positions_<seed>.dat`), and every measured sample's raw spin
configuration, bit-packed via the existing `QfeIsing::WriteField` format,
appended sequentially to one binary file (`<run_id>_configs_<seed>.bin`,
no header -- `n_sites` from the positions file gives the per-config byte
stride `ceil(n_sites/8)`). At n_refine=32, 1 bit/spin means 10000 configs
is only ~12.8MB -- cheap to store, fully re-analyzable offline without
rerunning the (much more expensive) MC.

New `scripts/analyze_saved_configs.py`: loads the packed configs
(`numpy.unpackbits`), recomputes `Y_lm` at every site exactly
(`scipy.special.sph_harm`, matching boost's `spherical_harmonic(l,m,theta,
phi)` convention via `sph_harm(m,l,phi_azimuth,theta_polar)`), and forms
`S_lm`/`M_l[m,m']`/`F_l` via a single vectorized matmul (`spins @
(w*conj(Y))`) for the *standard* (unrestricted) case. **Validation: this
exactly reproduces the campaign's long-standing anomaly** (l=1
`Delta_l_pair` -> 0.333, l>=2 cluster -> 0.5-0.8, on the n_refine=16 smoke
dataset, 1000 configs) -- confirms the from-scratch Python
reimplementation is correct (matches the trusted C++ pipeline's math
exactly), not just superficially similar.

**UV-pair exclusion, done exactly (not via lossy reprojection) for the
first time**: `compute_Ml_near` precomputes, once, every site pair within
`theta_cut` via a `scipy.spatial.cKDTree` chord-distance query (fast, no
per-configuration cost), computes each such pair's config-averaged
`<s_i s_j>` directly from the raw saved spins, and forms
`M_l_near[m,m'] = sum_{(i,j) near} w_i w_j Y*_lm(i) Y_lm'(j) <s_i s_j>`
(both `i<j` and `j<i` orderings plus the `i=j` self-contact term, all
included) -- then `M_l_far = M_l_full - M_l_near` is the exact,
UV-pair-excluded harmonic decomposition, no numerical reprojection
involved anywhere. First result (n_refine=16 smoke dataset,
`theta_cut=0.3`, only 1000 configs -- noisy): l=1 improves from 0.333 to
**0.174** (real movement toward 1/8), but l>=2 is still not a flat
constant -- `Delta_l_pair` now runs 2.58 (l=3) monotonically DOWN to 0.56
(l=8), a qualitatively different (and possibly more encouraging) pattern
than the flat-wrong-plateau seen without UV exclusion, but not yet a
clean confirmation given the small dataset. **Scaling up now**: submitted
a full n_refine=32, 10000-config `save_configs` run (matching the
production stats used throughout this session) to re-run this exact
UV-exclusion analysis with much better statistics -- see next entry for
result once it completes.

## 2026-08-25 (cont): n_refine=32 UV-exclusion result reproduces the
n_refine=16 pattern closely -- SETTLES the question: exact pair-exclusion
does not fix the recursion, and now there is a clean mathematical
explanation why no windowing scheme ever could

`compute_Ml_near`/`compute_Ml_full` re-run on the full n_refine=32
(10242 sites, 10000 configs, `theta_cut=0.3`) `save_configs` dataset
(12m41s wall time, mostly the 81x81 harmonic-pair loop over 1,175,070
near-pairs -- much slower than the n_refine=16 timing extrapolated,
`compute_Ml_near`'s inner loop has some avoidable redundant work,
flagged but not optimized this session). Result: **l=1 -> 0.18069**
(vs. n_refine=16's 0.174 -- consistent), **l=2 -> -3.45102** (garbage,
same as n_refine=16's -3.86), **l=3..8 -> 2.616, 1.752, 1.453, 1.205,
0.954, 0.633** (a real, monotonically decreasing sequence, closely
matching n_refine=16's 2.579, 1.761, 1.442, 1.185, 0.963, 0.560). Two
independent resolutions giving the *same* qualitative pattern rules out
statistical noise as the explanation for why this doesn't cleanly recover
1/8 -- the effect is real and reproducible, just not a fix.

**Why, understood properly now**: `M_l_far = M_l_full - M_l_near` is the
exact harmonic decomposition of the correlator *masked to theta>=0.3 and
hard-zeroed below that* -- mathematically the *same class of object* as
the earlier failed hard-cutoff real-space-reprojection attempt (see the
"windowed-Legendre-integral attempt... FAILS" entry above), just computed
exactly from raw configs instead of approximately from 100 coarse bins.
**The earlier diagnosis (numerical noise/ill-conditioning from crude
reprojection) was only a partial explanation.** The real, more
fundamental reason is a genuine Gibbs-phenomenon-type fact about Legendre
transforms: a discontinuous, zero-padded real-space window smears power
across *every* Legendre coefficient, regardless of how exactly the
partial-domain integral is computed -- exact pair-restriction from raw
configs is subject to the identical mathematical limitation as lossy
reprojection, because the problem was never precision, it's that **the
standard CFT recursion `F_l^cont = [(l-1+Delta)/(l+1-Delta)] F_{l-1}^cont`
is only valid for full-sphere integration of an UNMODIFIED correlator** --
there is no way to windowed-subtract a piece out and expect the remainder
to still obey the same recursion, no matter how exactly the subtraction
is performed.

**This settles the "modify the integral to exclude the UV" investigation
conclusively: it cannot work, for a structural mathematical reason, not
an implementation deficiency.** No further attempt at windowing/tapering/
exact-pair-restricting the Legendre integral is worth pursuing -- three
independent methods (hard cutoff, smooth taper, exact pair-restriction at
two resolutions) all fail in the same qualitative way (l=1 improves
somewhat, l>=2 becomes unstable/non-recursion-obeying), consistent with
this being the expected mathematical consequence of real-space windowing
rather than a fixable bug. **The direct real-space fit
(`real_space_2pt_test.cc`/`fit_real_space_2pt.py`, giving the clean
Delta~0.128 result) remains the correct and now the *only* validated way
to extract Delta_sigma from this data** -- it works precisely because it
never routes through the full-sphere Legendre recursion at all, just
fits the measured correlator's shape directly over a trusted angular
window.

**Architecture note**: `src/save_configs.cc`/`scripts/analyze_saved_configs.py`
(raw-config save + offline harmonic reprojection) remain useful
infrastructure regardless of this specific negative result -- e.g. for
computing the *standard* (unrestricted) `M_l`/`F_l` at full precision
without needing a new C++ production run each time, or for future
symmetry-test (`R_l`/`kappa2_r`) work on saved ensembles. Not deprecated,
just not the tool that resolves the Delta-extraction question.

## 2026-08-25 (cont): submitted first production-scale real-space Delta(a)
ladder (job-array 7311120), built shard-level jackknife infrastructure
since real_space_2pt_test has none in C++

Per PLAN.md's "next step" after the real-space-fit breakthrough: the two
existing Delta=0.128/0.1287 results were single-seed login-node spot
checks (n_refine=16,32 only, no error propagation beyond one chain's naive
per-bin std-error). Needed a real resolution ladder with proper error
bars before claiming a continuum-limit result. Asked the user how error
bars should work; user chose full multi-shard jackknife over the cheaper
single-run-per-point option.

`real_space_2pt_test.cc` has no in-C++ jack_block_size accumulator (unlike
`ising_s2_crit.cc`) -- rather than add one, used each independent-seed
shard itself as the jackknife unit (delete-one-shard resampling, done
offline in Python). This is arguably cleaner than same-chain
`jack_block_size` blocks since shards share zero thermalization history,
and required no C++ changes.

New `cluster/sge/real_space_shard_task.sh` (mirrors
`production_shard_task.sh`'s task-index layout: `shard_idx = idx %
N_SHARDS`, `l_idx = idx / N_SHARDS`; naive mesh_mode only, no mesh-cache
prewarm needed) and `cluster/sge/submit_real_space_ladder.sh`. Config
matches the validated spot checks exactly: `exact_sinh` coupling, naive
mesh, `n_therm=2000 n_traj=20000 n_skip=2 n_wolff=5 n_metropolis=4
n_bins=100 n_ref=300`. Ladder `n_refine={8,16,24,32,48,64,96,128}` (8
points, brackets the two already-validated points with room for a real
`1/n_refine^p` extrapolation; not extended down to n_refine=2-4 like
push7 -- too few sites for the `theta_min=0.3` window to leave enough
bins). `N_SHARDS=24`/point, 192 tasks total. `SEED_BASE=800000` (disjoint
from every prior push).

New `scripts/fit_real_space_ladder.py`: loads all shards per ladder point,
delete-one-shard jackknife -> `Delta_jk(n_refine)` with jackknife error,
then fits `Delta(a) = Delta_inf + c*n_refine^-p` (p=2 default, overridable)
across the ladder weighted by jackknife errors to get `Delta_inf` and its
error, compared directly against the exact 1/8.

Verified the whole pipeline end-to-end before submitting: ran
`real_space_shard_task.sh` locally for 4 fake tasks (tiny n_refine=4,8
ladder, 2 shards, n_traj=200) confirming task-index math and output paths
match, then ran `fit_real_space_ladder.py` against that toy output and
confirmed it loads/jackknifes/fits correctly (2-point ladder correctly
returns `nan` chi2/dof since dof=0 there, expected, not a bug). **Job
submitted**: SGE job-array 7311120, 192 tasks, `h_rt=2:00:00`, `-P qfe`.
Submitted `qw` (queued, waiting for slots) rather than immediately
running -- cluster was busy with a large unrelated array job (7309036,
`gb_array_o`, not this campaign). Check with `qstat -u misra | grep
s2prec_real_space_ladder` / `qacct -j 7311120`. **Next session**: once
complete, run `scripts/fit_real_space_ladder.py
campaign_runs/real_space_ladder_2026-08-25/` and read off `Delta_inf`.

**Note on this run's SGE behavior**: `qstat -u misra | grep s2prec_real_space`
intermittently returned zero matches for several minutes while ~90 of the
192 tasks were still genuinely `r` (running) -- confirmed via `qstat -j
7311120` (still showed active `usage` lines/cpu time) and via output file
counts still growing. Two separate false "job is done" reads happened this
session as a result. **Do not trust a single `qstat -u $USER | grep
<name>` snapshot as proof a job has finished** -- corroborate with `qstat
-j <jobid>` (errors only once the job is truly gone) and/or output file
counts before declaring an array complete.

## 2026-08-25 (cont): icosahedral-orbit-averaging test, per user idea --
confirms large (4-5x at a single fixed pair, ~2x when already pooling a
handful of pairs) free variance reduction, no bias

User proposed averaging the real-space two-point function over the exact
120-element icosahedral point group orbit (`QfeLatticeS2::G`, loaded from
`grp/elem/o3q5.dat` -- already computed internally by every `QfeLatticeS2`
via `ReadSymmetryData`, just not previously used by any real-space
measurement code) instead of `real_space_2pt_test.cc`'s z=cos(theta)
binning. Two distinct motivations, only the first tested so far: (1) exact
zero-extra-simulation-cost variance reduction -- for a fixed pair of sites
and a single measured configuration, apply every g in `G` to both site
positions, look up the resulting site indices, and average `s_i*s_j` over
all group images; since the `exact_sinh` action is exactly icosahedrally
invariant, every image has the identical ensemble-average expectation
value, so this is a rigorous, bias-free way to shrink the per-config
variance for free. (2) whether the *existing* z-binning method (which
pools pairs from different, icosahedrally-inequivalent orbits as long as
they share the same z) is itself a source of the l>=3 `kappa2_r`-style
contamination -- not yet tested, a different and harder question (concerns
different (ref,targ) pairs at the same z, not one pair's own orbit).

**New small standalone test, `src/test_orbit_avg.cc`** (built manually,
not wired into `cluster/build.sh` -- a feasibility/value probe, not
production code; brute-force O(n_sites) position->site lookup per group
element, fine at n_refine<=32, would need a KD-tree to scale further).
Reuses `QfeLatticeS2::G`/`r[]` (already loaded, no new group-theory code
needed -- `GrpElemO3::operator*(Vec3)` already existed) and the
`exact_sinh` coupling loop copied verbatim from `real_space_2pt_test.cc`.
For each of several reference sites, finds a target ~90 degrees away, then
computes the orbit of the (ref,targ) PAIR under all 120 elements of `G`
(apply g to both positions, brute-force-match to site indices, dedupe --
every ref/targ pair tested had a full, non-degenerate 120-element orbit,
i.e. trivial pair-stabilizer for all sites tried). Per measured config,
tracks both a "single-pair" estimator (raw `s_ref*s_targ`, or the average
of a handful of independent pairs) and an "orbit-averaged" estimator (same
config, same underlying pairs, averaged over their full 120-element
orbit).

**Result 1 (n_ref_sites=5-8, i.e. baseline already pools several
independent pairs)**: n_refine=8, 2500 configs: error ratio
single/orbit = **2.02x**. n_refine=32, 2500 configs: **2.35x**. Both
means agree with the un-averaged estimator within ~2 sigma (-0.20, -1.85
sigma) -- consistent with orbit-averaging being an unbiased variance
reduction, not a shift.

**Result 2, per user follow-up ("define the two point function from the
north pole and do the 120 averaging, see a larger error reduction") --
using a single fixed reference/target pair with no baseline pooling at
all (`n_ref_sites=1`) shows a much larger, cleaner improvement**:
n_refine=16, 4000 configs: error ratio **4.34x**. n_refine=32, 4000
configs: **5.24x**. Means still agree with the raw single-pair estimator
(-0.22, +0.84 sigma). The smaller ~2x ratio in Result 1 was diluted by the
baseline itself already averaging over multiple quasi-independent pairs;
the true value of orbit-averaging for a single measured pair is closer to
4-5x free variance reduction (still well below the naive
sqrt(120)~11x expected for fully independent images -- orbit images share
real correlations through the same global critical configuration, as
expected physically).

**Not yet done**: this only tests hypothesis (1) above (variance
reduction for a single fixed pair's own orbit) -- it does NOT yet test
whether z-binning's pooling of *different* pairs at the same z is a source
of the l>=3 symmetry-breaking-style contamination in the existing
real-space Delta fit (hypothesis (2)); that would require redoing
`real_space_2pt_test.cc`'s bin construction using exact pair-orbit
membership instead of z alone, a real driver change, not yet built.
**Practical implication for the running production ladder
(`campaign_runs/real_space_ladder_2026-08-25/`, job 7311120)**: that job
uses the existing (non-orbit) z-binned method and is still a valid
measurement of `Delta(a)` -- this orbit-averaging result is a candidate
efficiency/precision upgrade for a *future* push (could cut the
shard/traj budget needed for a given error target by ~4-5x if built into
`real_space_2pt_test.cc` properly with a KD-tree for production
n_refine), not a reason to distrust or redo the ladder currently running.

## 2026-08-25 (cont): built --orbit_avg into real_space_2pt_test.cc
(user direction: "build it in next"), found a real short-distance
disagreement vs. the z-binned method -- confirms hypothesis (2), but
outside the existing fit window the two methods agree

Wired icosahedral-orbit averaging (validated above via
`test_orbit_avg.cc`) directly into the production driver behind a new
`--orbit_avg` flag, default off (existing behavior, and the currently-
running `real_space_ladder_2026-08-25` job, untouched).

**New `SpatialHash` class** (grid-bucketed `unordered_map<int64_t,
vector<int>>` over quantized unit-sphere coordinates, cell size = 2x mean
edge length, 3x3x3-neighborhood nearest-site query) replaces
`test_orbit_avg.cc`'s brute-force O(n_sites) lookup -- O(1) per query,
needed to scale past n_refine~32.

**Orbit-avg precompute, done in the same O(n_ref * n_sites) pass the
non-orbit path already pays for** (no added asymptotic cost): for each
(ref, bin), track the single target site landing in that bin closest to
its z_center (replaces the non-orbit path's "pool every site in the bin"
with "one representative pair per bin"), then compute that one pair's
exact `|G|=120` orbit once via the spatial hash (dedup, since some
pairs have a nontrivial stabilizer -- smoke test found mean orbit size
118.2/120 at n_refine=16, n_ref=50). Per-config measurement then averages
`s_i*s_j` over each (ref,bin)'s precomputed orbit-pair list, then over
`ref`, mirroring the non-orbit path's two-stage pooling exactly so the
output file format (and every downstream script) is unchanged.

**Verification**: (1) default (non-`--orbit_avg`) path re-tested byte-
for-byte unchanged in structure (only touched code inside the new `else`
branch); (2) `--orbit_avg` smoke run at n_refine=16 completes cleanly,
reports `mean orbit size = 118.20 (|G|=120)`; (3) **direct comparison on
the identical MC trajectory** (same seed for both, so both analyze the
exact same sequence of spin configurations, only the offline binning
differs) at n_refine=16, n_ref=50, n_bins=20, 1250 measured configs: bins
0-18 (theta=0.55 to 2.82, i.e. everywhere except the single shortest-
distance bin) agree between the two methods to within 0.07-0.79 sigma,
with a small, monotonically growing (toward short theta) systematic
offset (orbit consistently a hair below default) that never exceeds 1
sigma in this range. **Bin 19 (theta=0.318, the single shortest-distance
bin, right at the edge of the campaign's existing fit-window floor
`theta_min=0.3`) disagrees at 6.04 sigma** (default=0.4795 vs
orbit=0.4574) -- direct evidence that the non-orbit z-binning method's
pooling of icosahedrally-inequivalent pairs really does introduce a
resolvable bias, concentrated exactly at short distances. **This
confirms hypothesis (2)** from the earlier orbit-averaging entry (z-
binning contamination is real, not just theoretical) but also shows it is
**already excluded by the existing `theta_min=0.3` fit window** -- the
currently-running `real_space_ladder_2026-08-25` job's Delta(a) results
should not need to be redone on this basis; this is a real, useful
independent confirmation of why that exclusion window is doing real work,
not (yet) a reason the ladder's `Delta_inf` estimate itself is wrong.

**Not yet done**: no production run with `--orbit_avg` yet (needs a fresh
push, own manifest/seed range — `--orbit_avg` was only smoke-tested at
n_refine=16/50 shards so far); the smaller-theta systematic trend seen in
bins 13-18 (still <1 sigma each but consistently one-directional) is
worth re-checking at full production statistics before concluding it's
negligible outside bin 19 specifically.

**Update (same session, immediately after)**: `cluster/sge/
real_space_shard_task.sh` gained an `ORBIT_AVG=1` passthrough and
`cluster/sge/submit_real_space_ladder_orbit.sh` (new submit script,
matched ladder/coupling/stats to `submit_real_space_ladder.sh`,
`SEED_BASE=810000` disjoint) was built and submitted as a second
production array (job 7311202) for direct comparison against the
z-binned job (7311120), both running concurrently. `CLAUDE.md` also
gained a full new "Real-space Delta_sigma extraction" section (previously
undocumented there despite being the campaign's primary Delta-extraction
method since the BREAKTHROUGH) plus a pivot notice pointing future
sessions to it first.

## 2026-08-25 (cont): confirmed F_0/F_1 CAN recover Delta correctly, with
the correct (matched) windowed integration scheme -- resolves why the
earlier windowed-Legendre attempt "failed"

User question: "are you able to get delta from F0 and F1 with the correct
integration scheme?" Re-examined why the earlier windowed-Legendre-
transform attempt (hard cutoff/taper/exact pair-restriction, all "FAILS"/
"SETTLES the question" entries above) concluded windowing doesn't work:
that attempt computed a windowed `F_l` from data, then tried to invert it
via the **unwindowed** closed-form recursion (`Owen_section_D`'s
`(l-1+Delta)/(l+1-Delta)` recursion and the shortcut identity
`Delta=2*F_1/(F_1+F_0)`) -- both are only valid for a full-sphere
integral of the pure power law. That is a scheme *mismatch* (windowed
data vs. unwindowed theory), not evidence that F_0/F_1 structurally
cannot carry `Delta` once windowed.

**New `scripts/windowed_f0_f1_delta.py`**: computes `F_0`, `F_1` from
data via a direct sum over included bins (`real_space_2pt_test` bins are
uniform in `z=cos(theta)`, so `sin(theta)dtheta = -dz` exactly and no
extra Jacobian is needed -- `F_l_data = sum_bins mean(bin)*P_l(z_center)*dz`),
and separately computes `F_0^model(Delta)`, `F_1^model(Delta)` by direct
numerical integration (`scipy.integrate.quad`, no recursion at all) of
`(2-2z)^{-Delta}*P_l(z)` over the **same** `[z_min,z_max]` window. Solves
for `Delta` by matching the amplitude-independent ratio
`2*F_1/(F_1+F_0)` between data and model (`scipy.optimize.brentq`).

**Result: this works.** n_refine=32 (`campaign_runs/real_space_2pt_2026-08-25/`),
`theta_min=0.3, theta_max=2.8`: **Delta=0.137** (vs. the naive full-sphere
`Delta_s(a)`'s ~0.33). Across 18 window choices (`theta_min` 0.2-0.8,
`theta_max` 2.6-2.9): **Delta ranges 0.104-0.140**, centered near the
exact 0.125. n_refine=16 smoke dataset: 0.119-0.154 across 3 windows.
Noisier/more window-sensitive than the full real-space bin fit
(`fit_real_space_2pt.py`'s 0.1275-0.1284, tight) -- expected, since this
uses only 2 numbers (F_0,F_1) per measurement instead of ~60-100
individual bins, so it has much less statistical leverage, but it is
**not** systematically biased the way the unwindowed full-sphere method
is. **Conclusion**: F_0/F_1 are not inherently incapable of encoding
`Delta` correctly -- the earlier "windowing doesn't work" finding was
specifically about the RECURSION (`F_l^cont(Delta)` for `l>=2` and the
shortcut identity), not about F_0/F_1 as quantities. This is a useful
independent cross-check (different quadrature/moment-based method,
same conclusion as the direct real-space fit) but not a replacement for
`fit_real_space_2pt.py` as the primary method (less statistical power,
single-run/no-jackknife so far -- has not been run against a full
production shard ladder).

## 2026-08-25 (cont): reconstructed push7's real-space correlators from
F_l(a), then applied the windowed-F0/F1 method to the RESUMMED curve --
only partial improvement (0.33 -> 0.27), confirming the l_max=8
truncation smears contamination across the whole curve rather than
localizing it

Per user requests ("reconstruct the real space correlators for the push7
data", then "get the delta_s(a) from F0 F1 from the corrected integration
method from the real space correlators"): built
`scripts/reconstruct_real_space_from_Fl.py` (resums
`C_resum(theta) = sum_{l=0}^{8} (2l+1)/2 * F_l(a) * P_l(cos(theta))` from
push7's pooled jackblocks, across the full ladder, normalized to 1.0 at a
reference angle for shape-only comparison against the two validated
direct real-space spot checks at n_refine=16,32 -- reproduces the
already-known shape mismatch visually) and
`scripts/delta_s_windowed_from_resummed.py` (applies
`windowed_f0_f1_delta.py`'s matched windowed-integration method directly
to the resummed curve instead of directly-measured bin data, across
push7's full 13-point ladder).

**Result: windowing the resummed curve only partially helps.**
`Delta_s(a)` (unwindowed) plateaus at 0.33-0.41 across the ladder (matches
`run_cft_ladder.py`'s earlier table exactly, sanity check passed); the
SAME `theta_min=0.3,theta_max=2.8` window applied to the resummed curve
gives **0.27-0.29** at every ladder point -- real, reproducible movement
toward 0.125, but nowhere near the ~0.10-0.15 the identical window
achieves on the *directly measured* real-space correlator
(`windowed_f0_f1_delta.py`'s result, earlier this session).

**Interpretation**: this is new, independent confirmation of a mechanism
already suspected (see the "two-operator fit"/"REJECTED" entries'
Gibbs-phenomenon discussion) -- a low-l_max (`l_max=8`) truncated
Legendre series cannot *localize* a short-distance feature the way a
finely-binned real-space measurement can. Every `P_l` oscillates across
the *entire* `[0,pi]` domain, so truncating the series at `l=8` smears
whatever short-distance contamination exists in the true correlator
across the whole reconstructed curve, not just near `theta~0` --
windowing out the small-`theta` region of the *resummed* curve therefore
only partially removes it (the smeared remainder still contaminates the
window), whereas windowing the *directly measured* correlator (100 fine
bins, no harmonic truncation at all) removes the actual localized
short-distance defect cleanly. **This is the clearest evidence yet for
why the real-space method (measure directly, then window) is fundamentally
different from, and superior to, any resum-then-window approach built on
top of a low-l_max harmonic measurement** -- not just a different
implementation of the same idea. Two plots sent to user: shape-normalized
resummed-vs-direct correlator overlay, and `Delta_s(a)` vs `1/n_refine`
(unwindowed vs. windowed-resummed) across the full ladder.

## CLEAN CONFIRMATION: Delta_sigma = 1/8 from the real-space continuum
extrapolation, both z-binned and orbit-averaged (2026-08-25/26)

Both production ladders submitted earlier this session finished cleanly:
job 7311120 (z-binned, `real_space_ladder_2026-08-25/`, 192/192 tasks,
`qacct` confirms `exit_status=0` on every task) and job 7311202
(orbit-averaged, `real_space_ladder_orbit_2026-08-25/`, 192/192 tasks, same
clean exit). Ran `scripts/fit_real_space_ladder.py` against both.

At the **previously-assumed `p=2`** (`Delta(a) = Delta_inf + c/n_refine^2`)
extrapolation ansatz, neither ladder cleanly confirms 1/8: z-binned gives
`Delta_inf = 0.12628(29)`, chi2/dof=3.06, **4.36 sigma** away from 0.125;
orbit-avg gives `0.12605(23)`, chi2/dof=1.60, **4.47 sigma** away. This
looked like a real, if small, discrepancy.

**The fix was the extrapolation power, not the data.** Re-ran with
`--power 1` (i.e. assume a *linear*, not quadratic, leading discretization
correction `Delta(a) = Delta_inf + c*n_refine^-1`) and both ladders snap
into excellent agreement with the exact CFT value:

| ladder | power | Delta_inf | chi2/dof | deviation from 1/8 |
|---|---|---|---|---|
| z-binned (7311120) | 1 | 0.12497(20) | 0.82 | -0.14 sigma |
| orbit-avg (7311202) | 1 | 0.12515(23) | 0.97 | +0.65 sigma |
| z-binned (7311120) | 2 | 0.12628(29) | 3.06 | +4.36 sigma |
| orbit-avg (7311202) | 2 | 0.12605(23) | 1.60 | +4.47 sigma |

`p=1` is decisively preferred on chi2/dof alone (0.82-0.97 vs 1.60-3.06)
independent of any prior on the exact answer, and the two independent
binning methods (z-binned vs. orbit-averaged, disjoint seed ranges,
same ladder/coupling/stats) agree with each other at every `n_refine`
point and both extrapolate to 1/8 within <1 sigma at `p=1`. **This is a
clean, cross-validated confirmation that `Delta_sigma = 1/8` in the
continuum limit via the real-space method** — the strongest positive
result the campaign has produced. A linear-in-`a` leading correction is
also physically sensible: this simplicial mesh/coupling construction has
no symmetry (like a hypercubic lattice's reflection symmetry) that would
forbid an O(a) term, so O(a) rather than O(a^2) is the generic expectation,
not a surprise.

New `scripts/plot_real_space_ladder.py` (reuses `fit_real_space_ladder.
analyze_point` directly, no duplicated fit logic) makes the `Delta_sigma(a)
vs 1/n_refine` overlay for both ladders + `p=1` fit lines + the exact-1/8
reference line: `campaign_runs/real_space_ladder_2026-08-25/plots/
delta_vs_invL_power1.png`. `fit_real_space_ladder.py`'s `--power` flag
already existed (default was `2`, now known to be the wrong default for
this data) — did not change the script's default since `--power 2` was
never validated as "the right answer with the wrong extrapolation," it's
simply a flag that needs to be set correctly per-call; documented the
`--power 1` finding in `CLAUDE.md` instead so future sessions don't
silently rerun with the stale default and mistake a bad ansatz for a real
discrepancy.

**Next steps for the campaign**: (1) this closes objective (2)'s central
question (Delta_sigma extraction) for the sigma operator — the campaign
could consider running the same real-space method for the epsilon
operator (Delta=1) next, not yet attempted. (2) objective (1) (spherical-
symmetry verification via `R_l`/`kappa2_r`) is still open per the
2026-08-24 gate-fail finding (l>=3 plateau) — unaffected by this real-space
work, still needs the `equal_rp` mesh fix or another explanation. (3)
repo still not committed to git (deferred by user 2026-08-24) — consider
doing this now that there's a clean positive result worth preserving.

## real_space_2pt_test.cc rewritten: dropped z=cos(theta) binning for exact
icosahedral-symmetry equivalence classes (2026-08-26)

User feedback, in direct response to being shown the driver's output
format: **"dont do binning like this"**, followed by two clarifications
when asked what should replace it — **"save the whole two point function
averaged over the 120 orbit"** and **"do this in an accumulated stats way
where the mean and variance of step n can be updated to step n+1"**.
Confirmed via `AskUserQuestion` that the intended replacement was "no
binning — per-pair/orbit-exact treatment" (as opposed to equal-angle bins
or adaptive bins).

**What changed**: both the original z=cos(theta) pooled-bin path and the
`--orbit_avg` bin-representative path it later grew into (built
2026-08-25, see "built --orbit_avg into real_space_2pt_test.cc" entry
above) are removed. Binning by z pools pairs that are only *approximately*
equivalent — exactly equivalent only in the continuum SO(3) limit, which
is precisely what the campaign is trying to test/extract from, not assume
— and the earlier same-seed z-binned-vs-orbit-avg comparison already found
a real 6-sigma disagreement at the shortest-distance bin from exactly this
effect (masked only because that bin was already excluded by the fit
window).

**New method**: for each of `--n_ref` random reference sites `i`, every
other site `j` is assigned to an exact equivalence class under `Stab(i)`
(the subgroup of the 120-element icosahedral point group `G` that fixes
site `i` exactly, found via the existing `SpatialHash` nearest-site
lookup) — two targets are pooled together if and only if some element of
`Stab(i)` maps one onto the other exactly, which is the complete set of
pairs the lattice's own discrete symmetry *guarantees* have identical
`<s_i s_j>`, no continuum assumption involved. For each resulting class, a
canonical pair's full 120-element `G`-orbit (both `i` and `j` moving
together, not just `Stab(i)`) is precomputed once; per configuration, this
driver now saves the two-point function averaged over that entire orbit —
directly satisfying "save the whole two point function averaged over the
120 orbit". Each class's per-config orbit-average feeds a `QfeMeasReal`
(`include/statistics.h`) exactly the way every other observable in this
campaign is measured — an online accumulator holding only running
sum/sum2/n, updated one configuration at a time, never buffering raw
per-config samples — which is exactly the "accumulated stats, mean/
variance of step n updated to step n+1" the user asked for; no new
accumulator machinery was needed, `QfeMeasReal` already worked this way.

Output format changed: `<data_dir>/<run_id>_real_space_2pt_<seed_hex>.dat`
now has one row per (reference site, class) with columns `ref_site
target_site theta mean err n orbit_size` — no more `bin z_lo z_hi
z_center` columns, since there is no fixed bin grid any more (`theta` is
the pair's own exact angular separation). **This is a breaking format
change** — old .dat files from before this rewrite are not readable by the
updated analysis scripts.

Added `--theta_min`/`--theta_max` (default: full sphere, `[0, pi]`) purely
for compute-cost control, not as a re-introduction of binning: it only
skips *building* classes outside the window, never pools pairs that
wouldn't otherwise be pooled. This exists because building a class costs
`O(|G|)=120` hash lookups, and a *generic* (non-high-symmetry) reference
site has trivial `Stab(i)` — so class count is close to `n_sites`, and at
large `n_refine` (e.g. 128, `n_sites~1.6e5`) with `n_ref~300` the
full-sphere setup phase alone is `O(1e10)` hash lookups, likely
impractical. Narrowing the window (e.g. to the campaign's usual `[0.3,
2.8]` fit range) is the practical mitigation, documented in the file's
header comment.

Compiled clean (`g++ -g -O3 -fopenmp ... src/real_space_2pt_test.cc`, same
flags as `cluster/build.sh` uses for `ising_s2_crit.cc`; this driver isn't
in `build.sh` itself, ad hoc compile as before). Smoke-tested end to end
at `n_refine=4` (login node, `--theta_min 0.3 --theta_max 2.8`, `--n_ref
20`): setup found 1333 classes with mean orbit size 109/120 (some
coincidental-image dedup, as expected on a coarse mesh), each class's
`n=1000` matches the expected measured-configuration count exactly,
`fit_real_space_2pt.py` ran against the output without error
(`Delta=0.147` at this tiny, untrustworthy resolution — not a validation,
just a pipeline smoke test). A two-shard ladder smoke test (single
`n_refine=4` point) also ran `fit_real_space_ladder.py` end to end
successfully.

**Downstream scripts updated to match**:
- `scripts/fit_real_space_2pt.py` — new column indices, and now weights
  the fit by each class's own `err` (previously unweighted OLS over bin
  means, since old bins had no natural per-row weight built into the fit).
- `scripts/fit_real_space_ladder.py` — the delete-one-shard jackknife
  previously worked by averaging bin means index-for-index across shards,
  which relied on every shard sharing a common bin grid (true for fixed
  z-bins). **That assumption is now false**: each shard draws its own
  random reference sites (seed-dependent), so shards no longer share a
  common theta grid at all. Replaced with: each jackknife fold does one
  `err`-weighted nonlinear least-squares fit over the *concatenation* of
  every included shard's (theta, mean, err) rows — pooling, not averaging,
  since different shards' classes are mostly different physical pairs, not
  repeated measurements of the same pair. `scripts/plot_real_space_ladder.py`
  needed no changes (calls `analyze_point` opaquely).
- `cluster/sge/real_space_shard_task.sh` / `submit_real_space_ladder.sh` —
  dropped `N_BINS`/`ORBIT_AVG` env vars and CLI flags (would otherwise
  silently pass a now-nonexistent `--n_bins`/`--orbit_avg` flag to the
  rebuilt binary), added `THETA_MIN`/`THETA_MAX` passthrough (defaulted to
  the existing validated fit window `0.3`/`2.8` in the submit script).
  `submit_real_space_ladder_orbit.sh` marked obsolete in a header comment
  (no more separate orbit-avg toggle to compare against) but left in place
  for historical reference to job 7311202.
- `CLAUDE.md`'s "Real-space Delta_sigma extraction" section rewritten to
  match (`--n_bins`/`--orbit_avg` bullets replaced).

**Not yet done**: re-run the two validated spot checks (n_refine=16,32,
`Delta=0.128`/`0.1287(3)`) against the rewritten driver to confirm the new
method reproduces them before trusting a new production ladder; the two
existing production ladders (jobs 7311120/7311202, the CLEAN CONFIRMATION
result) predate this rewrite and used the old binned/`--orbit_avg` format
— they remain valid as-is (nothing about their own data changed), but any
*new* ladder must use the rewritten driver and updated scripts above, not
mix formats.

## Offline M_l[m,m'] construction from MC-averaged two-point functions:
full_corr_test.cc, then storage-scalable configs-based alternative
(2026-08-26, same day)

**User direction**: "we will build the [spherical-symmetry] matrices from
the MC averaged two point functions. Building them in the simulator was a
mistake since the two point functions were not smoothed by high stats
yet" -- i.e. `M_l[m,m']` should be built OFFLINE from an already-converged
correlator, not accumulated per-raw-configuration inside
`ising_s2_crit.cc --l_max` the way it always has been. Follow-up: "the 120
averaging can always be done for free statistics" (apply the exact
icosahedral-orbit averaging from the `real_space_2pt_test.cc` rewrite here
too, since it's free variance reduction once every pair is genuinely
measured).

**First implementation: `src/full_corr_test.cc`** (new `bin/full_corr_test`).
Accumulates the FULL `n_sites x n_sites` `<s_i s_j>` matrix via a
symmetric rank-1 update per measured configuration
(`Eigen::selfadjointView<Upper>().rankUpdate`, row-major storage so
row-slices are contiguous for the output `fwrite`), then a post-processing
`SymmetrizeOverOrbits` pass replaces every upper-triangle entry with the
mean of its exact 120-element `G`-orbit (read-all-then-write-all per
orbit, safe because orbit membership is symmetric — an unvisited pair
implies its whole orbit is unvisited). Smoke-tested at n_refine=4: exact
diagonal=1.0, exact matrix symmetry, and three geometrically distinct
same-theta pairs collapsed to bit-identical values post-symmetrization —
confirms the free-orbit-averaging pass works. **Storage is dense
`O(n_sites^2)`** (~840MB at n_refine=32, ~215GB at n_refine=128) — flagged
in the header comment as impractical much past n_refine~32-48.

**`scripts/build_ylm_matrix_from_corr.py`**: contracts a `full_corr_test`
shard's correlator with `Y_lm` quadrature weights (`M_l = A^H @ corr @ A`,
`A[i,k]=w_i*conj(Y_lm(i))`, matching `CalcYlm`'s
`boost::math::spherical_harmonic(l,m,polar,azimuth)` convention via
`scipy.special.sph_harm(m,l,phi,theta)` — same mapping already used in
`analyze_saved_configs.py`). Turns each shard into exactly the `blocks`
structure `scripts/symmetry_test.py`'s `analyze_blocks()` already expects
(`blk['n']=n_meas`, `blk['sum'][(l,m,mp)] = n_meas*M_l[m,mp]` for `m<=mp`)
so R_l/chi2/kappa{2,3,4}_r and their delete-one-shard-jackknife errors are
computed by the existing, already-verified code — nothing new
reimplemented. Smoke-tested on 3 shards at n_refine=4, l_max=4: l=1,2 come
back consistent with zero (icosahedral group has no invariant below l=6,
matching the whole campaign's prior findings), l=3 shows a real non-zero
`kappa2_r` (7e-4 +/- 1.6e-4) — the same qualitative anomaly this campaign
already found via the in-simulator `--l_max` path at coarse resolution.
Good cross-check that the offline construction is correct.

**Then asked for options to fix the O(n_sites^2) storage problem before
scaling this up.** Presented: (1) skip the dense matrix entirely, reuse
already-existing `save_configs.cc` (bit-packed raw configs,
`O(n_meas*n_sites)` storage) + offline Y_lm projection — mathematically
IDENTICAL estimator to (a), by linearity of expectation
(`<S_lm S*_lm'> = sum_ij w_i w_j Y*_lm(x_i)Y_lm'(x_j) <s_i s_j>` regardless
of whether the config-average happens before or after the Y_lm
contraction) — so nothing is lost by choosing this over the dense-matrix
path, only computational order changes; (2) only store one representative
row per icosahedral site-orbit-type (~120x storage reduction, still
O(n_sites^2/|G|) asymptotically); (3) just cap `full_corr_test.cc`'s
ladder at small `n_refine` and lean on the real-space method (unaffected)
for the rest. **User chose option 1.**

**`scripts/build_ylm_matrix_from_configs.py`** (new): reuses
`analyze_saved_configs.py`'s existing `load_positions`/`load_configs`/
`build_ylm_matrix`/`compute_Ml_full` (already built 2026-08-25 for the
now-dead-end windowed-Legendre investigation, but the raw-config-loading +
`M_l` machinery itself was always correct and reusable) directly rather
than duplicating it, and wires shards into the same `symmetry_test.
analyze_blocks` pipeline as `build_ylm_matrix_from_corr.py` above — same
output table format, same delete-one-shard jackknife. Rebuilt
`bin/save_configs` (still compiled clean against current headers).
Smoke-tested end to end: 3 shards, n_refine=4, `--n_traj 4000 --n_skip 2`
-> 2000 configs/shard, 42000-byte `.bin` files (`ceil(162/8)*2000` bytes,
exactly matches the packed-bit format), and the resulting R_l/kappa2_r
table shows the same qualitative pattern (l=1,2 clean, l=3 elevated) as
the `full_corr_test.cc` path's smoke test (different seeds, so not a
bit-for-bit match, but the same physics signature).

**This is now the recommended path for the spherical-symmetry campaign
going forward, especially at large n_refine** where `full_corr_test.cc`'s
dense storage would be infeasible — e.g. n_refine=128, n_meas=10000:
`save_configs` needs ~205MB per shard vs. `full_corr_test`'s ~215GB dense
matrix at the same resolution. `full_corr_test.cc` remains built and
validated (kept for reference / cross-checks at small-to-moderate
n_refine where the dense-matrix path is convenient) but is not the
primary route for a real production ladder.

**Not yet done**: no production run yet with either offline-M_l path —
both are smoke-tested only (n_refine=4, tiny stats). Next step: pick a
resolution ladder, decide `--l_max`, and run a real production push
(`save_configs.cc` + `build_ylm_matrix_from_configs.py`) to actually
re-attempt the Phase 1 symmetry gate that stalled on the l>=3 plateau
(2026-08-24 finding) — this offline path doesn't change the *physics*
being tested, only how the same `M_l` estimator gets computed, so it is
not expected to resolve that anomaly on its own; it's primarily useful
for reaching larger `n_refine` than the in-simulator path's
`equal_rp`-focused pushes did.

## North-pole-centered two-point function, offline from configs; production
save_configs ladder launched and completed (2026-08-26, same day)

**User pushback, twice, both correct**: first, "why are you saving
configurations, just save the expectations of the two point function and
its variance" -- clarified that computing/storing the full two-point
function online (`full_corr_test.cc`'s per-config rank-1 update) pays real
`O(n_sites^2)` compute for no benefit over deferring to a single offline
matmul from cheaply-saved raw configs. Second (after being told even the
offline dense-matrix path is memory-bound): "how could the MC average be
more expensive than saving configurations, that's crazy" -- correct;
writing bit-packed spins is `O(n_sites)` I/O per config and should never
be the expensive part, only *materializing the full pairwise matrix*
(`O(n_sites^2)` storage) is unavoidably expensive, independent of timing.

**Resolution, per explicit final direction ("I only want the two point
function of the spins" / "all I want is the MC expectation value for the
two point function and its variance ... all of our checks can be done on
those data")**: don't materialize the full matrix at all for the
two-point-function deliverable. New `src/analyze_north_pole_configs.cc`
computes `<s_i s_j>` for exactly ONE reference row -- the north pole
(found robustly via `argmax(z)`, confirmed exactly `[0,0,1]`) -- offline
from `save_configs` shards. Reuses the same `Stab(i)`-dedup +
full-120-element-`G`-orbit-averaging construction as
`real_space_2pt_test.cc`/`full_corr_test.cc`, applied to one row only:
`O(n_sites*|G|)` setup, `O(n_meas*n_classes)` measurement, no
`O(n_sites^2)` anywhere. Shards combined via Welford's algorithm (mean +
variance across independent shards). Smoke-tested at n_refine=4:
`|Stab(north pole)|=10` (matches the C5v vertex-symmetry expectation),
161 targets collapse to 24 classes, and cross-checked directly against
`compute_two_point_stats.py`'s dense-matrix row 0: same physical values,
~5.6x smaller error from the orbit averaging (matches `test_orbit_avg.cc`'s
earlier ~4-5x finding).

**Production ladder launched and completed same session**: locked in
ladder `n_refine={4,8,16,32,64,128}` (matches `production_2026-08-23`
exactly, for direct comparability to its known l>=3 `kappa2_r` plateau),
`l_max=8`, and **96 shards/point** (960,000 configs/point -- above the
`n_blocks>=5*(2*l_max+1)=85` threshold needed for `chi2/dof` to actually
compute, 4x this campaign's usual shard count; short of the earlier "10M
configs" figure from before the storage-driven pivots, noted explicitly to
the user). New `cluster/sge/save_configs_shard_task.sh` /
`submit_save_configs_ladder.sh` (manifest-first, matches every other
production push's provenance pattern). **Job-array 7319121, 576 tasks,
submitted detached, completed cleanly**: `qacct` confirms
`exit_status=0` on all 576/576 tasks, zero error-indicating lines in any
log, exactly 96/96 shards present at every ladder point. Total storage
27 GB (matches the ~26 GB estimate), dominated by n_refine=128 (20 GB) and
n_refine=64 (5 GB). (Minor harmless bug noted and fixed in
`submit_save_configs_ladder.sh`: the pre-creation directory loop used
`seq -w` zero-padded shard indices piped through `printf '%03d'`, which
bash's builtin `printf` misparses as octal for "08"/"09" -- cosmetic only,
since each task independently creates its own directory with an
unambiguous plain integer; fixed for future runs via `$((10#$s))`.)

**Real bug found and fixed while running the M_l analysis against this
ladder**: `analyze_saved_configs.py`'s `compute_Ml_full` computed
`S = spins @ A` where `spins` is real float64 and `A` is complex128 --
numpy upcasts the entire real `spins` operand to complex128 for a
mixed-dtype matmul, which at n_refine=128 (`spins` is 13 GB real) means
allocating a ~26 GB complex *copy* of `spins` just to multiply by an
`(n_sites, n_lm)` matrix whose result is 1000x smaller than the copy.
Fixed by splitting into two real matmuls (`spins @ A.real`,
`spins @ A.imag`, recombined as `S = ... + 1j*...`), avoiding the giant
intermediate entirely. This affects both `analyze_saved_configs.py`'s own
CLI and `build_ylm_matrix_from_configs.py` (imports the same function).

**Separately discovered while timing this fix**: the interactive
session/login-node environment has `nproc=1` -- single-threaded numpy, so
the `M_l` matmul at n_refine=128 (`(10000,163842) @ (163842,81)`, ~2.7e11
flops) is compute-bound at an estimated ~2-3 minutes *per shard* even
after the memory fix, purely from lacking multi-core BLAS parallelism.
At 96 shards/point this is ~4-5 hours serially on this single core --
not something to run interactively in this session. **Next step**: submit
the `M_l`/`R_l`/`kappa2_r` and north-pole-two-point-function analyses as
SGE batch jobs (compute nodes, likely multi-core BLAS) rather than running
them on the login node, for the n_refine=64/128 points at least; smaller
ladder points (n_refine<=32) may be fast enough to run directly. Not yet
done as of this entry -- the production MC data (job 7319121) is complete
and validated, but neither analysis has been run to completion across the
full ladder yet.

## F_l normalization bug found: the l-space method DOES work, once a
spurious (2l+1) factor is removed (2026-08-26, same day)

**User request**: "can you also do a projection onto ylms and try the
delta = 2F0/F1 method" -- project the north-pole two-point function onto
Legendre space and compute `Delta_s(a) = 2*F_1/(F_1+F_0)`, the l-space
shortcut this campaign's earlier work (`cft_symmetry_test.py`,
`reference/Owen_section_D.tex`) already defines, but doing it from the
already-built north-pole data (single reference row) instead of the
expensive full `M_l` matmul pipeline.

Built `scripts/fl_delta_from_north_pole.py`: added an `n_targets` field to
`analyze_north_pole_configs.cc`'s output (count of target SITES per
equivalence class, distinct from the existing `orbit_size` which counts
PAIRS in the 120-element MC-averaging orbit) to serve as a per-class
quadrature weight, then computed `F_l = (2l+1)/(4*pi) * (1/n_sites) *
sum_classes n_targets * P_l(cos theta) * mean` -- carrying over the
`(2l+1)` prefactor from `cft_symmetry_test.py`'s `Trace(M_l)/(2*pi)`
convention, reasoning (at the time) that a single symmetric-exact
reference row should equal `Trace(M_l)/n_sites` up to an l-independent
constant. **First result: `Delta_s(a) = 0.32-0.34` at n_refine=4/8/16/32
-- flat, not converging, reproducing the exact plateau the campaign's
2026-08-25 "BREAKTHROUGH" finding attributed to UV contamination.**

**User pushed back twice, correctly**: "the F method should work for a
spherical CFT". This prompted actually checking the claim against data
rather than accepting the established narrative:
- Windowing test (excluding short-theta classes from the sum): `Delta_s`
  moved sharply (0.34 -> 0.29 -> 0.23 -> 0.18 -> 0.07 -> 0.01 as
  `theta_min` increased from 0 to 0.44) then blew up through a sign
  change past `theta_min~0.5` -- consistent with UV contamination *some*
  role, but the wild divergence (not a clean plateau) didn't fit a simple
  "just exclude the UV region" story.
- Direct point-by-point comparison of the measured correlator against the
  exact `A*(2-2cos theta)^{-1/8}` power law (amplitude fit from the bulk)
  showed **excellent agreement (within 1-3%) at every theta, including
  both poles** -- at n_refine=32, ratio measured/model stayed in
  [0.997, 1.030] across the ENTIRE range. This directly contradicted the
  "short-distance saturation" explanation: the correlator's SHAPE already
  matches Delta=1/8 almost everywhere, so a shape-driven l-space estimator
  should too.
- **Decisive test**: fed the theta grid + the EXACT synthetic `Delta=1/8`
  power law (not measured data) through `compute_Fl`. Result:
  `Delta_s(synthetic exact model) = 0.339` -- matching the plateau almost
  exactly, on data with NO lattice artifacts, NO statistics, NO UV
  contamination of any kind. Cross-checked against a true continuum
  integral (`scipy.integrate.quad` of `(2-2z)^{-1/8}*P_l(z)`, no lattice
  discretization at all): also `Delta_s = 0.33333` exactly. **This proved
  the discrepancy was in the formula, not the physics.**
- Read `reference/Owen_section_D.tex` directly (per `directive.md`'s
  "never hardcode a closed-form formula from memory, derive/verify
  first") rather than continuing to guess. Its `F_l^cont` is defined as
  `int_{-1}^1 dz (2-2z)^{-Delta} P_l(z)` -- **no `(2l+1)` factor anywhere**
  in that definition, unlike the separate `F_l(a) = (2l+1)/(2*pi*4*pi) *
  sum_ij sqrt(g_i g_j) P_l(...) <s_i s_j>` lattice formula (which needs
  the `(2l+1)` because it's built from `Trace(M_l)`, a sum over `2l+1`
  values of `m` reduced via the spherical-harmonic addition theorem).
  Numerically integrating the paper's own `F_l^cont` formula (no `2l+1`)
  gives `Delta_s = 0.125000` exactly, and `F_0` matches the paper's closed
  form `2^(1-2*Delta)/(1-Delta)` exactly.

**Root cause, confirmed**: a single reference row's Legendre projection of
the shape function `<s(theta)>` is the SAME kind of object as
`F_l^cont` (both are functions of one separation angle, integrated with
weight `P_l` and measure `dz`) -- it does NOT need the `(2l+1)` factor,
which only enters when reducing a sum over `m` (i.e. building the full,
bilinear `Trace(M_l)` from ALL site pairs) via the addition theorem. This
was a straightforward normalization bug in the new script, carried over
by incorrectly assuming the single-row and full-trace constructions share
identical prefactors.

**Fixed `compute_Fl`/`compute_Fl_jk` in `fl_delta_from_north_pole.py`
(removed the `(2l+1)` factor) and re-ran against real production data.
Clean result**: `Delta_s(a) = 0.1196 (n_refine=4), 0.1284 (8), 0.1289
(16), 0.1282 (32)` -- already within a few percent of 1/8 by n_refine=8,
essentially flat, no continuum extrapolation needed to get close (unlike
the real-space method, which needs the extrapolation to close the gap
from 0.139 at n_refine=4 down to ~0.125). `Delta_l_pair` at l=2..5 shows a
smaller but real residual trend (0.08-0.68 depending on l and resolution)
-- higher-l pairs are noisier/more resolution-sensitive, and l=6 remains
contaminated by the icosahedral invariant as in every earlier finding.
Continuum extrapolation of the (n_refine>=8) points gives
`Delta_inf~0.128` at both power=1 and power=2, ~2-3% high with poor
chi2/dof given how tight the individual-point errors now are (0.1284,
0.1289, 0.1282 aren't perfectly monotonic/smooth) -- a real, small,
not-yet-explained residual, not swept under the rug. Plot:
`campaign_runs/save_configs_ladder_2026-08-26/plots/
fl_vs_realspace_delta.png` (F_l method vs. real-space method, both vs.
1/n_refine).

**Open question this reopens, not yet resolved**: this specific bug was
in the NEW north-pole/single-row F_l script only. It says nothing
directly about whether the ORIGINAL 2026-08-25 finding (l-space
`Delta_l_pair`/two-operator fit built from `ising_s2_crit.cc`'s
IN-SIMULATOR, full bilinear `Trace(M_l)` -- a different code path,
correctly using the `(2l+1)` factor for that different construction) was
itself sound, or whether it has an analogous but distinct normalization
issue. The `M_l` SGE job (7322113, `build_ylm_matrix_from_configs.py`,
using `Trace(M_l)`-based `F_l`/`Delta_s`/`Delta_l_pair` code newly added
this session, reusing `cft_symmetry_test.py`'s formulas) is the natural
cross-check once it completes -- if IT also now gives `Delta_s~0.125`
rather than the historical `~0.33`, that would suggest the *original*
in-simulator conclusion (and the whole "l-space methods are fundamentally
unreliable due to UV contamination" pivot from 2026-08-25) may itself
need re-examination. If it still gives `~0.33`, that would mean the two
constructions have genuinely different failure modes and the original
finding stands on its own. **Do not assume either outcome -- check when
the job completes.**

**Partial cross-check done same session (n_refine=4, 8 shards, login
node -- ahead of the full SGE job 7322113, which was already too far
into its expensive n_refine=64/128 tasks by the time the F_l code landed
in `build_ylm_matrix_from_configs.py` to safely resubmit)**: ran the
UPDATED `build_ylm_matrix_from_configs.py` (true `Trace(M_l)`-based
`F_l`/`Delta_s`, correctly keeping the `(2l+1)` factor for that
construction) directly. **Result: `Delta_s(a) = 0.365 +/- 0.002` at
n_refine=4** -- still far from 1/8, NOT matching the fixed single-row
method's `0.1196` at the same resolution. **The two constructions
genuinely disagree on real (imperfect-symmetry) lattice data, despite
being mathematically identical in the exact-SO(3)-symmetric limit.** This
means: (1) the single-row bug fix above is real and specific to that
script, not a blanket fix for every l-space estimator; (2) the *original*
2026-08-25 in-simulator finding (full bilinear `Trace(M_l)`, `Delta_l_pair`
plateauing near ~0.33-0.6) is NOT explained by the `(2l+1)` bug and
appears to be a separate, still-standing anomaly; (3) full ladder data
for the `Trace(M_l)` method (n_refine=8..128, still finishing on job
7322113) is needed to see whether IT converges toward 1/8 with resolution
the way the fixed single-row method already does, or plateaus like the
original in-simulator method did. **Next step, not yet done**: get the
full-ladder `Trace(M_l)`-based `Delta_s(a)` trend once job 7322113
finishes, and investigate WHY the full bilinear sum (averaging over every
possible reference site, not just the north pole) disagrees with the
single well-chosen reference row -- candidate hypotheses not yet tested:
(a) some reference sites are much closer to mesh defects/non-uniformity
than the north pole (a maximally-symmetric vertex) and contaminate the
full trace; (b) a genuine remaining bug in how `compute_F_l`/
`diag_sums_per_l` reduce the block `sum` dict (worth an independent
hand-derivation check, not just trusting it reused correctly); (c) the
single-row method, by construction, only ever sees pairs related to one
reference point by the exact icosahedral orbit, which could bias it in a
way that happens to look good at small n_refine but doesn't reflect the
true bulk average -- not yet ruled out, though the clean, resolution-
independent flatness (0.1196 to 0.1289 across n_refine=4 to 64, tight
errors) argues against a large uncontrolled bias.

## RESOLVED: the north-pole-vs-Trace(M_l) disagreement above was a missing
(2l+1) division in Trace(M_l)'s F_l(a), not a genuine method disagreement
(2026-08-26, later same session)

**Context**: two side quests converged on this. First, per user direction
("do a small simulator run with these operators accumulating" / "expectation
values of (l+1)^2 operators plus powers of the magnetization to keep the
data volume lean"), built `src/lean_harmonic_stats.cc` + `bin/
analyze_lean_harmonic_stats.py` -- a new online-accumulating MC driver that
computes |S_lm|^2 for all (l_max+1)^2 modes and magnetization powers
directly during the sweep (no raw configs, no full M_l[m,m'] matrix ever
written), validated by reproducing `build_ylm_matrix_from_configs.py`'s
known Delta_s~0.365 baseline at n_refine=4 (got 0.3595+/-0.012, consistent).

Then, per user direction ("check the logic for the delta finder, there is a
factor of 3/8 or something messing up. the method should work exactly in
the continuum spherical cft, check the analytic derivation of this") --
re-derived the continuum limit of `cft_symmetry_test.py`'s `F_l(a) =
Trace(M_l)/(2*pi)` from scratch (three independent ways: direct P_l
double-integral, and a cleaner spherical-harmonic-orthonormality
derivation, both agreeing), rather than trusting the module's existing "(*)
verified 2026-08-25" derivation, which had only checked the ALGEBRAIC
identity `Trace(M_l) = (2l+1)/(4pi)*sum_ij w_iw_j P_l <s_is_j>` (an exact
lattice identity via the addition theorem) but never actually took that
lattice sum to its continuum limit and compared against Owen_section_D's
F_l^cont.

**Finding**: `F_l(a) := Trace(M_l)/(2*pi)` is a DOUBLE sum/integral over
site pairs (unlike the single-row north-pole method's F_l, a single
integral). Expanding the two-point function in the same Y_lm basis and
using spherical-harmonic orthonormality to collapse the double integral
shows `F_l(a) -> (2l+1) * F_l^cont` in the continuum limit, NOT `F_l^cont`
itself. `delta_s_formula`/`delta_l_pair_formula` were derived from (and
only correctly invert) the recursion `F_l^cont = [(l-1+Delta)/(l+1-Delta)]
* F_{l-1}^cont` -- feeding them `Trace(M_l)/(2*pi)` directly multiplies the
l=1/l=0 ratio by an unaccounted `(2*1+1)/(2*0+1) = 3`. At Delta=1/8: 3 x
0.125 = 0.375 -- matching the ~0.33-0.4 plateau this whole "l-space is
unreliable" pivot was built on, almost exactly. l=0's own factor is
`(2*0+1)=1` (no correction there), which is why `F_0` alone never looked
anomalous.

**This is a DIFFERENT bug from the single-row north-pole F_l's earlier
(2l+1) issue** (that script's error was an ERRONEOUS EXTRA (2l+1) that
needed REMOVING, since a single-reference-row projection is already
literally F_l^cont, no addition-theorem double-sum step involved). This
one is the OPPOSITE: the full bilinear Trace(M_l) construction is missing
a REQUIRED division by (2l+1). Both are real, and they don't cancel each
other or point to the same root cause -- confirmed by working through the
continuum limit for each construction separately, not by pattern-matching.

**Fix**: `cft_symmetry_test.py`'s `compute_F_l` now divides by `(2l+1)`
(full derivation in its updated docstring); `analyze_lean_harmonic_stats.py`
's local `F_l_from_abs2` was given the same fix (it doesn't import
`compute_F_l`, re-implements the same reduction locally).
`build_ylm_matrix_from_configs.py`/`build_ylm_matrix_from_corr.py` import
`compute_F_l` directly, so they're fixed automatically -- no separate edit
needed there, though their docstrings/print statements still say "EXPECTED
to be unreliable" (stale, not yet corrected).

**Confirmed numerically, two ways**:
- `lean_harmonic_stats` n_refine=4 test run (2500 configs): `Delta_s`
  0.3595+/-0.012 -> 0.1362+/-0.005 after the fix.
- `build_ylm_matrix_from_configs.py` rerun against the FULL
  `save_configs_ladder_2026-08-26` n_refine=4 data (96 shards, 960,000
  configs, high precision): `Delta_s(a) = 0.13955 +/- 0.00024`. This
  **matches the real-space method's independently-derived n_refine=4
  point (~0.139) to 3 decimal places** -- two structurally unrelated
  methods (real-space windowed fit vs. full bilinear Trace(M_l), now
  correctly normalized) now agree tightly. Right direction too (above 1/8,
  same as every other correctly-implemented estimator in this campaign).

**Reopens the 2026-08-25 "l-space is fundamentally unreliable" pivot**:
that conclusion was built on the buggy (missing (2l+1)) `Delta_l_pair`/
`Delta_s(a)` plateau. With the fix, l-space and real-space now agree at
n_refine=4 -- the original UV-contamination story needs re-examination
across the rest of the ladder (n_refine=8/16/32/64) before trusting either
framing. **Not yet done**: rerun the full ladder (n_refine=8/16/32/64, and
job 7322113's n_refine=64/128 once it finishes) through the fixed
`compute_F_l` and compare against the real-space continuum extrapolation;
update CLAUDE.md's "Real-space Delta_sigma extraction" pivot section,
which currently states the l-space method is unreliable as settled fact.
`OptimizeIntegrator`'s separate negative-weight crash bug (n_refine>=32,
see the north-pole entries above) is still unresolved and unrelated to
this fix.

## Ladder-wide confirmation of the (2l+1) fix using bin/lean_harmonic_stats
(2026-08-26, later same session)

Per user direction ("run it for 4,6,8,10,12,14,16,20,24,28,32 and make 1/L
plots for the various observables"), ran the new online-accumulating
`bin/lean_harmonic_stats` driver (q=5, l_max=8, `exact_sinh` coupling,
`UpdateWeights` quadrature, single seed=42 run per point, n_therm=2000,
n_traj=20000, n_skip=2 -> 10,000 measurements/point, jack_block_size=200
-> 50 blocks/point) across n_refine={4,6,8,10,12,14,16,20,24,28,32}
(162 to 10,242 sites). Total wall time ~5.5 minutes on the login node
(largest point, n_refine=32, took ~94s) -- confirms the lean driver is
fast enough for exploratory ladders without SGE.

**Result**: `Delta_s(a)` (post-(2l+1)-fix) decreases monotonically-ish
from `0.1418+/-0.0020` (n_refine=4) toward `0.125` as resolution
increases: 0.1342, 0.1365, 0.1322, 0.1292, 0.1320, 0.1313, 0.1264, 0.1235,
0.1249, 0.1276 at n_refine=6/8/10/12/14/16/20/24/28/32 respectively --
right sign (above 1/8, shrinking with resolution) and right rough
magnitude for a leading finite-lattice-spacing correction, consistent
with a 1/8 continuum limit within ~1-2 jackknife errors across the whole
ladder. This is a SECOND independent confirmation of the `compute_F_l`
(2l+1) fix above (the first being the single high-statistics n_refine=4
point, 0.13955+/-0.00024) -- now spanning the resolution range rather
than one point, with the same expected trend direction/shape as the
real-space method's ladder.

Per-l `Delta_l_pair(a)` (l=1..8) scatter around ~0.11-0.15 with no sign of
the old ~0.33-0.6 plateau at any l or resolution (see
`plot_lean_ladder.py`'s `delta_l_pair_vs_invL.png`) -- the anomaly this
whole 2026-08-25 pivot was built on does not reappear anywhere in this
fixed pipeline.

New script: `scripts/plot_lean_ladder.py` -- reuses
`analyze_lean_harmonic_stats.analyze()` (now also returns a results dict,
not just prints) per ladder point, produces `delta_s_vs_invL.png`,
`delta_l_pair_vs_invL.png`, `F_l_vs_invL.png` (small multiples, l=0..8),
`magnetization_vs_invL.png` (<M^2>, Binder U4). Output:
`lean_harmonic_stats/ladder_2026-08-26/plots/`.

**Not yet done**: no continuum-extrapolation fit (power-law in 1/L,
jackknife-propagated) run on this ladder yet -- `fit_real_space_ladder.py`
has the pattern to follow. This run also used `UpdateWeights`, not the
(separately still-broken) `OptimizeIntegrator` quadrature. Per user
direction, a higher-statistics rerun of this same ladder is queued next
(see below/subsequent entry).

## Overnight lean_harmonic_stats push: ladder extended to n_refine=128,
stats doubled (2026-08-26, later same session)

Per user direction ("submit an array job" then "schedule a much larger
stats ladder with larger lattices for an overnight run"): moved
`lean_harmonic_stats` output from top-level `lean_harmonic_stats/` to
`campaign_runs/` (matches every other production run's convention --
`ladder_2026-08-26`'s exploratory run and its plots moved there too, see
`CLAUDE.md`). Added `cluster/sge/lean_harmonic_stats_task.sh` (one array
task per `n_refine` point, no shards -- the driver's own jackknife blocks
already give a per-run error estimate) and two submit scripts following
this campaign's "copy, don't edit history" convention:

- `submit_lean_harmonic_stats.sh` -- ladder 4..32 (11 points),
  `n_traj=100000` (5x the exploratory login-node run). Submitted as job
  7328454, cancelled before running (still `qw`, cluster congested) once
  superseded by the overnight push below -- never produced output, safe
  to disregard its job ID.
- `submit_lean_harmonic_stats_overnight.sh` -- ladder extended to
  `4:6:8:10:12:14:16:20:24:28:32:40:48:56:64:80:96:112:128` (19 points,
  matching this campaign's standard `save_configs_ladder_2026-08-26`
  upper range), `n_traj=200000` (2x the previous push, 100,000
  measurements/point), `jack_block_size=1000` (100 native blocks),
  `h_rt=12:00:00`. Runtime estimate from login-node timing (n_refine=32,
  10000 meas, 94s -> 9.4ms/meas, cost scales ~linearly with n_sites):
  n_refine=128 (163,842 sites) is the long pole at ~4.2 hours; tasks run
  in parallel across compute nodes so total wall time is bounded by the
  single longest task. Submitted as **job 7328515** (19 tasks),
  `campaign_runs/lean_ladder_overnight_2026-08-26/` (manifest + per-task
  output), `campaign_runs/lean_ladder_overnight_2026-08-26_logs/`.

**Not yet done**: job was still `qw` (queued, cluster congested -- ~5,268
jobs system-wide at submission time) as of this entry, not yet started.
`OptimizeIntegrator`'s negative-weight bug (n_refine>=32) is irrelevant
here since this driver uses `UpdateWeights`, not `OptimizeIntegrator`.

**Superseded same session, per user direction "use more parallelism"**:
job 7328515 above was cancelled before it ever started (still `qw`) and
resubmitted with shard-level parallelism added. `cluster/sge/
lean_harmonic_stats_task.sh` now takes an `(n_refine, shard_index)` pair
per task (task-index layout: `shard_idx = global_idx % N_SHARDS`,
`l_idx = global_idx / N_SHARDS`, `seed = SEED_BASE + global_idx`),
mirroring `production_shard_task.sh`'s convention -- multiple
independent-seed shards per ladder point run as separate concurrent SGE
tasks (spread across more compute nodes) rather than one long serial task
per point. Output filenames are automatically collision-free across
shards (`lean_harmonic_stats.cc` names files by `<run_id>_lean_<seed_hex>.dat`,
and each shard gets a distinct seed), so no directory restructuring was
needed. `scripts/plot_lean_ladder.py` updated to group and pool all shard
files sharing an `n_refine` before calling `analyze()` (previously assumed
exactly one file per point) -- `analyze_lean_harmonic_stats.py`'s
`combine()` already supported multi-file pooling, no change needed there.

`submit_lean_harmonic_stats_overnight.sh` updated: `N_SHARDS=16` (matches
`submit_production.sh`'s convention), same ladder (19 points, 4..128) and
per-shard `n_traj=200000`, giving `19*16=304` total tasks -- 16x the
concurrent parallelism AND 16x the total statistics per point vs. the
un-sharded submission, at the same per-shard wall time (~4.2h estimated
for n_refine=128's shards, unchanged since N_SHARDS doesn't split
n_traj). Resubmitted as **job 7328520** (304 tasks),
`campaign_runs/lean_ladder_overnight_2026-08-26/` (same output directory,
manifest now also records `n_shards`/`seed_base` instead of a single
`seed`). Check `qstat -u misra` / `qacct -j 7328520` for status; once
complete, run `scripts/plot_lean_ladder.py
campaign_runs/lean_ladder_overnight_2026-08-26` for the 1/L plots.

**Note**: the prior session's `CLAUDE.md` accumulated substantial further
work (real off-diagonal `kappa2_r` added to `lean_harmonic_stats.cc`,
`--mesh_mode` ported into it, naive-vs-equal_area mesh-compare jobs, an
n_refine=512 ladder script, two SGE bugs found/fixed) that was never
back-filled into this journal — see `CLAUDE.md`'s "2026-08-27" entries for
that detail; not reproduced here to avoid duplication.

## 2026-08-27 (`/init` pass, session 6) — job status check + real kappa2_r reduction

Ran `/init`. Checked `qstat -u misra` and `qacct` against every job ID
`CLAUDE.md` said to check next session:

- **Job 7340503** (real off-diagonal `kappa2_r` at n_refine=64,128,
  naive mesh) — confirmed finished, `exit_status=0` on all 192/192 tasks.
  Ran `scripts/reduce_ml_blocks.py` against both resolutions' 96-shard
  pickle sets (`campaign_runs/save_configs_ladder_2026-08-26/ml_blocks/
  q5k{64,128}/`). **Result: `l=1,2` are noise-consistent
  (`kappa2_r<~1e-5`), `l>=3` plateaus at `kappa2_r~3e-5` to `7e-5` and
  does NOT shrink from n_refine=64 to 128** (l=3: `5.23(78)e-5` at 64 vs
  `6.82(103)e-5` at 128; chi2/dof 2.7-9.0 across l=3..8 at both
  resolutions) — the first confirmation of the l>=3-4 symmetry-breaking
  plateau using the campaign's actual designated gate (full off-diagonal
  `M_l`, not `symmetry_check_lean.py`'s diagonal-only proxy), at the
  largest `n_refine` ever reached with a full `M_l` matrix. `Delta_s(a)`
  from the same reduction: `0.12564(24)` (n_refine=64), `0.12573(26)`
  (n_refine=128) — both ~2.4 sigma of 1/8, consistent with objective (2)'s
  separate confirmation. This run is naive-mesh-only (job 7340503 predates
  `lean_harmonic_stats.cc` getting `--mesh_mode`) — does not yet test
  whether `equal_area` changes the plateau.
- **Jobs 7341057/7341058** (16-shard naive/eqarea mesh-compare,
  `campaign_runs/lean_mesh_compare_2026-08-27/`) — running, ~63/304 and
  ~102/304 tasks currently `r` respectively, rest presumably `qw`/pending
  — not complete, not yet analyzable.
- **Jobs 7341059/7341060** (100-shard naive/eqarea mesh-compare,
  `campaign_runs/lean_mesh_compare_100shards_2026-08-27/`) — 7341059 just
  starting (5/1900 tasks `r`), 7341060 still fully `qw` (0/1900 started)
  — far from complete.
- **`cluster/sge/submit_lean_ladder_512.sh`** (n_refine up to 512 ladder)
  — **discovered NOT actually submitted**, despite `CLAUDE.md` describing
  it as "submitted": no `s2prec_lean_ladder_512*` job in `qstat`, and
  `campaign_runs/lean_ladder_512_2026-08-27/` contains only the
  `manifest.json` the script writes before its first `qsub` call — no
  mesh cache, logs, or `.dat` files. The script itself looks complete and
  correct; the invocation that runs it appears to have been dropped at
  the end of the prior session. Corrected `CLAUDE.md` to say so. **Not
  resubmitted this session** (large job, ~30h wall time per shard,
  960 total tasks — deferred to explicit user direction rather than
  launched unilaterally).

Did not touch git — repo still has the same pre-existing uncommitted
working-tree state (`CLAUDE.md`/`PLAN.md`/`directive.md`/`journal.md`
modified, several `campaign_runs/` dirs untracked).

## 2026-08-28 -- naive-vs-equal_area real kappa2_r comparison, n_refine<=48 (job 7341057/7341058, partial)

Jobs 7341057 (naive, 16 shards) and 7341058 (equal_area, 16 shards) --
the mesh-compare production arrays matching `lean_ladder_overnight_
2026-08-26`'s stats -- are still running (task IDs are processed in
ladder-point order, one point's full 16 shards at a time), but every
ladder point up to and including n_refine=48 had already finished
(exit_status=0, all 16/16 shards) on both mesh modes as of this check;
n_refine=64 and above are still in flight. Built filtered symlink dirs
(`campaign_runs/lean_mesh_compare_2026-08-27/analysis_upto48/{naive,
eqarea}/`) containing only the confirmed-complete points and ran
`scripts/symmetry_check_lean_kappa2.py` (the real off-diagonal kappa2_r,
not the diagonal-only proxy) against each, plus a per-(n,l) ratio/sigma
table.

**Result: equal_area's l>=3 kappa2_r plateau is real and measurably lower
than naive's, but it is still a plateau, not a shrinking trend.**

- At coarse resolution (n_refine=4-10) equal_area crushes naive by
  2x-64x (tens of sigma) -- reconfirms push7's original finding at much
  higher statistics.
- From n_refine=12 out to 48 (the newly available range), naive/eqarea
  ratio settles to roughly 2.5-4x for l=3-5 and l=7-8 (4-7 sigma,
  consistently significant), l=6 a bit smaller (~2-3x, 3-5 sigma). This
  does NOT shrink toward 1 as n_refine increases across this range --
  equal_area stays systematically better, not converging with naive.
- **Neither mesh mode's own l>=3 plateau is shrinking with resolution
  in this range.** naive l=6: 2.90e-5 (n=12) -> 2.39e-5 (n=48), roughly
  flat. eqarea l=6: 1.09e-5 (n=12) -> 0.94e-5 (n=48), also flat. Same
  pattern at l=3,4,5,7,8 for both modes -- values wobble within a factor
  of ~2 but show no continuum-limit-consistent downward trend out to the
  largest common resolution checked so far.

**This contradicts the earlier (lower-stats) push7-session conclusion
that "both mesh modes converge to the same plateau by n_refine>=12-16"**
(2026-08-25 entry above) -- with real production statistics (16 shards x
200k traj, not push7's spot-check runs), equal_area's plateau is clearly
and persistently lower than naive's over this whole range, not identical.
Revises objective (1)'s open question: mesh non-uniformity (area at
least) is not a red herring after all -- it measurably helps -- but it
does not fully explain the l>=3 symmetry breaking, since even the
better-behaved equal_area mesh still has a nonzero, non-shrinking
plateau. The circumradius+perimeter (`equal_rp`) mesh, or a genuine
residual-lattice-artifact explanation, remain live candidates.

**Next step, not yet done**: once jobs 7341057/58 finish (n_refine up to
128) and the 100-shard companions (7341059/60) and the n_refine-up-to-512
ladder (7342648 chain) complete, rerun this same comparison at higher
n_refine to see whether either mesh mode's plateau eventually shrinks, or
whether both saturate at a genuine nonzero floor -- the current n<=48
window is not yet resolving enough of the ladder to distinguish "shrinks
very slowly" from "saturates."

## 2026-08-28 (cont) -- full n_refine=4..128 naive-vs-equal_area comparison
(jobs 7341057/7341058 finished), plus new delta_l_pair/delta_l plots

Jobs 7341057 (naive, 16 shards) and 7341058 (equal_area, 16 shards) finished
completely between the last check and this one (`qacct`: 304/304 tasks each,
`exit_status=0`) -- all 19 ladder points, n_refine=4..128, both mesh modes,
16/16 shards. Reran `scripts/symmetry_check_lean_kappa2.py` on the full
range (previously only checked up to n_refine=48) and ran the existing
`scripts/plot_mesh_compare_full.py` for the first time against the complete
dataset: `campaign_runs/lean_mesh_compare_2026-08-27/analysis_full/plots/`
(`mesh_compare_{R_l,kappa2_r,kappa3_r,kappa4_r}_vs_invL.png`,
`mesh_compare_scalar_vs_invL.png`, `mesh_compare_ratio_all_metrics_vs_invL.png`).

**kappa2_r plateau confirmed out to n_refine=128 for both mesh modes** --
extends the n<=48 finding from earlier today, no new behavior at higher
resolution: naive l=3..8 kappa2_r sits in the `2e-5`-`6e-5` range from
n_refine=64 to 128 (e.g. l=6: 2.42e-5 at 64, 2.38e-5 at 128 -- flat);
equal_area stays ~3-5x lower (l=6: 7.4e-6 at 64, 1.05e-5 at 128) but
likewise flat, not shrinking. Neither mode's plateau moves between
n_refine=64 and 128 -- the "does this eventually shrink or saturate"
question from the earlier entry is still open in the sense that 128 isn't
qualitatively different from 48, but it does rule out "n_refine=48 was
just not far enough yet" as an easy resolution -- 128 shows the same flat
floor. `Delta_s(a)`/Binder `U4` (the scalar CFT-data/bulk metrics, from
`plot_mesh_compare_full.py`'s second figure) are naive/eqarea-indistinguishable
at every resolution, unlike the symmetry-breaking metrics -- mesh mode
only matters for the l>=3 kappa_r family, not for the scalar observables.

**New script, per user request ("define what R_l is, kappa2 is the
important thing" then "make plots for all of the other observables ...
the delta from all pairs of harmonics from the recursion relation")**:
`scripts/plot_delta_pairs_mesh_compare.py`. Reuses
`cft_symmetry_test.py`'s `compute_F_l`/`delta_l_formula`/
`delta_l_pair_formula`/`jackknife_err` verbatim (the already-fixed,
(2l+1)-corrected versions -- see the 2026-08-26 "RESOLVED" entries above),
swapping only the block-loading step for the lean-data equivalent
(`analyze_lean_harmonic_stats.combine` + `to_analyze_blocks_format`,
the same swap `symmetry_check_lean_kappa2.py` already does for
R_l/kappa2_r) since `cft_symmetry_test.analyze()` itself expects raw
`save_configs`-style jackblocks files, not this lean-driver format.
Produces two new plots:
`mesh_compare_delta_l_pair_vs_invL.png` (Delta_{l,pair}(a), the
single-recursion-step scaling-dimension estimate, one panel per l=1..8
pair, naive vs eqarea vs 1/L) and `mesh_compare_delta_l_vs_invL.png`
(delta_l(a), the dimension-independent conformal-symmetry test, l=2..8).

**Result: every l pair's Delta_{l,pair}(a) converges cleanly to the exact
1/8 by n_refine~64-128, both mesh modes, l=1 through 8** -- e.g. l-1=0->l=1:
naive 0.1252 at n_refine=128 (vs 0.1398 at n_refine=4); l-1=3->l=4: naive
0.1229 at 128 (vs 0.1849 at 4); l-1=7->l=8: naive 0.1265 at 128 (vs 0.4345
at 4). All panels show a monotonic-ish approach to 0.125 with decreasing
1/L, naive and eqarea overlapping almost everywhere (mesh mode doesn't
matter for this diagnostic either, consistent with the Delta_s(a)/U4
finding above). This is a clean, independent confirmation of objective (2)
(Delta_sigma=1/8) using the recursion relation at every harmonic level, not
just l=1 (which is what Delta_s(a) alone checks) -- extends the existing
closure of objective (2) rather than reopening it. Genuinely striking
contrast with the kappa2_r result above: the SAME (l,l-1) harmonic pairs
that show a persistent, non-shrinking spherical-symmetry-breaking plateau
in kappa2_r show a cleanly-converging, textbook-looking CFT scaling
dimension in Delta_l_pair -- the l>=3 anomaly is specific to the
off-diagonal/eigenvalue-spread structure of M_l, not a general breakdown
of the harmonic-decomposition method at those l values.

Cluster context checked while this ran: ~9,100 jobs system-wide
(5,219 qw, 3,691 r) -- explains why 7341059/7341060 (the 100-shard
companions) and 7342648 (n_refine-up-to-512 chain) are moving slowly/not
started; not a stuck-job problem, genuine queue congestion plus the
inherent n_refine=128 cost (163,842 sites, ~9.4ms/measurement scaling
~linearly with n_sites).

**Not yet done**: delta_l(a) itself (the dimension-independent test, as
opposed to Delta_l_pair) has not been read/interpreted here, only plotted
-- worth checking whether it also shows clean convergence to 0 or has its
own l-dependent structure. The 100-shard/n_refine-512 jobs are still not
analyzable. No new production run was needed for any of this session's
work -- all reuses existing completed job output.

## 2026-08-28 (cont 2) -- kappa2 linear plot restricted to n_refine>32, new
Delta_{l,pair}(a) recursion-extrapolation plot (all l), delta_l(a) interpreted

**Per user request ("only plotted for L greater than 32")**:
`scripts/plot_kappa2_linear.py` gained `--min_n_refine` (default 32,
strict `>`) and now drops the coarse ladder points from
`mesh_compare_kappa2_r_vs_invL_linear.png` -- the plateau (naive
~2-6e-5, eqarea ~3-5x lower, both flat) is easier to read without the
large-1/L points compressing the y-axis.

**New `scripts/plot_delta_s_recursion_fit.py`**, per user request ("make
a delta_s(a) plot that shows the extrapolated predictions from the
recursion relation" then "for all harmonics"): one panel per l=1..8
pair, each showing naive/eqarea `Delta_{l,pair}(a)` data plus its own
`Delta_inf + c*n_refine^-p` fit (best of p=1/p=2 by chi2/dof, reusing
`fit_lean_ladder.fit_power`) extrapolated to `1/L=0`, marked with an
open-diamond point, against the exact-1/8 reference line. Output:
`campaign_runs/lean_mesh_compare_2026-08-27/analysis_full/plots/
delta_l_pair_recursion_extrapolation.png`. Confirms the existing
per-l-family-fit finding (2026-08-26) on the current full-statistics
naive/eqarea ladder: l=1 alone prefers p=1 (naive/eqarea `Delta_inf`
0.12487(7)/0.12492(7), chi2/dof 1.6); every l=2..8 prefers p=2, with
`Delta_inf` scattered 0.121-0.124 (mostly 1-3 sigma from 1/8, l=6 naive
an outlier at 27.6 chi2/dof/2.4 sigma low, likely a bad p=2 fit for that
level specifically rather than new physics -- consistent with the
2026-08-26 finding that no single power fits every l cleanly). Mesh
mode barely matters here, matching every other recursion-relation
diagnostic in this campaign.

**delta_l(a) interpreted for the first time** (previously only plotted,
per the prior entry's "not yet done" item): `mesh_compare_delta_l_vs_invL.png`
shows, for every l=2..8, `delta_l(a)` (the dimension-independent
consistency check that the F_{l-2},F_{l-1},F_l recursion actually closes,
not a Delta estimate) is statistically consistent with 0 across most of
the ladder (naive/eqarea overlapping, error bars mostly covering 0) and
only peels away sharply -- with an amplitude that grows steeply with
l (l=2 : ~0.01 max; l=8: ~-0.5) -- at the largest 1/L (coarsest,
smallest n_refine) values. This is the signature of a short-distance
lattice-discretization artifact that heals as the continuum limit is
approached, not a resolution-independent symmetry-breaking floor: unlike
`kappa2_r`'s l>=3 plateau (flat all the way to n_refine=128, never
shrinking), `delta_l(a)`'s deviation is concentrated at small n_refine
and vanishes well before n_refine=128. Higher l levels need finer
resolution before `delta_l` settles near 0 (their curves stay away from
0 out to larger n_refine than l=2/3's do) -- consistent with a UV/
short-distance-dominated correction that gets proportionally worse for
higher angular resolution on a fixed mesh, the same qualitative pattern
as the l>=4 continuum-extrapolation power (`p~2-2.7`) being steeper than
l=1's. **Net read: the recursion relation itself is well satisfied in
the continuum limit at every l** -- `delta_l(a)` closing to 0 both
confirms internal consistency of the Owen-section-D recursion (not just
that individual `Delta_l_pair(a)` extrapolate near 1/8) and reinforces
that the still-open `kappa2_r` l>=3 anomaly is a distinct phenomenon
(basis-dependent eigenvalue-spread of `M_l`, not a general recursion or
CFT-data-extraction breakdown -- same conclusion as the 2026-08-28
entry above, now independently supported by `delta_l` rather than just
`Delta_l_pair`).

## 2026-08-28 (cont 3) -- delta_l_pair recursion-extrapolation plot updated
to show p=2 and floated-power fits, per user request

Per user follow-up ("try using quadratic fits and also floating power
fits"): `scripts/plot_delta_s_recursion_fit.py` reworked to drop the
p=1/p=2-by-chi2 selection and instead show BOTH a fixed p=2 fit and a
floated-power fit (reusing `fit_lean_ladder_family.fit_single_power`/
`fit_floated_power` directly rather than reimplementing) on every panel,
naive and eqarea each getting their own two curves + extrapolated
`1/L=0` diamond. l=1's floated power still comes out near 1
(naive/eqarea p=1.10/1.07, matching the established l=1-prefers-linear
finding even though p=2 is also shown for comparison); l>=3 floated
powers cluster 2.2-3.0, consistent with the existing family-fit result.
**l=2 is degenerate**: its floated-power fit pins at the upper bound
`p=5.00` for both mesh modes (`fit_floated_power`'s bounds are
`[0.1, 5.0]`) -- l=2's `Delta_l_pair(a)` is nearly flat/non-monotonic
over most of the ladder (see its panel: values wander within ~0.002 of
1/8 with no clear power-law curvature until the very largest 1/L), so
the floated fit has no real power-law signal to lock onto and just
saturates the bound rather than converging -- not a real p~5 discretization
effect, flag as a fit artifact specific to this level, not a physical
finding.

## 2026-08-28 (cont 3) -- equal_rp added to 100-shard lean_harmonic_stats mesh-compare, per user direction

User asked to schedule a large stats run using "equalized area, circumcenter
and perimeter" -- clarified this maps to the existing `equal_rp` mesh mode
(joint circumradius+perimeter+area equalization via
`QfeLatticeS2::EqualizeCircumPerim`, arXiv:2407.00459 Appendix B.1/B.2),
which had only ever been run through `ising_s2_crit.cc` at low
resolution/stats (push8, n_refine<=32, inconclusive) -- never through the
high-stats `lean_harmonic_stats` driver used for the current naive/eqarea
100-shard and 512-ladder production. User confirmed scale via AskUserQuestion:
100-shard ladder to n_refine=128 (matching the existing naive/eqarea
100-shard run), not the 512-ladder.

Wrote `cluster/sge/submit_lean_mesh_compare_100shards_eqrp.sh` (clone of
`submit_lean_mesh_compare_100shards.sh`, third mesh mode added alongside
existing naive/eqarea in the same `lean_mesh_compare_100shards_2026-08-27/`
stage dir rather than a new stage -- own `eqrp/` output subdir, own
`_eqrp_logs/` dir, own `manifest_eqrp.json`, `SEED_BASE=820000` to stay
independent of naive(800000)/eqarea(810000)). Prewarm task script
(`prewarm_mesh_cache_task.sh`) and lean production task script
(`lean_harmonic_stats_task.sh`) already supported `MESH_MODE=equal_rp`
unchanged (added 2026-08-25/27 respectively) -- no code changes needed,
submit-script-only.

Submitted (background, prewarm blocks on `-sync y` up to 2h before
autolaunching the 1900-task production array
`s2prec_lean_mesh_compare_100sh_eqrp`, `h_rt=12:00:00`/task). Targets the
still-open Objective (1) spherical-symmetry gate (`kappa2_r` plateau flat
and non-shrinking for both naive and equal_area out to n_refine=128) --
equal_rp is the untested-at-scale hypothesis for whether joint
circumradius/perimeter/area equalization (vs area-only) closes that gap.

**Correction, 2026-08-29**: despite the above, this job **never actually
reached SGE** -- `qacct -j s2prec_lean_mesh_compare_100sh_prewarm_eqrp`/
`s2prec_lean_mesh_compare_100sh_eqrp` both come back "job name not
found", `qstat` showed nothing queued/running under either name, and
`campaign_runs/lean_mesh_compare_100shards_2026-08-27/eqrp/` plus its
`mesh_cache/*eqrp*` entries were completely empty. Only `manifest_eqrp.json`
and an empty `_eqrp_logs/` dir existed (both dated 2026-08-28 16:13) --
the submit script ran far enough to write its manifest and `mkdir -p`
its output dirs, then apparently died before `qsub` ever successfully
contacted the scheduler (no error trace survived to explain why). Same
"described as submitted but wasn't" failure mode as the L=512 job above.
**Resubmitted 2026-08-29**: re-ran
`cluster/sge/submit_lean_mesh_compare_100shards_eqrp.sh` directly
(nohup'd, detached). Prewarm array landed as **job 7357536** (19 tasks) --
confirmed actually running this time. Once its `-sync y` wait completes
and all 19 cache files verify, the same script auto-launches the
1900-task `s2prec_lean_mesh_compare_100sh_eqrp` production array with no
further action needed. Check `qstat -u misra` / `qacct -j 7357536` before
assuming either step is done.

## 2026-08-29 -- naive-vs-equal_area 100-shard kappa2_r analysis: crash root cause found, higher-stats plateau confirmed through n_refine=112

The naive/eqarea 100-shard production (job 7341060 + naive counterpart)
finished cleanly days ago (1900/1900 tasks each, verified complete and
uniform on disk at every ladder point including n_refine=128), but the
one analysis job run against it (`s2prec_plot_mesh_compare_backup`,
7353395, via `scripts/plot_mesh_compare_full.py`) crashed with a
`KeyError: (5, -1, 1)` while processing n_refine=128/eqarea, after
successfully printing real `kappa2_r`/`R_l` for every point from
n_refine=4 up through 112. **Root cause**: job 7353395 happened to run
its analysis while the eqarea production job's slowest shards (at
n_refine=128) were still being written -- a known torn-write hazard
(`analyze_lean_harmonic_stats.load()`'s docstring already documents a
2026-08-28 instance of this against job 7341060 specifically), but this
particular manifestation wasn't the malformed-float case that function's
`except` clause catches: a line flushed mid-write can be cut exactly at
a field boundary, giving fewer than the expected 975 whitespace-
separated columns while every individual token still parses as a valid
float -- no exception is raised, `to_analyze_blocks_format`'s
`zip(offdiag_list, blk["offdiag_re"], blk["offdiag_im"])` just silently
stops early, and the dict is missing whichever `(l,m,mp)` keys came
after the truncation point. Confirmed the data itself is NOT corrupted:
independently re-checked all 100 current `q5k128_eqarea_lean_*.dat`
files, every one has exactly 100 blocks x 975 fields, uniform, complete
-- the race was real but transient, and production has been done for
days now.

**Fix applied**: just reran the analysis (no code change yet) now that
production can't be mid-write. Two versions launched, both directly on
the login node (not via SGE, so no repeat of the original race):
`scripts/plot_mesh_compare_partial.py` (new script, same plotting code
as `plot_mesh_compare_full.py` imported verbatim, just adds a `--max_n`
cutoff so a full quartet of plots -- R_l/kappa2_r/kappa3_r/kappa4_r +
Delta_s/U4 scalar + naive/eqarea ratio-summary, all with real jackknife
error bars -- can be produced from n_refine<=112 without waiting on the
expensive n_refine=128 point) and a plain rerun of
`plot_mesh_compare_full.py` (includes n_refine=128, much slower --
~10,000 blocks at that point vs ~100-1000 at the smaller ones). Both
were still running as of this entry; check
`campaign_runs/lean_mesh_compare_100shards_2026-08-27/plots_partial/`
and `.../analysis_full/` for output.
**Not yet done**: harden `load()` to also drop a line with fewer than
the expected column count (not just a literal float-parse failure) --
the exact edge case that caused this crash is still unfixed in the
shared library code, just avoided this time by not reading during a
write. Worth doing before the equal_rp 100-shard job (job 7357536,
above) lands, since its analysis will hit the same shared code path.

**Result (n_refine=4..112, confirmed clean, no jackknife-covariance
caveats beyond what analyze_blocks already provides)**: extends the
"obj 1" finding to much higher resolution and statistics than before.
naive `kappa2_r` at l=3 does not shrink with n_refine -- it visibly
*rises*: ~2.3e-5 (n=32) -> ~3.7e-5 (n=48) -> ~4.7e-5 (n=96) -> ~4.5e-5
(n=112), a genuine plateau/mild growth, not a coarse-mesh artifact.
equal_area's l=3 stays essentially flat around ~1e-6 across the same
range -- roughly **30-40x lower than naive** (higher than the earlier
"2.5-4x" figure quoted at n_refine<=48; the gap between the two mesh
modes widens at higher resolution, it doesn't narrow). Every l=3..8
level shows the same qualitative pattern for both mesh modes. Neither
mesh's plateau is shrinking toward zero at these resolutions --
equal_area suppresses the symmetry-breaking amplitude by over an order
of magnitude but does not eliminate it, matching (and now
higher-statistics-confirming) the standing conclusion that mesh choice
alone doesn't close the spherical-symmetry gate. The n_refine=128 point
and the equal_rp comparison (job 7357536/its still-pending production
array) are the two open pieces still needed for the complete picture.

## Extended mesh-compare plots to all magnetization moments + refreshed per-l Delta plots (2026-08-29)

Per user direction ("there are more observables to plot: deltas(a) from
all harmonics delta_l as well as magnetization and powers of the
magnetization and all"): `analyze_lean_harmonic_stats.load()` was already
parsing `M_sum`/`M3_sum` into `blk["M"]`/`blk["M3"]` per block, but
`analyze()` never used them -- only `<M^2>`, `<M^4>`, and Binder `U4` were
ever exposed, despite the first and third moments and susceptibility being
free (already-accumulated) data. `analyze()` now also computes `<M>`,
`<M^3>`, and susceptibility per site `chi = (<M^2>-<M>^2)/n_sites`, each
with the same delete-one-block jackknife error as the existing moments.

New **`scripts/plot_mesh_compare_magnetization.py`** (mirrors
`plot_mesh_compare_full.py`'s structure/style): one small-multiples figure
with all six magnetization observables (`<M>`, `<M^2>`, `<M^3>`, `<M^4>`,
`chi`, `U4`), naive vs equal_area vs `1/L`, plus a naive/eqarea ratio
figure for the three strictly-positive moments (`<M^2>`, `<M^4>`, `chi`).
Also reran the pre-existing `scripts/plot_delta_pairs_mesh_compare.py`
(`Delta_l_pair(a)`/`delta_l(a)`, all l=1..8) against the same dataset --
it existed but had not yet been regenerated against the 100-shard run.

Both run against `campaign_runs/lean_mesh_compare_100shards_2026-08-27/
{naive,eqarea}` (confirmed no production job was writing to that directory
at the time, avoiding the torn-write race documented in the entry just
above this one) -> output under `.../analysis_full/`:
`mesh_compare_magnetization_vs_invL.png`,
`mesh_compare_magnetization_ratio_vs_invL.png`,
`mesh_compare_delta_l_pair_vs_invL.png`, `mesh_compare_delta_l_vs_invL.png`.

**Result (n_refine=128, 100 shards x 100,000 meas = 10M measurements/mesh,
the highest-stats point in this dataset)**: `<M> = -1.44e-04 +/- 1.41e-04`
and `<M^3> = -3.42e-05 +/- 3.83e-05` on eqarea (naive: `<M> = -2.01e-04`) --
both consistent with zero on both meshes, as expected from the model's
exact Z2 symmetry (a real, not just assumed, null check). `<M^2>~0.200`,
`<M^4>~0.0521`, `chi = 1.222e-06 +/- 2.8e-10`, `U4 = 0.5671 +/- 0.0001` --
naive and eqarea agree to ~3-4 significant figures on every one of these,
no resolvable mesh-mode discrepancy anywhere in the ladder. `Delta_l_pair`
at l=1..8 clusters at 0.124-0.126 for both meshes except **l=7, which
reads low on both** (naive 0.1224, eqarea 0.1242) -- consistent with the
existing `fit_lean_ladder_family.py` finding (2026-08-26 entry, "l=7 is a
genuine anomaly") that l=7 is an outlier across multiple independent
diagnostics, not new here.

**Conclusion**: unlike `kappa2_r`/`R_l` (real, large, resolution-
independent naive-vs-eqarea gap at l>=3, per the entry above), the bulk
magnetization moments and the diagonal-only per-l CFT-scaling-dimension
estimators (`Delta_l_pair`, `delta_l`) show no resolvable naive-vs-
equal_area difference at this statistics level. The mesh-mode-dependent
symmetry-breaking signal this campaign has been chasing lives specifically
in the off-diagonal `M_l[m,m']` channel (`kappa2_r`/`R_l`/`kappa3_r`/
`kappa4_r`), not in these diagonal/bulk observables -- a useful negative
control, not just an incidental extra plot.

## 2026-08-29: eqrp 100-shard production landed; found and recovered a
## real stuck-job bug in the ladder_512 resubmission; kappa2_r 3-way
## comparison (naive/eqarea/eqrp) submitted, not yet analyzed

**eqrp production (job 7358337, 1900 tasks) finished clean at 11:08am** --
the third mesh mode's data for the `lean_mesh_compare_100shards_2026-08-27`
ladder is now fully on disk (`.../eqrp/`, 1900 `.dat` files), alongside the
already-analyzed naive/eqarea pair. Not analyzed same-session as landing;
see below.

**Real bug found: `resubmit_ladder512_n512_prewarm_and_production.sh`'s
`qsub -sync y` for the missing n_refine=512 prewarm point silently hung
for 10h41m after its job actually finished.** The prewarm job itself
(7357578) completed cleanly at 06:02am (`qacct: exit_status=0`, cache file
`q5k512_eqarea_step0.300.dat` present and correct) -- but the `-sync y`
client process on the login node never returned (`qstat -j 7357578` said
"job does not exist" while the client was still blocking), so **step 2 of
that script -- submitting the actual ladder_512 naive/eqarea production
arrays -- never ran**, despite the script's own log implying it might have
(only the header line of `_resubmit_2026-08-29.log` had been written,
consistent with a hang immediately after step 1's `-sync y` call, not a
crash). Root cause not fully diagnosed (a repeat of the same qmaster
instability documented in `submit.log` from the original 2026-08-28
attempt is the leading suspect, not confirmed). **Recovery**: verified all
15 mesh-cache files were genuinely present (`ls`, not just trusting the
script), killed the two hung processes (`bash resubmit_...sh` and its
child `qsub -sync y`, PIDs 2338246/2338258), then ran step 2's two `qsub`
calls directly by hand (same flags the script would have used) --
**submitted jobs 7363384 (naive) and 7363385 (eqarea)**, 480 tasks each
(15 ladder points x 32 shards), `h_rt=30:00:00`. Both `qw` (queued,
cluster congestion) as of this entry, not yet started.

**Lesson for next session, extending the existing 2026-08-27 "verify
independently after -sync y" lesson**: that lesson covered the case where
`-sync y` returns but downstream logic is silently wrong; this is the
complementary failure mode -- `-sync y` can also simply never return even
though the job it's waiting on has already finished server-side. Any
`-sync y`-gated recovery script should be checked with `ps`/`qstat -j
<id>` if it's been running suspiciously long, not assumed to still be
faithfully blocking.

**New `scripts/compare_three_mesh_kappa2r.py`** (thin wrapper around
`symmetry_check_lean_kappa2.point_kappa2`, no reimplementation) extends
the existing naive-vs-eqarea `kappa2_r` comparison to all three mesh
modes now that eqrp's data exists. First attempt run interactively on the
login node -- resource-killed (SIGKILL, exit 137) after ~11 minutes on the
n_refine=128 point (~10,000 blocks/mesh, the same cost that made the
naive/eqarea pair slow) -- moved to a proper detached SGE job per the
workspace-wide "prefer qsub over long interactive runs" convention. New
**`cluster/sge/compare_three_mesh_task.sh`**, submitted as **job 7363395**
(single task, `h_rt=4:00:00`), `qw` as of this entry. **Not yet run/no
numbers yet** -- check `qacct -j 7363395` and the job's log under
`campaign_runs/lean_mesh_compare_100shards_2026-08-27/three_mesh_compare_logs/`
next session; this is the concrete next step to actually resolve whether
`equal_rp` shrinks the l>=3 `kappa2_r` plateau further than `equal_area`
already does (open since the 2026-08-25 `equal_rp` mesh mode was built
specifically to test this).

## 2026-09-02: lean_ladder_512 (n_refine up to 512) analyzed for the first
## time -- gate still FAILS at l>=3, plateau confirmed out to the largest
## mesh this campaign has ever reached with the full off-diagonal M_l

Jobs 7363384 (naive, 480 tasks) and 7363385 (eqarea, 480 tasks) had
finished (per the 2026-08-29 entry) but were never actually analyzed --
`campaign_runs/lean_ladder_512_2026-08-27/` had only raw shard files, no
plots. Ran `scripts/symmetry_check_lean_kappa2.py` and
`scripts/plot_mesh_compare_full.py` against both directories for the
first time this session.

**Job completion, checked via `qacct`**: naive 479/480 tasks exit_status=0
(1 task exit 137, walltime kill) but **all 32 n_refine=512 shard files are
present and complete** (39 lines each = 4 header + 35 full jackknife
blocks) -- the walltime-killed task apparently finished writing before the
kill landed. eqarea: 461/480 exit_status=0, 19 tasks exit 137, all
concentrated in the n_refine=512 point (tasks 460-480) -- **3 of the 32
n_refine=512 shard files had a truncated final block** (mid-`fprintf` cut
by the walltime kill: 859/851/763 fields instead of 975 on the last line,
every other line and every other shard's line count correct). Fixed by
`sed -i '$d'` on just those 3 files (dropped the one incomplete trailing
block each, keeping 33/33/27 complete blocks respectively instead of the
full 35) -- `scripts/symmetry_check_lean_kappa2.py`/`analyze_lean_harmonic_stats.load()`
still don't harden against this (the 2026-08-27 "not yet done" TODO), so
without this fix the eqarea n_refine=512 point crashes with `KeyError` on
whichever `(l,m,mp)` pair happens to fall past the truncation point.
**None of this touches the campaign's earlier n_refine<=384 results** --
those files were already complete.

**Result**: `Delta_s(a)` at n_refine=512 (largest point, ~1.1M
measurements/mesh): naive `0.125326+/-0.000247` (+1.3 sigma), eqarea
`0.125007+/-0.000249` (+0.03 sigma) -- both consistent with the exact 1/8,
extending objective (2)'s existing confirmation to the largest lattice
this campaign has run.

**Real off-diagonal `kappa2_r` (the actual objective-1 gate, not a
diagonal proxy) at n_refine=512**:

```
l        naive        eqarea      naive/eqarea ratio
1     9.5e-06      7.4e-06        ~1.3 (noise level, both consistent w/ 0)
2     8.4e-06      5.9e-06        ~1.4 (same)
3     5.2e-05      1.2e-05        ~4.2
4     4.0e-05      1.2e-05        ~3.4
5     3.9e-05      1.2e-05        ~3.2
6     3.2e-05      1.8e-05        ~1.8
7     4.3e-05      1.9e-05        ~2.3
8     4.6e-05      1.6e-05        ~2.9
```

**Gate does NOT pass at n_refine=512 on either mesh mode.** l=1,2 remain
consistent with zero on both meshes (as at every smaller resolution
tested). l>=3 is still a real, nonzero plateau on both -- naive's l=3
value at 512 (5.2e-5) is actually the *largest* naive l=3 value in the
whole n_refine=4..512 ladder (previous high was ~4.7e-5 at n_refine=96),
continuing the "rises, does not shrink" pattern first seen at n_refine=128
(2026-08-24 entry) out to 4x that resolution. eqarea's l=3 at 512 (1.2e-5)
is within its own noise band relative to n_refine=384's 5.7e-6 and
128's 9.1e-6 -- flat, not shrinking either. The naive/eqarea suppression
factor at n_refine=512 (~3-4x for l=3-5) is smaller than the ~30-40x seen
at n_refine=96-112 in the 2026-08-29 100-shard entry -- **not necessarily
a real narrowing of the gap** given n_refine=512 has ~30x fewer
measurements/shard than the 100-shard ladder (32 shards x ~35k meas here
vs 100 shards x 100k there) and correspondingly larger error bars; not
enough statistics at 512 yet to say the mesh-mode gap itself is closing.

Plots: `campaign_runs/lean_ladder_512_2026-08-27/plots/
mesh_compare_{R_l,kappa2_r,kappa3_r,kappa4_r,ratio_all_metrics,scalar}_vs_invL.png`.

**Separately, confirmed the naive-vs-eqarea-vs-eqrp 3-way `kappa2_r`
comparison (job 7363395, flagged "not yet analyzed" in the 2026-08-29
entry) had actually already completed under a different job ID**:
7363395 itself died instantly (`.venv_plot/bin/activate` path typo'd as
`/var/spool/sge/.venv_plot/...` in the task script, `exit_status=1`,
0 wallclock) but a corrected resubmission, **job 7364152**, ran
successfully the same day and its output was sitting unread in
`campaign_runs/lean_mesh_compare_100shards_2026-08-27/three_mesh_compare_logs/
s2prec_compare_three_mesh.o7364152`. Result (n_refine up to 128, the
100-shard ladder): `equal_rp`'s l=3 `kappa2_r` is comparably small to
`equal_area`'s at every resolution (both ~1e-6 range by n_refine>=16,
vs. naive's ~2-5e-5) -- `equal_rp` does not obviously do better than
`equal_area` at this statistics level, contrary to the paper-motivated
hope that circumradius+perimeter equalization (not just area) would be
needed. Neither smoothed mesh's plateau shrinks toward zero across the
ladder. Full per-l table in the log file above; not yet turned into a
plot.

**Bottom line, updated**: objective (1) (spherical-symmetry / `kappa2_r`
gate) still does not pass, now confirmed at the largest mesh this
campaign has built (n_refine=512, ~2.6M sites) and across all three mesh
modes tried (naive, equal_area, equal_rp) at up to n_refine=128. No mesh
smoother tried so far eliminates the l>=3 signal or makes it shrink with
resolution -- it remains an open, unexplained finding.

## 2026-10-04 (new session, new machine) -- repo copied off the cluster
## into a standalone checkout; environment rebuilt from scratch; pipeline
## re-validated end to end; n_refine=768 ladder prepared

New user, new working copy. This checkout was copied out of the original
cluster project directory and is **missing every gitignored/bulk
artifact**: no `bin/`, no `.venv_plot/`, no `campaign_runs/` (so no
production data and none of the plots any prior entry refers to), and no
`.gitignore` itself. No `qsub`/`qstat`/`qacct`, no `module`, no
`/projectnb` -- this session had no cluster access at all, so nothing
below was run on SCC.

**Environment rebuilt and verified on a plain Linux box** (useful as a
record that the package really is self-contained, per the 2026-08-22
entry):

- `boost` from the distro (`libboost-math-dev`) rather than
  `module load boost/1.83.0`; Eigen already vendored.
- Both drivers compile clean with the `cluster/build.sh` flag set, minus
  the SCC-specific `-I$SCC_BOOST_INCLUDE`:
  `g++ -g -O3 -fopenmp -Wall -Wno-deprecated-declarations
  -Wno-sign-compare -I include -I include/unsupported
  -DGRP_DIR='"<root>/grp"' src/<driver>.cc -o bin/<driver>`
  (warnings are all pre-existing `fscanf` `warn_unused_result` noise from
  `lattice.h`).
- `.venv_plot` recreated: `python3 -m venv .venv_plot` +
  `pip install numpy scipy matplotlib sympy`. Note the original venv was
  Python 3.10; this one is 3.11 and every analysis script ran unmodified.

**Full `mesh_compare_kappa2_r_vs_invL.png` pipeline re-validated end to
end** on a deliberately tiny local ladder (`campaign_runs/
local_smoke_2026-10-04/`, naive+eqarea, n_refine={4,6,8,12,16}, 4
shards x 40000 traj -- 80k measurements/point vs. production's 1.1M).
Mechanism confirmed: `lean_harmonic_stats` -> `symmetry_check_lean_kappa2.
point_kappa2` -> `symmetry_test.analyze_blocks` -> `plot_mesh_compare_full.
py` produces all six figures including the target filename. **The numbers
in that smoke run are noise, not physics** -- kappa2_r ~1e-4 flat across
every l with naive/eqarea indistinguishable, exactly as expected at 1/14th
the statistics and 1/32nd the top resolution. Raw `.dat` deleted after
plotting (175MB); the plots are kept as the record.

**Real fix, not just setup: hardened `analyze_lean_harmonic_stats.load()`
against a field-boundary-truncated final block.** This is the "not yet
done" TODO carried since 2026-08-27 and the thing that forced the manual
`sed -i '$d'` on 3 eqarea shards in the 2026-09-02 entry. The existing
`try/except ValueError -> break` only catches a tear *mid-token*; a
walltime kill that cuts `fprintf` at a field boundary yields a line of
859/851/763 clean floats instead of 975, parses without raising, and
appends a short block whose missing `(l,m,mp)` entries surface much later
as a `KeyError` in `to_analyze_blocks_format`. Now the line width is
checked explicitly against `6 + n_lm + 2*n_offdiag` and a short final
block is dropped with a warning. Reproduced the exact 975->859 failure on
a smoke shard and confirmed the fix drops one block (200 -> 199) and
converts cleanly.

**Prepared `cluster/sge/submit_lean_ladder_768.sh`** (copied from
`submit_lean_ladder_512.sh` per the "copy, don't edit in place"
convention). Extends the ladder one rung to n_refine=768 (5,898,242
sites, 2.25x the sites of 512) with **every statistics/physics parameter
held identical to the 512 push** -- l_max=8, exact_sinh, n_traj=70000,
n_skip=2, n_wolff=5, n_metropolis=4, jack_block_size=1000, 32 shards --
so the new rung is comparable to the existing ladder rather than
confounded. Four deliberate infrastructure changes, each tied to a
documented failure:

1. **Production split into a short-rung and a long-rung array** (4..128
   at `h_rt=12:00:00`, 192..768 at `h_rt=120:00:00`) instead of one array
   under a single h_rt. The 512 push ran everything at `h_rt=30:00:00`
   and lost 19 eqarea tasks to walltime kills, all at its top rung.
2. **`h_rt=120:00:00` on the long array.** Cost model calibrated against
   this campaign's own ~9.4ms/measurement at n_refine=32 scaling linearly
   in n_sites (which correctly predicted 512's ~24h) puts 768 at ~53h/
   shard; 512's kills show that model runs optimistic, hence the margin.
3. **`-pe omp 4` on the long and prewarm arrays, purely for SCC per-core
   memory** (same rationale as `submit_gd_stress_1024.sh`). Measured here:
   1.96 kB/site for `lean_harmonic_stats`, 1.48 kB/site for
   `ising_s2_crit`, i.e. ~11.6GB and ~8.7GB at n_refine=768 -- past a
   default slot, which is why 512 (~5.1GB) got away without a reservation.
4. **Seed bases 920000/921000/930000/931000**, disjoint per arm and from
   the 512 push's 900000/910000.

Prewarm keeps the `-sync y` + independent-file-verification pattern (the
`-sync y` ban in CLAUDE.md is about *production* arrays, and the explicit
check, not the sync's return, is the real gate). `LADDER_SHORT=""
LADDER_LONG="768"` runs only the new rung if the 512 data is reachable.
**Dry-ran the whole script against a stubbed `qsub`**: normal path, the
only-new-rung override, the prewarm-failure abort (correctly never
reaches step 2), and the task-index -> (n_refine, shard, seed) layout for
the long arm all verified. **Not submitted** -- no cluster access from
this session.

Also added a `.gitignore` (`bin/`, `.venv_plot/`, `__pycache__/`,
`campaign_runs/**/*.dat`, `core.*`) since the copied checkout had none,
and untracked 19 stale `.pyc` files (cpython-36/310) that had been
committed.

## 2026-10-04 (cont): group-theory result -- l<=2 CANNOT show kappa2_r
## breaking, so half the step-6 gate is vacuous; plus stale-file cleanup

**New analytical result, `scripts/icosahedral_harmonic_decomposition.py`**
(committed as a script per directive.md's "derive it, don't quote it from
memory" rule; the character table is checked for orthonormality at runtime
so a typo can't propagate silently).

A q=5 mesh is exactly icosahedrally symmetric -- machine precision for
`naive` and `equal_area`. The action, quadrature weights and MC measure
all inherit that symmetry, so `M_l[m,m']` commutes with the
representation of the icosahedral group `I` on the (2l+1)-dim Y_lm space.
By Schur, `M_l` is a scalar on each irreducible component, so it has at
most as many distinct eigenvalues as there are inequivalent irreps in the
restriction. Decomposing:

```
l=1  (dim 3)  = T1          -> ONE irrep -> M_1 = c*Identity exactly
l=2  (dim 5)  = H           -> ONE irrep -> M_2 = c*Identity exactly
l=3  (dim 7)  = T2 + G      -> two irreps -> kappa2_r may be nonzero
l=4  (dim 9)  = G + H       -> two irreps
l=6  (dim 13) = A + T1 + G + H   (first l>0 containing the trivial rep)
```

**Therefore kappa2_r, kappa3_r and kappa4_r vanish IDENTICALLY at l=1 and
l=2, at any mesh resolution, however coarse.** Those two levels are not
capable of showing SO(3) breaking in this observable.

**Consequence: the "l<=2 clean, l>=3 plateau" split that this campaign has
reported since the 2026-08-24 gate evaluation is not a physical boundary.**
The two levels that pass PLAN.md step 6 are exactly the two that cannot
fail it, and l=3 is simply the first level where a nonzero kappa2_r is
kinematically allowed at all. Every measured l=1,2 value (9.5e-6 / 8.4e-6
at n_refine=512, "consistent with zero") is pure MC noise around an exact
zero, as it must be -- it is a useful null test of the pipeline, but it is
NOT evidence that symmetry is restored at low l.

**What this does NOT explain**: why the l>=3 plateau fails to *shrink*
under refinement. Icosahedral symmetry is exact at every n_refine, so this
argument constrains *which* l can break SO(3), not how fast the breaking
dies away. The open question is unchanged; its framing is not.

**Validation**: the same computation reproduces, analytically, the l=6
fact the campaign previously found only empirically via the deterministic
free-scalar FEM run (`src/fem_scalar_test.cc`, journal 2026-08-25) -- l=6
is the first l>0 whose restriction contains the trivial rep A, and the
script's invariant levels (6, 10, 12, ...) match PLAN.md's empirically
noted "l=6, 10, 12, 15, ..." exactly. That agreement is the cross-check
that the character table and decomposition are right.

**Latent hazard found while auditing the mesh cache**: `ReadPositions`
validates ONLY the `n_sites` header, and the cache filename
(`q5k<n>[_eqarea|_eqrp]_step<%.3f>.dat`) records the step size but **not
`--equal_area_iters`**. A cache written by a push using 1000 iterations is
therefore silently accepted by a later push asking for 20000 -- same
n_sites, same step, different relaxation state, no warning. Not hit in
practice (every current submit script prewarms its own job-local cache),
but it is a real silent-provenance hole in a repo whose whole SGE
convention is provenance-first. Worth putting `iters` in the cache key.

**Stale files deleted** (per user direction):
- `mesh_cache/q5k{2,3,4,6,8,12,16,24,32,48,64,96,128}_step0.300.dat` (13
  files) -- `equal_area` caches written by `submit_push7_dual.sh` under
  the OLD naming, before the `_eqarea`/`_eqrp` suffix was added. Current
  code looks for `q5k<n>_eqarea_step0.300.dat`, so these were unreachable
  dead weight. Deliberately NOT renamed into the new scheme: given the
  iters hole above, silently reviving a 6-week-old cache of uncertain
  relaxation provenance is exactly the wrong move, and regenerating is
  ~30s at n_refine=128. The 7 `_eqrp` files are still reachable under
  current naming and were left alone.
- `grp/.venv/` -- a committed Python 3.6 virtualenv, non-functional here
  (its interpreter isn't even executable in this checkout).
- `scripts/rebin_jackblocks` -- a committed compiled x86-64 ELF binary;
  `rebin_jackblocks.c` and `.awk` remain, so nothing is lost.

## 2026-10-04 (cont): cluster/build.sh now builds every driver, not just
## ising_s2_crit

`cluster/build.sh` built exactly one target (`ising_s2_crit`), and
`CLAUDE.md` told you to "compile by hand with the same flags, see each
file's header comment for its exact g++ invocation" for everything else --
but **9 of the 16 `src/*.cc` files have no such header comment**,
including `lean_harmonic_stats.cc`, the driver every current production
job actually runs. So the documented path to building the workhorse did
not exist, and a new user following CLAUDE.md could not submit a job.

Verified first that all 16 targets compile with one identical flag set --
there is nothing per-file about any of them -- then rewrote the script as
a loop over a target list:

```
cluster/build.sh                      # 6 production drivers (default)
cluster/build.sh --all                # + 10 diagnostics
cluster/build.sh lean_harmonic_stats  # named target(s)
cluster/build.sh --list
```

Default behaviour is a superset of the old script's (still builds
`ising_s2_crit`), so nothing that called it before changes.

Two portability fixes folded in: boost is loaded from the module system
only when one is present (`type module`), falling back to a system-wide
boost otherwise, and `-I "$SCC_BOOST_INCLUDE"` is omitted when that
variable is unset. The same script now works unchanged on the SCC and
off-cluster -- the off-cluster path is what this session used to build and
validate the whole pipeline. Per-target build logs land in
`bin/.<target>.buildlog` and are removed on success; a bad target name or
a compile failure exits nonzero naming the target.

CLAUDE.md's "Build" section rewritten to match, and corrected on two
points: that the diagnostics are safe to re-run (verified: every `test_*`
is pure stdout with zero file writes; `fem_scalar_test` touches only a
mesh cache, `diag_optimize_integrator`/`dump_face_sa_range` only a path
you pass them), and that a header edit does NOT trigger a rebuild since
there is no dependency tracking.

## 2026-10-04 (cont): commutant cross-check added, and the q=3/4/5
## comparison -- why the icosahedron buys an extra clean harmonic level

Extended `scripts/icosahedral_harmonic_decomposition.py` with a second,
**character-table-free** computation of the same result, as a guard
against the stored table being subtly wrong. The commutant of the spin-l
rep -- the space of all operators commuting with every group element,
which is exactly where `M_l` is forced to live -- has dimension
`sum_a n_a^2 = (1/|G|) sum_g |chi_l(g)|^2`, needing only the group's class
sizes and rotation angles. `M_l` is forced proportional to the identity
iff that dimension is 1. The script now asserts this agrees with the
character-table decomposition for every l, and it does.

Also added the `q` comparison, which answers "why q=5" quantitatively:

```
 q  group                |G|   first l where kappa2_r is symmetry-allowed
 3  T  (tetrahedron)      12   l = 2
 4  O  (octahedron)       24   l = 2
 5  I  (icosahedron)      60   l = 3
```

The icosahedral group is the largest of the three admissible base
polyhedra (`S2.h`: q=3/4/5 = tetrahedron/octahedron/icosahedron), so it
keeps the harmonics degenerate longest: q=5 forces BOTH l=1 and l=2 to be
exactly scalar, where q=4 and q=3 force only l=1. Concretely, on an
octahedral mesh `l=2` would already split (into E + T2) and would show a
nonzero kappa2_r plateau of its own -- so the campaign's "l=1,2 clean"
observation is specific to the icosahedral mesh family, not a general
property of the method. Worth knowing if a q=4 cross-check is ever run as
a control: its gate would have one fewer null level, and an l=2 signal
there would be expected rather than alarming.

## 2026-10-04 (cont): convention-free numerical proof of the l<=2 result;
## job_progress.sh; quadrature-weight check

**`scripts/verify_harmonic_rep_numerically.py`** -- independent
confirmation of the previous entry's Schur argument that uses no
character table, no Wigner D matrix, and no Y_lm phase convention, and
reads the campaign's OWN group data (`grp/elem/o3q5.dat`, the same file
the simulator reads). For each of the 60 proper icosahedral rotations it
*discovers* the matrix `T(R)` acting on the `Y_lm` by least squares
rather than asserting it, then checks the fit residual, unitarity, the
representation property, and the commutant dimension. Results at 400
sample points:

```
 l   max fit resid   max |TT^H - I|   max rep error   commutant dim
 0      2.8e-16         2.0e-15         1.0e-15             1   <- forced scalar
 1      3.6e-15         1.8e-15         1.6e-15             1   <- forced scalar
 2      8.0e-15         2.0e-15         1.2e-15             1   <- forced scalar
 3      1.3e-14         2.7e-15         1.3e-15             2
 4      1.9e-14         1.8e-15         1.3e-15             2
 5      2.5e-14         5.3e-15         3.2e-15             3
 6      3.2e-14         5.2e-15         4.5e-15             4
 7      3.9e-14         3.3e-15         2.3e-15             4
 8      4.6e-14         5.6e-15         3.7e-15             6
```

Every commutant dimension matches `sum_a n_a^2` from the character-table
decomposition exactly (l=8's 6 = 1+1+2^2 from `T2 + G + 2H`). The
`max fit resid` column is the one that matters conceptually: it *measures*
rotation-invariance of `span{Y_lm}` at fixed l rather than assuming it --
had the span not been closed, no `T` could fit and the residual would be
O(1).

**`cluster/sge/job_progress.sh`** -- `qstat` answers "what is the
scheduler doing", which is usually not the question; this also counts
shard files actually on disk per `(n_refine, mesh_mode)` against the
expected shard count, and re-checks every shard for a ragged final line
(the walltime-kill truncation of journal.md 2026-09-02). Tested against a
synthetic tree including a deliberately truncated file.

**Quadrature-weight check (throwaway, not committed)** -- measured what
`UpdateWeights()` actually produces and what it buys, since the campaign
has already been bitten once by a normalization bug in this area:

- `sum_i w_i = n_sites` **exactly** (642.0000000000 at n_refine=8,
  10242.0000000001 at 32), confirming CLAUDE.md's corrected statement and
  *not* `4*pi`. So `w_i` is dimensionless with mean 1 -- a relative area
  share, converted to a true solid angle by `*4*pi/n_sites`.
- `w_i` spans a factor of ~1.8 (n_refine=8) to ~2.0 (32) between its
  smallest and largest value, and that spread does NOT shrink with
  refinement -- the mesh stays genuinely non-uniform, which is why a
  weight is needed at all.
- Discrete orthogonality `sum_i w_i (4pi/n_sites) Y_lm Y*_l'm'` vs the
  same sum with the weights dropped:

```
  pair             weighted (k=8)  weighted (k=32)   UNweighted (k=8)  UNweighted (k=32)
  (3,1)-(3,1)        -5.1e-03        -3.3e-04          -5.7e-02          -5.6e-02
  (2,1)-(4,1)         5.0e-03         3.2e-04           5.6e-02           5.5e-02
  (0,0)-(6,0)         8.1e-03         5.2e-04           9.1e-02           8.9e-02
```

With weights the error falls by ~16x from n_refine 8->32, i.e. **O(a^2)**,
exactly as a convergent quadrature should. **Without weights it does not
improve at all** (5.7e-2 -> 5.6e-2) -- the unweighted sum converges to the
wrong answer, so the weights are not a refinement, they are what makes the
projection converge in the first place.
- Side observation, consistent with the group theory above: `(1,0)-(1,0)`,
  `(2,0)-(2,0)` and `(1,0)-(3,0)` come out exact to machine precision even
  UNWEIGHTED, while `(0,0)-(6,0)` does not -- the same icosahedral-
  invariance pattern (low l protected, l=6 the first leak) showing up in
  the quadrature rather than in `M_l`.

## 2026-10-04 (cont): measured the equal_area prewarm cost; added a
## SKIP_PREWARM recovery path after `-sync y` blocked a live submission

First real submission of `submit_lean_ladder_768.sh` (by the user, on
SCC) sat in step 1 with 3/16 prewarm tasks running. Expected -- `qsub
-sync y` blocks by design -- but it exposed a design miss: **`-sync y`
was the wrong choice for THIS push.** Every earlier push topped out at
n_refine<=128, where the equal_area relaxation is ~10 min, so blocking
cost nothing. At n_refine=768 the relaxation is 36x the work of 128. The
house prewarm pattern was copied from `submit_lean_ladder_512.sh` without
re-deriving its cost at the new top rung.

**Measured relaxation cost** (`EqualizeFaceAreas`, 20000-iter cap,
step=0.3), locally, rather than extrapolating CLAUDE.md's note:

```
n_refine   wall      iters_used
  16       0.26 s    (early stop, well under cap)
  32       8.2 s     9908   (early stop)
  48      37.0 s
  64      68.5 s     20000  (hits cap)
  96     157.5 s     20000
```

Two regimes. Below ~n=48 the rel_tol/patience early stop fires, so cost
grows faster than n^2 (the iteration count is still rising). **From n=64
up the cap is always hit, so iterations are constant and cost is clean
n^2**: 157.5/68.5 = 2.30 vs (96/64)^2 = 2.25. An earlier read of the
32-vs-64 numbers looked like n^3 -- that was the iteration count doubling
across the early-stop boundary, not a real cubic.

Extrapolating the capped regime and anchoring to CLAUDE.md's cluster
datapoint (~30s at n=128 for 1000 iters => 600s for 20000, vs 274s
locally, so the SCC core is ~2.2x slower here):

```
  n_refine    192     256     384     512     768
  SCC prewarm  23m     40m    1.5h    2.7h    6.0h
```

So the prewarm array's long pole is ~6h, comfortably inside its
`h_rt=24:00:00`. No walltime risk -- but ~6h is far too long to hold a
bare login shell.

**Fixes, all tested against a stubbed qsub and a synthetic tree:**

1. `SKIP_PREWARM=1` skips the prewarm qsub but still runs the cache
   verification (which was always the real gate, not the `-sync` return).
   This is the recovery path if the submitting shell dies: the prewarm
   array survives -- it belongs to SGE, not the shell -- but the
   production arrays never get submitted. Re-running the script as-is
   would submit a SECOND prewarm array racing the first on the same cache
   files, the exact race the prewarm exists to prevent.
2. **`OUT_ROOT_BASE` is now overridable.** This was a real bug in the
   recovery path: the directory is stamped `$(date +%F)`, so a
   SKIP_PREWARM re-run the next day -- the normal case, given a 6h
   prewarm that can cross midnight -- would compute a fresh dated
   directory, find an empty mesh_cache, and refuse.
3. The manifest is no longer overwritten on a re-run; it is this push's
   provenance record and a rerun's sha/timestamp would silently replace
   the ones that actually produced the prewarmed meshes.
4. `job_progress.sh` now reports prewarm progress (cached ladder points
   vs. still relaxing, read from the manifest's own ladder). Without it
   the tool says "no output yet" for six hours and looks stuck -- which
   is exactly what the user saw.

**Still open for a future push**: `-hold_jid <prewarm_jobid>` on the
production arrays would let the submit script return in seconds instead
of blocking at all. Not adopted here because `-hold_jid` releases on
prewarm *completion* regardless of exit status, so it would need the
production task script to handle a missing cache safely first.

## 2026-10-04 (cont): mesh-cache verification now checks exact file size

The prewarm verification used `[[ -s "$f" ]]` -- "exists and non-empty".
That is not enough. A prewarm task killed mid-`fwrite` (walltime limit,
or a `qdel`) leaves a TRUNCATED cache file, which passes `-s` but is
correctly rejected by `ReadPositions` at runtime as a cache miss -- at
which point all 32 shards of that ladder point relax the mesh from cold
simultaneously and race to rewrite the same file. That is exactly the
failure the prewarm stage exists to prevent, surfacing in production
instead of at submit time.

`WritePositions` writes an int32 `n_sites` header plus `n_sites * Vec3`
(3 doubles), so a complete q=5 cache is exactly

    4 + 24*(10*n^2 + 2)  =  240*n^2 + 52  bytes

verified against real cache files at n_refine=4/8/16 (3892 / 15412 /
61492 bytes, exact). The submit script now checks that size per ladder
point and refuses to submit production otherwise, naming the offending
point and the byte delta. Tested by truncating a synthetic n_refine=768
cache by 100 bytes out of 141 MB -- caught, production not submitted.

Relevant immediately: the user is `qdel`-ing the first 768 submission's
prewarm array mid-flight, which can leave exactly such a partial file.

## 2026-10-05: SKIP_PREWARM no longer needs a filepath passed by hand

Per user preference ("without injecting any filepaths to account for the
change of day"), `submit_lean_ladder_768.sh` now resolves its output
directory asymmetrically:

- a NEW run (no `SKIP_PREWARM`) is stamped with today's date, as before;
- a `SKIP_PREWARM=1` recovery run adopts the **most recent existing**
  `campaign_runs/lean_ladder_768_*` directory instead of today's date,
  and refuses with a clear message if none exists rather than inventing
  an empty one;
- `OUT_ROOT_BASE` still overrides both, for resuming a specific older run
  when several are on disk.

This closes the day-boundary trap recorded in the 2026-10-04 entry: the
prewarm runs for hours, so recovery is normally attempted the next day,
and a plain date stamp would then invent a fresh directory, find an empty
mesh_cache, and refuse purely because of the clock. All three paths
tested against a stubbed qsub (fresh run, recovery adopting a
past-dated dir with no env var, and refusal with nothing on disk ->
exit 1).

Also removed `campaign_runs/lean_ladder_768_2026-10-04/manifest.json`,
which was committed by accident in 7dddcf0. It was generated by a
stubbed-qsub dry run on this container, not by a real submission, and its
`simulator_sha256` is this container's local build -- i.e. it was
misleading provenance for a run that never happened. The first real
submission writes its own manifest.
