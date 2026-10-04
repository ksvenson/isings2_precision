# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A high-precision cluster campaign to verify spherical symmetry and extract
CFT data for the 2D Ising model simulated directly on S^2 (as a simplicial
lattice), by decomposing spin two-point functions into spherical harmonics.
`PLAN.md` is the scientific plan (phases, method, statistics policy, gates)
— read it first. `directive.md` holds standing cross-cutting conventions.
`journal.md` is the full running log (session-by-session findings, dead
ends, bug fixes, job IDs) — check its tail for the current state before
starting work; this file only covers what's stable enough to guide new
code, not a history of how it got that way.

**Current state of the two campaign objectives** (see `PLAN.md` for full
detail/status, `journal.md` for how each conclusion was reached):

- **Objective 2 (extract Delta_sigma, Delta_epsilon): Delta_sigma = 1/8
  CONFIRMED**, via two independent methods that agree to <1.2 sigma once
  correctly normalized/extrapolated — the direct real-space fit (`Delta_inf
  = 0.125` region, "Real-space Delta_sigma extraction" below) and the
  l-space `Trace(M_l)` method (`lean_harmonic_stats`) once a `(2l+1)`
  normalization bug was found and fixed. Delta_epsilon is not yet
  addressed (no energy-operator definition built).
- **Objective 1 (verify spherical/SO(3) symmetry via `R_l`/`kappa2_r`):
  still open.** `l<=2` are clean (consistent with zero, shrink or stay flat
  under refinement). **`l>=3` shows a real, many-sigma, resolution-
  independent `kappa2_r` plateau** that does not shrink out to the largest
  meshes tested (n_refine up to 128 with the full off-diagonal `M_l`
  matrix; up to 512 in flight for a diagonal-only proxy). `equal_area`
  mesh smoothing reduces this plateau by ~30-40x vs. the naive mesh but
  does not eliminate it or make it shrink with resolution; `equal_rp`
  (circumradius+perimeter smoothing, motivated by arXiv:2407.00459) is the
  next candidate, comparison job submitted but not yet analyzed as of the
  last journal entry. Bulk magnetization moments and diagonal-only per-l
  `Delta_l_pair`/`delta_l` show no naive-vs-equal_area discrepancy — the
  symmetry-breaking signal lives specifically in the off-diagonal `M_l`
  channel, not in these observables (a useful negative control).

## Relationship to `IsingS2/`

**Fully self-contained as of 2026-08-22** — everything needed to build and
run the simulator lives under this directory: `src/`, `include/` (headers +
vendored Eigen), `grp/` (icosahedral group-theory data), `scripts/`
(analysis). `../IsingS2/` remains the *origin* of this code (an older,
unstructured research repo) but is no longer read by any build or run path
here, and **does not auto-track** — a fix landing there must be manually
re-copied into `include/`/`src/`/`grp/` here and logged in `journal.md`.

- **`src/ising_s2_crit.cc`** is the driver this campaign actually uses (not
  `isingS2.cc` — stale, uses the base `QfeLattice` type instead of
  `QfeLatticeS2`, missing methods the S2 subclass needs). Builds a
  `QfeLatticeS2` mesh, runs Metropolis+Wolff updates, and measures the
  spherical-harmonic projections this campaign's method needs.
- `include/S2.h` (`QfeLatticeS2`) — mesh construction, quadrature weights
  (`UpdateWeights()`, the cotangent-Laplacian/lumped-Voronoi-dual-area FEM
  scheme normalized to `sum wt_i = n_sites`, **not** `4*pi` — a
  historical doc bug in this file, now corrected; `ising_s2_crit.cc` uses
  `UpdateWeights()`, not `OptimizeIntegrator(l_max)`, despite the latter
  existing), `UpdateYlm(l_max)`/`GetYlm(s,l,m)` (cached `Y_lm`, used via
  `CalcYlm` — returns `Y_lm` not `Y*_lm`, l-independent so it cancels in
  every ratio-based diagnostic), and the mesh-smoothing routines
  `EqualizeFaceAreas`/`EqualizeCircumPerim` (see "Mesh construction modes"
  below).
- `include/ising.h` (`QfeIsing`) — spin field / update code
  (`HotStart`/`ColdStart`, `Metropolis`, `WolffUpdate`,
  `ReadField`/`WriteField`, `Action`).

## Build

```
cluster/build.sh
```

Builds `src/ising_s2_crit.cc` into `bin/ising_s2_crit` (gitignored, rebuild
after any change to that file). `module load boost/1.83.0` (for
`boost/math/special_functions/spherical_harmonic`), Eigen vendored
header-only at `include/Eigen` (no module needed), system `g++ 8.5.0` (no
`module load gcc`), `-DGRP_DIR` pointed at this package's own `grp/`. The
original `IsingS2/src/Makefile` is Homebrew/macOS-flavored and does not
apply here.

Other `src/*.cc` binaries (`lean_harmonic_stats`, `save_configs`,
`full_corr_test`, `real_space_2pt_test`, `analyze_north_pole_configs`, and
the ad hoc `test_*.cc` diagnostics) are not in `cluster/build.sh` — compile
by hand with the same flags, see each file's header comment for its exact
g++ invocation.

## Running the measurement driver

```
bin/ising_s2_crit --q 5 --n_refine 2 --l_max 4 \
  --coupling_rule exact_sinh --data_dir campaign_runs/<name> \
  --n_therm 2000 --n_traj 20000 --seed 1234
```

Key flags (`getopt_long` table near the top of `main`): `--q`/`--n_refine`
(icosahedral mesh family + refine level → `run_id = q<q>k<n_refine>`,
`n_sites = 10*n_refine^2 + 2`), `--l_max` (harmonic cutoff),
`--coupling_rule {duality,exact_sinh}`, `--mesh_mode
{naive,equal_area,equal_rp}`, `--cold_start`,
`--n_therm`/`--n_traj`/`--n_skip`/`--n_wolff`/`--n_metropolis`,
`--data_dir`, `--seed`, `--wall_time`.

**The driver does not `mkdir -p` its output directory** — `<data_dir>/<run_id>/`
must already exist before invoking it (SGE task scripts must create it).

### Coupling rules (`--coupling_rule`)

- `duality` (default) — pre-existing dual-link-length formula,
  `K = 0.5*asinh(1/sinh(2L))`, `L` from a cross-ratio of face-angle
  cosines/edge lengths around the link.
- `exact_sinh` — Twist's flat-triangle star-triangle rule
  `K(theta) = 0.5*asinh(cot(theta))`, applied per adjacent face and
  averaged across the link's two faces (a judgment call, since the S2 mesh
  has no single-triangle-per-edge convention). **This is the coupling rule
  production runs use** — the first precision pass used `exact_sinh`
  per 2026-08-21 user direction; `duality` remains a free, non-blocking
  comparison.

### Mesh construction modes (`--mesh_mode`, default `naive`)

Motivated by the objective-1 finding that `kappa2_r` for `l>=3` plateaus at
a nonzero value on the naive mesh instead of shrinking with refinement —
naive-mesh triangle-area non-uniformity is a live suspect (though not
confirmed sufficient, see "Current state" above):

- `naive` — unmodified constructor mesh (flat-subdivide-then-radially-
  project). Exactly icosahedrally symmetric to machine precision.
- `equal_area` — gradient descent + backtracking line search minimizing
  `E = sum_f (FlatArea(f) - mean)^2` (`QfeLatticeS2::EqualizeFaceAreas`,
  `--equal_area_iters`/`--equal_area_step`, default 1000 iters). True
  per-vertex triangle-area gradient, projected onto the tangent plane and
  renormalized to the unit sphere each step. Icosahedral-symmetry
  preservation verified to ~1e-15 (`src/test_gd_symmetry.cc`). Negligible
  cost (~30s at n_refine=128, login node, single-threaded, for 1000
  iters). Reduces the `l>=3` `kappa2_r` plateau by ~30-40x vs. naive but
  does not eliminate it or make it shrink with resolution.
- `equal_rp` — circumradius+perimeter mesh smoother implementing
  arXiv:2407.00459 (Brower & Owen) Appendix B.1/B.2: their coupling rule
  is only exactly critical when every triangle has equal circumradius
  **and** equal perimeter, not equal area (their own "basic" mesh and an
  equal-area smoother both fail for that reason — this is what motivated
  building this mode). `QfeLatticeS2::EqualizeCircumPerim` minimizes the
  joint objective `E = E_R + E_P + E_A` (area kept in the objective
  despite the paper only requiring `E_R+E_P`, since area equalization is
  an independently confirmed win). Gradient via exact local finite
  differences (`LocalEnergyWithVertexMoved`), not a hand-derived closed
  form. Symmetry preservation good but not machine-precision (max
  mismatch 7.6e-9 to 5.6e-6 across n_refine 4-32) due to FD noise. Reuses
  `--equal_area_iters`/`--equal_area_step` (no separate flags), cache
  suffix `_eqrp`. Not yet compared against `equal_area` at production
  scale for the `kappa2_r` gate (in flight, see "Current state" above).

Two modes were tried and dropped from the CLI (code remains in `S2.h` for
reference, not exposed via `--mesh_mode`): `equal_area_analytic`
(Snyder-style analytic construction — exact symmetry but worse area
uniformity than `equal_area`) and `--orbit_path`/`ReadOrbits` (removed as
needlessly complicated).

### Output files (`<data_dir>/<run_id>/<run_id>_*_<seed_hex>.dat`)

- `*_bulk_*.dat` — action/magnetization/mag^2/U4/susceptibility.
- `*_ylm_2pt_*.dat` — diagonal-only `|S_lm|^2/vol^2` per `(l, 0<=m<=l)`.
- `*_ylm_2pt_full_*.dat` — full upper-triangle (`m<=m'`) of the Hermitian
  `(2l+1)x(2l+1)` matrix `M_l[m,m'] = <S_lm S*_lm'>` the symmetry test
  needs (negative-`m` reconstructed via `(-1)^m * conj()`).
- `*_legendre_2pt_*.dat` — an m-summed `C_l` proxy; not confirmed to equal
  `C_l` from the addition theorem, not used by any current analysis path.

## Other drivers (`src/`)

- **`lean_harmonic_stats.cc`** (`bin/lean_harmonic_stats`) — the current
  workhorse for both objectives at scale. Online-accumulating driver
  (during the same Metropolis+Wolff sweep as `ising_s2_crit`, no raw
  config or full `M_l` matrix ever written) computing `|S_lm|^2` (all m)
  and, as of the second build pass, the 444 off-diagonal `S_lm*conj(S_lmp)`
  terms too — i.e. it now has everything needed for the real `kappa2_r`
  gate, not just a diagonal proxy. Also accumulates `M`, `M^2`, `M^3`,
  `M^4` per jackknife block. Output size `O(n_blocks*n_lm)`, independent of
  `n_sites` — this is what makes n_refine>=128 ladders affordable. Uses
  `UpdateWeights` (not `OptimizeIntegrator`, sidesteps a negative-weight
  bug in the latter at n_refine>=32). Supports `--mesh_mode` (ported from
  `ising_s2_crit.cc`). Analyzed via `scripts/analyze_lean_harmonic_stats.py`
  (moments/`Delta_l_pair`/`Trace(M_l)`-based `Delta_s`) and
  `scripts/symmetry_check_lean_kappa2.py` (real off-diagonal `kappa2_r`,
  reuses `symmetry_test.py`'s `analyze_blocks` via a format converter —
  the diagonal-only `symmetry_check_lean.py` is superseded for new data,
  still needed to read old diagonal-only files).
- **`save_configs.cc`** (`bin/save_configs`) — dumps raw site positions +
  bit-packed spin configurations for offline re-analysis, `O(n_meas*n_sites)`
  bits (storage-scalable, the recommended path for reconstructing `M_l`
  offline at large `n_refine`). Paired with
  **`scripts/build_ylm_matrix_from_configs.py`** (offline `Y_lm`
  contraction: `M_l = S^H @ S / n_meas`, `S = spins @ A`). For large
  `n_refine` the per-shard matmul is compute-bound (not memory-bound after
  a real-vs-complex-matmul fix) — run via
  `cluster/sge/compute_ml_block_shard.py` (one shard → one pickled block)
  + `scripts/reduce_ml_blocks.py` (combine block pickles), not a serial
  per-point loop (the latter hit `h_rt` walltime kills at n_refine>=64).
- **`full_corr_test.cc`** (`bin/full_corr_test`) — accumulates the full
  dense `n_sites x n_sites` `<s_i s_j>` matrix, `O(n_sites^2)` memory
  (~840MB at n_refine=32, ~215GB at n_refine=128 — do not run at large
  `n_refine`). Paired with `scripts/build_ylm_matrix_from_corr.py`.
  Mathematically identical `M_l` estimator to the `save_configs` path (by
  linearity of expectation), just a different computational order/memory
  tradeoff. Useful only at small-to-moderate `n_refine`.
- **`real_space_2pt_test.cc`** (`bin/real_space_2pt_test`) — the primary
  Delta_sigma-extraction driver, see "Real-space Delta_sigma extraction"
  below.
- **`analyze_north_pole_configs.cc`** (`bin/analyze_north_pole_configs`) —
  offline, single-reference-row `<s_i s_j>` from `save_configs` shards
  (`O(n_sites*|G|)` setup, no `O(n_sites^2)` anywhere). Used by
  `scripts/fit_north_pole_ladder.py`.
- **Ad hoc diagnostics** (`test_ico_symmetry.cc`, `test_gd_symmetry.cc`,
  `test_gd_equal_area.cc`, `test_gd_convergence.cc`, `test_gd_stress.cc`,
  `test_eqrp_symmetry.cc`, `test_orbit_avg.cc`, `fem_scalar_test.cc`,
  `diag_optimize_integrator.cc`, `dump_face_sa_range.cc`) — one-off
  correctness/symmetry/performance checks, not production paths. See each
  file's header comment for build/run instructions.

## Real-space Delta_sigma extraction

The primary Delta-extraction method (see "Current state" above for why —
every l-space estimator built from `S_lm`/Legendre recursion integrates
over the *whole* sphere, so short-distance lattice contamination leaks
into every coefficient and cannot be excluded by any ratio built from
them; a windowed Legendre transform was tried three ways and fails for a
structural mathematical reason — see `journal.md`'s "SETTLES the
question" entry. The l-space `Trace(M_l)` method was later found to also
give the correct answer once a normalization bug was fixed — both methods
now agree; the real-space method remains the one to use for new
Delta-extraction work since it doesn't depend on getting that
normalization right).

- **`src/real_space_2pt_test.cc`**, own `bin/real_space_2pt_test`. Reuses
  the coupling-assignment/mesh-mode loop verbatim from
  `ising_s2_crit.cc`. For each of `--n_ref` random reference sites,
  targets are partitioned into exact equivalence classes under `Stab(i)`
  (the icosahedral-group subgroup fixing site `i`), each class's full
  120-element `G`-orbit is averaged per configuration and fed into an
  online mean/variance accumulator (`include/statistics.h`'s
  `QfeMeasReal`) — no raw-sample buffering, no arbitrary binning window
  (optional `--theta_min`/`--theta_max` only for compute-cost control).
  Output: `<data_dir>/<run_id>_real_space_2pt_<seed_hex>.dat`, one line
  per (reference site, class): `ref_site target_site theta mean err n
  orbit_size`.
- **`scripts/fit_real_space_2pt.py`** — fits `<s(theta)>` to
  `A*(2-2cos(theta))^{-Delta}` over a window, `err`-weighted, single-run
  (no jackknife) — quick spot checks. `theta_min>=0.3` excludes the
  contaminated short-distance region.
- **`scripts/fit_real_space_ladder.py`** — production continuum
  extrapolation: each independent-seed shard is its own delete-one-shard
  jackknife unit (no in-C++ jackknife accumulator for this driver), fits
  `Delta_inf + c*n_refine^-p` across the ladder. **`--power` defaults to 2
  but this is not validated for this mesh/coupling — always also try
  `--power 1`.** Production ladders were ~4.4 sigma off 1/8 at `p=2`
  (visibly bad fit) and <1 sigma at `p=1` — a linear leading
  discretization correction, not quadratic (no symmetry forbids an O(a)
  term here). This `p=1`-preferred finding is specific to `l<=1`/the
  real-space aggregate — the analogous l-space per-l study
  (`fit_lean_ladder_family.py`) found `l=1` alone prefers `p~1` but every
  `l>=2` decisively prefers `p~2` — do not assume one power applies to
  every harmonic-level diagnostic (e.g. the `kappa2_r` gate).
- **`scripts/plot_real_space_ladder.py`** — overlays `Delta_sigma(a)` vs
  `1/n_refine` plus extrapolation fits and the exact-1/8 line.
- **`cluster/sge/real_space_shard_task.sh`** / `submit_real_space_ladder.sh`
  — provenance-first per-`(n_refine,shard)` SGE array, same pattern as
  `production_shard_task.sh`.

## Python environment

```
.venv_plot/bin/python scripts/<script>.py ...
```

Local venv (gitignored, Python 3.10.12, `numpy`/`scipy`/`matplotlib`/
`sympy`, no `requirements.txt` — install ad hoc with `.venv_plot/bin/pip`).
SGE tasks that run analysis scripts `source .venv_plot/bin/activate`
rather than invoking the interpreter path directly. There is no test suite
or linter — correctness is checked via the ad hoc diagnostic binaries
(`test_ico_symmetry`, `test_gd_symmetry`, etc.) and validated production
runs, not an automated suite.

## Cluster conventions

SGE (`qsub`/`qstat`/`qacct`), not Slurm — see the workspace-wide
convention in `/projectnb/qfe/misra/CLAUDE.md`. Prefer detached `qsub`
over `qsub -sync y` for production arrays — `-sync y` has twice failed in
this campaign in ways that are hard to detect (once silently breaking
downstream script logic right after a successful prewarm, once hanging
10+ hours after the job it waited on had already finished server-side).
**After any `-sync y` step, verify the actual output landed with an
independent check (`ls`/`qacct`), don't trust the wrapping script's own
internal verification.**

`cluster/sge/` holds one array task per `(n_refine, shard_index)`
(provenance-first, task identity fixed by the tuple — see
`/projectnb/qfe/misra/CLAUDE.md`'s cross-project convention):

- `production_shard_task.sh` — the per-task `ising_s2_crit` invocation,
  fully parameterized via env vars set by `qsub -v` from the submit
  script (`ROOT LADDER_CSV N_SHARDS OUT_ROOT L_MAX COUPLING_RULE N_THERM
  N_TRAJ N_SKIP N_WOLFF N_METROPOLIS SEED_BASE JACK_BLOCK_SIZE`, plus
  optional `MESH_MODE`/`EQUAL_AREA_ITERS`/`EQUAL_AREA_STEP`/
  `MESH_CACHE_DIR`).
- `lean_harmonic_stats_task.sh` — the analogous per-task invocation for
  `lean_harmonic_stats`, shard-aware `(n_refine, shard_index)`.
- `prewarm_mesh_cache_task.sh` — one array task per ladder point,
  pre-populates a mesh-mode's position cache (`--mesh_cache_dir` in
  `S2.h`) serially before a production array with many shards/point runs
  against the same `(n_refine, mesh_mode)` — avoids concurrent shards
  racing to write the same cache file cold.
- `submit_*.sh` — one submit script per production push, each hardcoding
  that push's ladder/`l_max`/stats/mesh_mode/`jack_block_size`/`h_rt` and
  writing a `manifest.json` (ladder, stats, coupling rule, mesh mode,
  simulator sha256) into `campaign_runs/<job_tag>/` before submitting.
  Non-`naive` pushes run the prewarm array first (`qsub -sync y`) and
  verify every cache file exists before submitting production. **Copy an
  existing `submit_<name>.sh` as the template for a new push rather than
  editing one in place** — keeps the config that produced a given
  `campaign_runs/<job_tag>/` reconstructible from git history.

**`jack_block_size` is not a fixed constant** — pick it so the native
block count (`n_traj / jack_block_size`) lands in the few-hundred range
regardless of `n_traj`, targeting `PLAN.md`'s stats gate `n_blocks >=
10*(2*l_max+1)`. Recompute per push, don't reuse the last push's literal
value — too small a block size at large `n_traj` reproduces a multi-GB
jackblocks-file problem this campaign has hit before.

## Method summary (spherical-symmetry objective)

See `PLAN.md` for the full method and open items. In brief: for each mesh
resolution, generate independent thermalized configurations, project the
spin field onto spherical harmonics per configuration (`S_lm`), form the
`(2l+1)x(2l+1)` correlation matrix `M_l[m,m'] = <S_lm S*_lm'>` at each
harmonic level `l`, and check it is diagonal/degenerate in `m` to
statistical precision via the basis-independent reduced spectral cumulants
`kappa2_r, kappa3_r, kappa4_r` (`scripts/symmetry_test.py`'s
`analyze_blocks`, built from `Tr(M_l^k)` for `k=1..4`) — unlike
`R_l`/diagonal-chi2, these don't depend on the mesh's arbitrary
quantization axis. **Never hardcode a from-memory closed-form formula**
(e.g. the harmonic-level eigenvalue `C_l(l, Delta)`) without it first
being symbolically derived/verified against known special cases (see
`directive.md`; `scripts/derive_cl_closed_form.py`).

`q=5` meshes have exact icosahedral (not full SO(3)) symmetry with its own
invariant subspace starting at `l=6` (multiplicity 1) — `l=6` (and, at
higher `l_max`, `l=10,12,15,...`) should be excluded from any "does
SO(3)-breaking shrink under refinement" claim, confirmed via a
deterministic free-scalar FEM test (`src/fem_scalar_test.cc`) where `l=6`'s
off-diagonal `R_l` never shrinks under refinement/relaxation unlike its
neighbors. This is a separate, understood effect from the `l>=3` `kappa2_r`
plateau (l=6 is not an outlier in that diagnostic).
