#!/usr/bin/env bash
# One array task = one (n_refine, shard_index) real-space two-point-function
# measurement (real_space_2pt_test), for the Delta(a) continuum-
# extrapolation ladder (PLAN.md "next step" after the 2026-08-25
# breakthrough finding that real-space fits, not F_l/Delta_l_pair, are the
# correct way to extract Delta_sigma). Mirrors production_shard_task.sh's
# task-index layout exactly, but drives real_space_2pt_test instead of
# ising_s2_crit -- no l_max/jack_block_size (that binary has neither; it
# has no in-C++ jackknife accumulator). Each shard is a fully independent
# seeded MC chain; error bars across shards are computed offline in
# scripts/fit_real_space_ladder.py via a delete-one-shard jackknife, not by
# anything in this script.
#
# naive mesh_mode only (matches the two validated login-node spot checks
# in journal.md's "BREAKTHROUGH" entry) -- no mesh-cache prewarm step
# needed.
#
# Required environment (set by submit script via `qsub -v`):
#   ROOT LADDER_CSV N_SHARDS OUT_ROOT COUPLING_RULE
#   N_THERM N_TRAJ N_SKIP N_WOLFF N_METROPOLIS N_REF SEED_BASE
# Optional: THETA_MIN/THETA_MAX (default full sphere 0/pi) to restrict which
# reference-target pairs get built into equivalence classes -- see
# real_space_2pt_test.cc's header comment (2026-08-26 rewrite dropped
# --n_bins/--orbit_avg entirely in favor of exact icosahedral-symmetry
# equivalence classes; no binning parameter exists on the driver any more).
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
GLOBAL_IDX=$((SGE_TASK_ID - 1))
SHARD_IDX=$((GLOBAL_IDX % N_SHARDS))
L_IDX=$((GLOBAL_IDX / N_SHARDS))
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$L_IDX]}"

RUN_ID="q5k${NREFINE}"
DATA_DIR="$OUT_ROOT/${RUN_ID}/shard_$(printf '%02d' "$SHARD_IDX")"
mkdir -p "$DATA_DIR"

EXE="$ROOT/bin/real_space_2pt_test"
test -x "$EXE"

SEED=$((SEED_BASE + GLOBAL_IDX))

THETA_ARGS=()
if [[ -n "${THETA_MIN:-}" ]]; then THETA_ARGS+=(--theta_min "$THETA_MIN"); fi
if [[ -n "${THETA_MAX:-}" ]]; then THETA_ARGS+=(--theta_max "$THETA_MAX"); fi

echo "task=$SGE_TASK_ID n_refine=$NREFINE shard=$SHARD_IDX seed=$SEED data_dir=$DATA_DIR theta=[${THETA_MIN:-0},${THETA_MAX:-pi}]"

"$EXE" --q 5 --n_refine "$NREFINE" --coupling_rule "$COUPLING_RULE" \
  --mesh_mode naive --data_dir "$DATA_DIR" \
  --n_therm "$N_THERM" --n_traj "$N_TRAJ" --n_skip "$N_SKIP" \
  --n_wolff "$N_WOLFF" --n_metropolis "$N_METROPOLIS" \
  --n_ref "$N_REF" --seed "$SEED" "${THETA_ARGS[@]}"
