#!/usr/bin/env bash
# One array task = one (n_refine, shard_index) production measurement for
# the spherical-symmetry-test ladder (PLAN.md Phase 1), 16 shards/point per
# user direction 2026-08-23. Mirrors Twist/cluster/sge/phase2_production_shard_task.sh's
# task-index layout.
#
# Task-index layout:
#   global_idx = SGE_TASK_ID - 1
#   shard_idx  = global_idx % N_SHARDS
#   l_idx      = global_idx / N_SHARDS
#   n_refine   = LADDER[l_idx]
#
# N_TRAJ/N_THERM/N_SKIP are per-shard (not split, unlike Twist's per-point
# total-split convention) since each shard here is an independent seeded
# run, not a fraction of one point's trajectory budget.
#
# Required environment (set by submit_production.sh via `qsub -v`):
#   ROOT LADDER_CSV N_SHARDS OUT_ROOT L_MAX COUPLING_RULE
#   N_THERM N_TRAJ N_SKIP N_WOLFF N_METROPOLIS SEED_BASE
# Optional environment (default to naive-mesh behavior if unset, so
# existing pushes' -v lists that predate --mesh_mode keep working
# unchanged):
#   MESH_MODE (default naive; naive|equal_area|equal_rp) EQUAL_AREA_ITERS
#   (default 20000, used for both equal_area and equal_rp -- the driver
#   reuses the same --equal_area_iters/--equal_area_step flags for
#   equal_rp's relaxation loop, see ising_s2_crit.cc) MESH_CACHE_DIR
#   (default unset -- caching disabled; only used when MESH_MODE=equal_area
#   or equal_rp. Submit scripts that parallelize many shards per
#   (n_refine, mesh_mode) point should pre-populate this cache serially
#   before submitting the array -- concurrent shards hitting a cold cache
#   would race to fopen(...,"wb") the same file. See
#   submit_push7_dual.sh's prewarm step.)
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
GLOBAL_IDX=$((SGE_TASK_ID - 1))
SHARD_IDX=$((GLOBAL_IDX % N_SHARDS))
L_IDX=$((GLOBAL_IDX / N_SHARDS))
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$L_IDX]}"

MESH_MODE="${MESH_MODE:-naive}"
# Must mirror the driver's own mesh_mode_suffix logic (ising_s2_crit.cc)
# exactly -- it builds run_id from this suffix and does not mkdir -p its
# output dir, so a mismatch here means the run silently writes nowhere the
# driver can find/reopen.
RUN_ID_SUFFIX=""
if [[ "$MESH_MODE" == "equal_area" ]]; then
  RUN_ID_SUFFIX="_eqarea"
elif [[ "$MESH_MODE" == "equal_rp" ]]; then
  RUN_ID_SUFFIX="_eqrp"
fi
RUN_ID="q5k${NREFINE}${RUN_ID_SUFFIX}"

# Driver appends /<run_id>/ itself and does not mkdir -p it, so create
# that nested dir; --data_dir is the shard's base.
DATA_DIR="$OUT_ROOT/${RUN_ID}/shard_$(printf '%02d' "$SHARD_IDX")"
mkdir -p "$DATA_DIR/${RUN_ID}"

EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

SEED=$((SEED_BASE + GLOBAL_IDX))

echo "task=$SGE_TASK_ID n_refine=$NREFINE shard=$SHARD_IDX seed=$SEED data_dir=$DATA_DIR mesh_mode=$MESH_MODE"

MESH_ARGS=(--mesh_mode "$MESH_MODE")
if [[ "$MESH_MODE" == "equal_area" || "$MESH_MODE" == "equal_rp" ]]; then
  MESH_ARGS+=(--equal_area_iters "${EQUAL_AREA_ITERS:-20000}")
  if [[ -n "${EQUAL_AREA_STEP:-}" ]]; then
    MESH_ARGS+=(--equal_area_step "$EQUAL_AREA_STEP")
  fi
  if [[ -n "${MESH_CACHE_DIR:-}" ]]; then
    MESH_ARGS+=(--mesh_cache_dir "$MESH_CACHE_DIR")
  fi
fi

"$EXE" --q 5 --n_refine "$NREFINE" --l_max "$L_MAX" \
  --coupling_rule "$COUPLING_RULE" --data_dir "$DATA_DIR" \
  --n_therm "$N_THERM" --n_traj "$N_TRAJ" --n_skip "$N_SKIP" \
  --n_wolff "$N_WOLFF" --n_metropolis "$N_METROPOLIS" --seed "$SEED" \
  --jack_block_size "${JACK_BLOCK_SIZE:-100}" "${MESH_ARGS[@]}"
