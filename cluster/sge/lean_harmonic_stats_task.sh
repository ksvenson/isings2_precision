#!/usr/bin/env bash
# One array task = one (n_refine, shard_index) point of the
# lean_harmonic_stats ladder (bin/lean_harmonic_stats -- online |S_lm|^2 +
# magnetization-moment accumulation, no raw configs / no full M_l matrix,
# see src/lean_harmonic_stats.cc's header). N_SHARDS independent-seed
# shards per ladder point, matching production_shard_task.sh's task-index
# layout/parallelism convention -- more shards means more concurrent SGE
# tasks (spread across more compute nodes) AND more total statistics, at
# the same per-shard n_traj. Shards for the same n_refine are pooled at
# analysis time (scripts/analyze_lean_harmonic_stats.py's combine() -- and
# scripts/plot_lean_ladder.py's ladder driver -- already accept multiple
# files per point, no separate combining step needed).
#
# Task-index layout:
#   global_idx = SGE_TASK_ID - 1
#   shard_idx  = global_idx % N_SHARDS
#   l_idx      = global_idx / N_SHARDS
#   n_refine   = LADDER[l_idx]
#
# Required environment (set by the submit script via `qsub -v`):
#   ROOT LADDER_CSV N_SHARDS OUT_ROOT L_MAX COUPLING_RULE N_THERM N_TRAJ
#   N_SKIP N_WOLFF N_METROPOLIS JACK_BLOCK_SIZE SEED_BASE
# Optional: MESH_MODE (default naive), EQUAL_AREA_ITERS, EQUAL_AREA_STEP,
#   MESH_CACHE_DIR -- added 2026-08-27 once lean_harmonic_stats.cc gained
#   --mesh_mode support (see its header comment); unset MESH_MODE keeps
#   this task script's existing naive-only behavior unchanged.
set -euo pipefail
MESH_MODE="${MESH_MODE:-naive}"

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
GLOBAL_IDX=$((SGE_TASK_ID - 1))
SHARD_IDX=$((GLOBAL_IDX % N_SHARDS))
L_IDX=$((GLOBAL_IDX / N_SHARDS))
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$L_IDX]}"
SEED=$((SEED_BASE + GLOBAL_IDX))

EXE="$ROOT/bin/lean_harmonic_stats"
test -x "$EXE"

mkdir -p "$OUT_ROOT"

echo "task=$SGE_TASK_ID n_refine=$NREFINE shard=$SHARD_IDX seed=$SEED out_root=$OUT_ROOT"

EXTRA_ARGS=(--mesh_mode "$MESH_MODE")
if [[ "$MESH_MODE" != "naive" ]]; then
  EXTRA_ARGS+=(--equal_area_iters "${EQUAL_AREA_ITERS:-20000}" --equal_area_step "${EQUAL_AREA_STEP:-0.3}")
  if [[ -n "${MESH_CACHE_DIR:-}" ]]; then
    EXTRA_ARGS+=(--mesh_cache_dir "$MESH_CACHE_DIR")
  fi
fi

"$EXE" --q 5 --n_refine "$NREFINE" --l_max "$L_MAX" \
  --coupling_rule "$COUPLING_RULE" --data_dir "$OUT_ROOT" \
  --n_therm "$N_THERM" --n_traj "$N_TRAJ" --n_skip "$N_SKIP" \
  --n_wolff "$N_WOLFF" --n_metropolis "$N_METROPOLIS" --seed "$SEED" \
  --jack_block_size "$JACK_BLOCK_SIZE" "${EXTRA_ARGS[@]}"
