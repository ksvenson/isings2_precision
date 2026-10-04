#!/usr/bin/env bash
# One array task = one (n_refine, shard_index) bin/save_configs run, for
# the 2026-08-26 offline-M_l / north-pole-two-point-function production
# ladder (see journal.md "Offline M_l[m,m'] construction..." entries and
# "pin down ladder, shard count, and l_max" push). Mirrors
# real_space_shard_task.sh's task-index layout exactly, but drives
# save_configs instead -- no l_max/jack_block_size (save_configs has
# neither; all analysis happens offline afterward via
# scripts/build_ylm_matrix_from_configs.py and
# bin/analyze_north_pole_configs).
#
# naive mesh_mode, exact_sinh coupling_rule -- matches production_2026-08-23
# exactly, so the offline-built M_l/kappa2_r result is directly comparable
# to that run's known l>=3 plateau finding, not confounded by a mesh/
# coupling change.
#
# Required environment (set by submit script via `qsub -v`):
#   ROOT LADDER_CSV N_SHARDS OUT_ROOT COUPLING_RULE
#   N_THERM N_TRAJ N_SKIP N_WOLFF N_METROPOLIS SEED_BASE
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
GLOBAL_IDX=$((SGE_TASK_ID - 1))
SHARD_IDX=$((GLOBAL_IDX % N_SHARDS))
L_IDX=$((GLOBAL_IDX / N_SHARDS))
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$L_IDX]}"

RUN_ID="q5k${NREFINE}"
DATA_DIR="$OUT_ROOT/${RUN_ID}/shard_$(printf '%03d' "$SHARD_IDX")"
mkdir -p "$DATA_DIR"

EXE="$ROOT/bin/save_configs"
test -x "$EXE"

SEED=$((SEED_BASE + GLOBAL_IDX))

echo "task=$SGE_TASK_ID n_refine=$NREFINE shard=$SHARD_IDX seed=$SEED data_dir=$DATA_DIR"

"$EXE" --q 5 --n_refine "$NREFINE" --coupling_rule "$COUPLING_RULE" \
  --mesh_mode naive --data_dir "$DATA_DIR" \
  --n_therm "$N_THERM" --n_traj "$N_TRAJ" --n_skip "$N_SKIP" \
  --n_wolff "$N_WOLFF" --n_metropolis "$N_METROPOLIS" --seed "$SEED"
