#!/usr/bin/env bash
# One array task = one (n_refine, shard_index) M_l block computation
# (scripts/compute_ml_block_shard.py), for the two ladder points
# (n_refine=64, 128) whose serial per-point analysis
# (analyze_ml_ladder_task.sh, job 7322113 tasks 5/6) hit SGE's h_rt
# wall-clock limit before finishing all 96 shards -- see
# CLAUDE.md/journal.md 2026-08-27 "kappa2_r method, not the diagonal
# proxy" entry. Splits the same per-shard work
# (build_ylm_matrix_from_configs.shard_to_block) into independently
# schedulable tasks instead of one 96-shard serial loop, mirroring every
# other production push's shard-then-reduce pattern in this campaign.
#
# Required environment (set by submit script via `qsub -v`):
#   ROOT LADDER_CSV OUT_ROOT N_SHARDS L_MAX BLOCK_OUT_DIR
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
GLOBAL_IDX=$((SGE_TASK_ID - 1))
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}

L_IDX=$((GLOBAL_IDX / N_SHARDS))
SHARD_IDX=$((GLOBAL_IDX % N_SHARDS))
NREFINE="${LADDER[$L_IDX]}"

RUN_ID="q5k${NREFINE}"
SHARD_STR=$(printf '%03d' "$SHARD_IDX")
SHARD_DIR="$OUT_ROOT/${RUN_ID}/shard_${SHARD_STR}"
POS=$(ls "$SHARD_DIR"/*_positions_*.dat)
CFG=$(ls "$SHARD_DIR"/*_configs_*.bin)

mkdir -p "$BLOCK_OUT_DIR/${RUN_ID}"
OUT="$BLOCK_OUT_DIR/${RUN_ID}/blk_${SHARD_STR}.pkl"

echo "task=$SGE_TASK_ID n_refine=$NREFINE shard=$SHARD_STR pos=$POS cfg=$CFG out=$OUT"

source "$ROOT/.venv_plot/bin/activate"
python3 "$ROOT/scripts/compute_ml_block_shard.py" "$POS" "$CFG" --l_max "$L_MAX" -o "$OUT"
