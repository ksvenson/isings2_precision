#!/usr/bin/env bash
# Submits per-shard M_l block computation (ml_block_shard_task.sh) for
# n_refine=64 and 128 -- the two save_configs_ladder_2026-08-26 points
# whose serial per-point analysis (job 7322113 tasks 5/6) hit the 8h
# h_rt wall-clock limit before finishing (see that task script's header
# comment). 96 shards/point x 2 points = 192 independent tasks, each
# doing exactly one shard's matmul (the actual expensive step) --
# reduce_ml_blocks.py combines the resulting pickles afterward (cheap,
# runs fine on the login node).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LADDER_CSV="64:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_SHARDS=96
N_POINTS=${#LADDER[@]}
N_TASKS=$((N_SHARDS * N_POINTS))

OUT_ROOT="$ROOT/campaign_runs/save_configs_ladder_2026-08-26"
BLOCK_OUT_DIR="$OUT_ROOT/ml_blocks"
L_MAX=8
H_RT=${H_RT:-4:00:00}

JOB_TAG="ml_block_shards_$(date +%Y-%m-%d)"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR" "$BLOCK_OUT_DIR"

echo "submitting M_l per-shard block array: $N_TASKS tasks ($N_POINTS points x $N_SHARDS shards)"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",OUT_ROOT="$OUT_ROOT",N_SHARDS="$N_SHARDS",L_MAX="$L_MAX",BLOCK_OUT_DIR="$BLOCK_OUT_DIR" \
  "$ROOT/cluster/sge/ml_block_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
echo "once complete, reduce with:"
echo "  .venv_plot/bin/python scripts/reduce_ml_blocks.py $BLOCK_OUT_DIR/q5k64/blk_*.pkl > $OUT_ROOT/analysis/ml_q5k64.txt"
echo "  .venv_plot/bin/python scripts/reduce_ml_blocks.py $BLOCK_OUT_DIR/q5k128/blk_*.pkl > $OUT_ROOT/analysis/ml_q5k128.txt"
