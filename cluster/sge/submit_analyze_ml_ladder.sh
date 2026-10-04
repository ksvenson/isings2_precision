#!/usr/bin/env bash
# Submits the M_l/R_l/kappa2_r analysis (build_ylm_matrix_from_configs.py)
# for each point of save_configs_ladder_2026-08-26 as its own SGE array
# task -- see analyze_ml_ladder_task.sh header comment for why this can't
# run on the login node at large n_refine.
#
# One task per ladder point (not per shard): each task loads all 96
# shards for its n_refine and does the per-shard M_l matmul + delete-one-
# shard jackknife itself (build_ylm_matrix_from_configs.py's existing
# per-shard loop) -- no need for a separate reduction step. h_rt is set
# generously per the observed >3min/shard at n_refine=128 on a 1-core
# session (worst case ~96*3min=4.8h; cluster compute nodes may do better,
# but budget for the pessimistic case).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LADDER_CSV="4:8:16:32:64:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}

OUT_ROOT="$ROOT/campaign_runs/save_configs_ladder_2026-08-26"
ANALYSIS_OUT_DIR="$OUT_ROOT/analysis"
N_SHARDS=96
L_MAX=8
H_RT=${H_RT:-8:00:00}

JOB_TAG="analyze_ml_ladder_$(date +%Y-%m-%d)"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR" "$ANALYSIS_OUT_DIR"

echo "submitting M_l analysis array: $N_POINTS tasks (one per ladder point)"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -t "1-$N_POINTS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",OUT_ROOT="$OUT_ROOT",N_SHARDS="$N_SHARDS",L_MAX="$L_MAX",ANALYSIS_OUT_DIR="$ANALYSIS_OUT_DIR" \
  "$ROOT/cluster/sge/analyze_ml_ladder_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
