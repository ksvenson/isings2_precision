#!/usr/bin/env bash
# One array task = one ladder point's M_l/R_l/kappa2_r analysis
# (build_ylm_matrix_from_configs.py over all N_SHARDS of a save_configs
# production run), for the 2026-08-26 offline-M_l analysis. Split off the
# login node: compute_Ml_full's matmul is single-core-bound and takes
# minutes per shard at large n_refine (>3 min/shard at n_refine=128 on a
# 1-core interactive session) -- not something to run interactively.
#
# Required environment (set by submit script via `qsub -v`):
#   ROOT LADDER_CSV OUT_ROOT N_SHARDS L_MAX ANALYSIS_OUT_DIR
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
L_IDX=$((SGE_TASK_ID - 1))
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$L_IDX]}"

RUN_ID="q5k${NREFINE}"
POINT_DIR="$OUT_ROOT/${RUN_ID}"

FILES=()
for s in $(seq -f '%03g' 0 $((N_SHARDS - 1))); do
  SHARD_DIR="$POINT_DIR/shard_${s}"
  POS=$(ls "$SHARD_DIR"/*_positions_*.dat)
  CFG=$(ls "$SHARD_DIR"/*_configs_*.bin)
  FILES+=("$POS" "$CFG")
done

mkdir -p "$ANALYSIS_OUT_DIR"
OUT="$ANALYSIS_OUT_DIR/ml_${RUN_ID}.txt"

echo "task=$SGE_TASK_ID n_refine=$NREFINE n_shards=${#FILES[@]} (files) out=$OUT"

source "$ROOT/.venv_plot/bin/activate"
python3 "$ROOT/scripts/build_ylm_matrix_from_configs.py" "${FILES[@]}" --l_max "$L_MAX" > "$OUT"
echo "wrote $OUT"
