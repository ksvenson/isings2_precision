#!/usr/bin/env bash
# One array task = one n_refine point of the EqualizeFaceAreas performance
# sweep (bin/test_gd_stress, ad hoc diagnostic binary, not the production
# driver). Points are fully independent (no mesh cache, no shared state),
# so this runs as a normal parallel SGE array -- unlike the mesh-cache
# prewarm, there is nothing here that needs `-sync y` or serialization.
#
# Required environment (set by submit_gd_stress_1024.sh via `qsub -v`):
#   ROOT LADDER_CSV N_ITER STEP OUT_DIR
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$((SGE_TASK_ID - 1))]}"

EXE="$ROOT/bin/test_gd_stress"
test -x "$EXE"

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/gd_stress_k${NREFINE}.log"

echo "task=$SGE_TASK_ID n_refine=$NREFINE n_iter=$N_ITER step=$STEP" | tee "$LOG"
/usr/bin/time -v "$EXE" "$NREFINE" "$N_ITER" "$STEP" >>"$LOG" 2>&1
echo "task=$SGE_TASK_ID n_refine=$NREFINE done" | tee -a "$LOG"
