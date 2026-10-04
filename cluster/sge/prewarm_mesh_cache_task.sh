#!/usr/bin/env bash
# One array task = one ladder point's equal_area/equal_rp mesh-position
# cache entry. Split out of submit_push7_dual.sh's original *serial,
# login-node* prewarm loop 2026-08-25 after that loop stalled indefinitely
# at n_refine=64 -- root cause was login-node CPU contention (node-wide
# load average ~11-25 against a 1-core-limited interactive session), not a
# bug in EqualizeFaceAreas itself (confirmed by reproducing the stall with
# a GD_DEBUG trace that never got past the initial gradient pass). Each
# ladder point is independent (writes its own cache file), so this runs as
# a normal parallel SGE array instead of a serial login-node loop -- no
# race, since each task only ever touches its own n_refine's cache path.
#
# Generalized 2026-08-25 (same day) to also cover --mesh_mode equal_rp
# (was hardcoded to equal_area) -- suffix/flags must mirror
# production_shard_task.sh's MESH_ARGS logic exactly, same reasoning as
# that script's own comment about run_id mismatches.
#
# Required environment (set by submit script via `qsub -v`):
#   ROOT LADDER_CSV L_MAX EQUAL_AREA_ITERS EQUAL_AREA_STEP MESH_CACHE_DIR
# Optional: MESH_MODE (default equal_area; equal_area|equal_rp)
set -euo pipefail

: "${SGE_TASK_ID:?must run as an SGE array task (qsub -t)}"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
NREFINE="${LADDER[$((SGE_TASK_ID - 1))]}"
MESH_MODE="${MESH_MODE:-equal_area}"

RUN_ID_SUFFIX="_eqarea"
if [[ "$MESH_MODE" == "equal_rp" ]]; then
  RUN_ID_SUFFIX="_eqrp"
fi

EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

SCRATCH_DIR="$MESH_CACHE_DIR/_prewarm_scratch_k${NREFINE}"
mkdir -p "$SCRATCH_DIR/q5k${NREFINE}${RUN_ID_SUFFIX}"

echo "task=$SGE_TASK_ID n_refine=$NREFINE mesh_mode=$MESH_MODE prewarming $MESH_CACHE_DIR/q5k${NREFINE}${RUN_ID_SUFFIX}_step${EQUAL_AREA_STEP}.dat"

"$EXE" --q 5 --n_refine "$NREFINE" --l_max "$L_MAX" --mesh_mode "$MESH_MODE" \
  --equal_area_iters "$EQUAL_AREA_ITERS" --equal_area_step "$EQUAL_AREA_STEP" \
  --mesh_cache_dir "$MESH_CACHE_DIR" --data_dir "$SCRATCH_DIR" \
  --n_therm 1 --n_traj 1 --n_skip 1 --n_wolff 1 --n_metropolis 1 \
  --seed 1 --jack_block_size 1

rm -rf "$SCRATCH_DIR"
echo "task=$SGE_TASK_ID n_refine=$NREFINE done"
