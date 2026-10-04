#!/usr/bin/env bash
# Recovery for submit_lean_ladder_512.sh: task 15 (n_refine=512) of the
# original 15-point prewarm array (job 7342648) was killed by its
# h_rt=4:00:00 limit (qacct: exit_status=137, ru_wallclock=14401s, right
# at the wall-clock cap, maxvmem only 2.1GB -- a walltime kill, not OOM).
# The other 14 points' cache files are already on disk and untouched.
# This script only redoes the missing n_refine=512 point (much larger
# h_rt=24:00:00, generous margin since the failure mode was purely
# walltime), then re-verifies all 15 cache files and runs the exact same
# step-2 production submission as the original script.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LADDER_CSV="4:6:8:12:16:24:32:48:64:96:128:192:256:384:512"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=32
L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=70000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=1000
EQUAL_AREA_ITERS=20000
EQUAL_AREA_STEP=0.3

OUT_ROOT_BASE="$ROOT/campaign_runs/lean_ladder_512_2026-08-27"
OUT_ROOT_NAIVE="$OUT_ROOT_BASE/naive"
OUT_ROOT_EQAREA="$OUT_ROOT_BASE/eqarea"
MESH_CACHE_DIR="$OUT_ROOT_BASE/mesh_cache"
mkdir -p "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA" "$MESH_CACHE_DIR" \
         "${OUT_ROOT_BASE}_naive_logs" "${OUT_ROOT_BASE}_eqarea_logs"

echo "=== step 1: prewarm ONLY the missing n_refine=512 point (task 15/15) ==="
PREWARM_JOB=$(qsub -terse \
  -P qfe -N s2prec_lean_ladder_512_prewarm_n512 -j y -o "$OUT_ROOT_BASE/" \
  -sync y -t "15-15" -l h_rt=24:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",L_MAX="$L_MAX",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR",MESH_MODE=equal_area \
  "$ROOT/cluster/sge/prewarm_mesh_cache_task.sh")
echo "prewarm(n=512) done: $PREWARM_JOB"

STEP_FMT=$(printf '%.3f' "$EQUAL_AREA_STEP")
for n in "${LADDER[@]}"; do
  test -f "$MESH_CACHE_DIR/q5k${n}_eqarea_step${STEP_FMT}.dat" || {
    echo "ERROR: missing cache file for n_refine=$n, aborting" >&2
    exit 1
  }
done
echo "all $N_POINTS cache files verified present"

echo "=== step 2: submit naive + equal_area production arrays ==="
N_TASKS=$((N_POINTS * N_SHARDS))

qsub \
  -P qfe -N s2prec_lean_ladder_512_naive -j y -o "${OUT_ROOT_BASE}_naive_logs/" \
  -t "1-$N_TASKS" -l h_rt=30:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_NAIVE",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=900000,MESH_MODE=naive \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

qsub \
  -P qfe -N s2prec_lean_ladder_512_eqarea -j y -o "${OUT_ROOT_BASE}_eqarea_logs/" \
  -t "1-$N_TASKS" -l h_rt=30:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_EQAREA",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=910000,MESH_MODE=equal_area,EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR" \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_lean_ladder_512"
echo "once complete, compare with:"
echo "  scripts/symmetry_check_lean_kappa2.py $OUT_ROOT_NAIVE"
echo "  scripts/symmetry_check_lean_kappa2.py $OUT_ROOT_EQAREA"
