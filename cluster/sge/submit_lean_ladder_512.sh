#!/usr/bin/env bash
# Naive-vs-equal_area comparison ladder extended to n_refine=512, per
# explicit user direction ("I want a job that goes up to 512" / "make it
# L powers of 2 and powers of 2 times 3 up to 512"). Ladder:
# {4,6,8,12,16,24,32,48,64,96,128,192,256,384,512} (2^k for k=2..9, and
# 3*2^k for k=1..7).
#
# n_traj reduced from the 4..128 ladder's 200000 to 70000 -- per-sweep
# cost scales with n_sites (~(n/128)^2 relative to the empirically
# measured n_refine=128 timing, see CLAUDE.md/journal.md 2026-08-27), so
# keeping n_traj=200000 at n_refine=512 would need ~67h/shard, well past
# any reasonable h_rt. n_traj=70000 keeps the longest point (512) to
# ~24h/shard (h_rt=30:00:00 gives safety margin). This IS an untested
# extrapolation of the linear-in-n_sites scaling past the largest
# n_refine ever run in this campaign (128) -- treat the actual per-point
# wall times as a live calibration check of that assumption, not just
# infrastructure. 32 shards/point (per user direction "I want over a
# million configs for this": n_traj=70000/n_skip=2 = 35000 meas/shard,
# x32 shards = 1,120,000 meas/point) -- 15 points x 32 shards x 2 mesh
# modes = 960 total tasks.
#
# Uses the off-diagonal-capable bin/lean_harmonic_stats build (real
# kappa2_r, not just the diagonal proxy -- see CLAUDE.md 2026-08-27).
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
# SGE's -o log directory must exist before qsub runs, or every task fails
# immediately (Eqw, "can't open output file ... Is a directory") -- hit
# this exact bug earlier in the session on the sibling mesh-compare jobs
# (manual qsub calls that skipped this mkdir), so it's explicit here too.
mkdir -p "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA" "$MESH_CACHE_DIR" \
         "${OUT_ROOT_BASE}_naive_logs" "${OUT_ROOT_BASE}_eqarea_logs"

cat > "$OUT_ROOT_BASE/manifest.json" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "lean_ladder_512_2026-08-27",
  "ladder": [$(IFS=,; echo "${LADDER[*]}")],
  "n_shards": $N_SHARDS, "l_max": $L_MAX, "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "mesh_modes_compared": ["naive", "equal_area"],
  "note": "n_traj reduced from 200000 (the 4..128 ladder's value) to 70000 to keep n_refine=512 within ~24h/shard"
}
EOF

echo "=== step 1: prewarm equal_area mesh cache ($N_POINTS points) ==="
PREWARM_JOB=$(qsub -terse \
  -P qfe -N s2prec_lean_ladder_512_prewarm -j y -o "$OUT_ROOT_BASE/" \
  -sync y -t "1-$N_POINTS" -l h_rt=4:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",L_MAX="$L_MAX",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR",MESH_MODE=equal_area \
  "$ROOT/cluster/sge/prewarm_mesh_cache_task.sh")
echo "prewarm done: $PREWARM_JOB"
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
