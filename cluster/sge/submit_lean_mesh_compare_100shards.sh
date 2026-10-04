#!/usr/bin/env bash
# 100-shards/point sibling of submit_lean_mesh_compare.sh, submitted as a
# SECOND, separate job per explicit user direction, alongside (not
# replacing) the 16-shard version -- same ladder/stats otherwise (own
# output dir, own mesh cache, own seed_base range to keep both jobs
# fully independent). 100 shards/point = 6.25x the total statistics
# (10M vs 1.6M measurements/point) at the SAME per-task wall time (each
# shard still runs the same n_traj=200000 regardless of shard count,
# since shards run in parallel, not serially) -- see CLAUDE.md/journal.md
# 2026-08-27 for the wall-clock/compute-cost reasoning that motivated
# this. Uses the same off-diagonal-capable bin/lean_harmonic_stats build
# as the 16-shard job (rebuilt earlier this session, validated against
# an independent code path) -- scripts/symmetry_check_lean_kappa2.py
# gives the real R_l/kappa2_r from this data once it completes.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LADDER_CSV="4:6:8:10:12:14:16:20:24:28:32:40:48:56:64:80:96:112:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=100
L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=200000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=1000
EQUAL_AREA_ITERS=20000
EQUAL_AREA_STEP=0.3

# Separate output directories per mesh mode -- lean_harmonic_stats.cc
# writes flat into <data_dir>/<run_id>_lean_<seed>.dat with only a
# filename suffix (run_id="q5k<n>_eqarea" vs "q5k<n>") distinguishing
# mesh modes, and the existing analysis scripts' globs
# (symmetry_check_lean.py, plot_lean_ladder.py: "q5k*_lean_*.dat") are NOT
# anchored enough to reliably exclude the other mode's files if both sat
# in the same directory (glob "*" freely skips over "_eqarea") -- would
# either crash (regex n_refine extraction fails on the eqarea filename) or
# silently pool both mesh modes together. Keeping them in separate
# directories sidesteps this rather than hardening every downstream glob.
OUT_ROOT_BASE="$ROOT/campaign_runs/lean_mesh_compare_100shards_2026-08-27"
OUT_ROOT_NAIVE="$OUT_ROOT_BASE/naive"
OUT_ROOT_EQAREA="$OUT_ROOT_BASE/eqarea"
MESH_CACHE_DIR="$OUT_ROOT_BASE/mesh_cache"
mkdir -p "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA" "$MESH_CACHE_DIR"

cat > "$OUT_ROOT_BASE/manifest.json" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "lean_mesh_compare_100shards_2026-08-27",
  "ladder": [$(IFS=,; echo "${LADDER[*]}")],
  "n_shards": $N_SHARDS, "l_max": $L_MAX, "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "mesh_modes_compared": ["naive", "equal_area"],
  "note": "naive results here should match lean_ladder_overnight_2026-08-26 (same ladder/stats, different seed_base)"
}
EOF

echo "=== step 1: prewarm equal_area mesh cache ($N_POINTS points) ==="
PREWARM_JOB=$(qsub -terse \
  -P qfe -N s2prec_lean_mesh_compare_100sh_prewarm -j y -o "$OUT_ROOT_BASE/" \
  -sync y -t "1-$N_POINTS" -l h_rt=2:00:00 \
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

echo "=== step 2: submit naive + equal_area production arrays (separate output dirs) ==="
N_TASKS=$((N_POINTS * N_SHARDS))

qsub \
  -P qfe -N s2prec_lean_mesh_compare_100sh_naive -j y -o "${OUT_ROOT_BASE}_naive_logs/" \
  -t "1-$N_TASKS" -l h_rt=12:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_NAIVE",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=800000,MESH_MODE=naive \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

qsub \
  -P qfe -N s2prec_lean_mesh_compare_100sh_eqarea -j y -o "${OUT_ROOT_BASE}_eqarea_logs/" \
  -t "1-$N_TASKS" -l h_rt=12:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_EQAREA",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=810000,MESH_MODE=equal_area,EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR" \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_lean_mesh_compare_100sh"
echo "once complete, compare with:"
echo "  scripts/symmetry_check_lean.py $OUT_ROOT_NAIVE"
echo "  scripts/symmetry_check_lean.py $OUT_ROOT_EQAREA"
