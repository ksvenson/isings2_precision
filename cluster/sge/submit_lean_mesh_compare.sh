#!/usr/bin/env bash
# Naive-vs-equal_area comparison ladder for the diagonal-degeneracy
# spherical-symmetry check (scripts/symmetry_check_lean.py), submitted
# 2026-08-27 per explicit user correction ("we should be comparing the
# naive to the equal areas ALWAYS") after the entire overnight
# lean_harmonic_stats ladder turned out to have used the naive mesh only
# -- lean_harmonic_stats.cc had no --mesh_mode flag at all before this
# session (ported from ising_s2_crit.cc, see that file's header comment).
#
# Matches lean_ladder_overnight_2026-08-26's ladder/stats exactly (same
# 19-point ladder, 16 shards/point, n_traj=200000) so the two mesh modes'
# results are directly comparable at matched statistical power, and the
# new equal_area run can be compared point-for-point against the existing
# naive results already analyzed this session.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LADDER_CSV="4:6:8:10:12:14:16:20:24:28:32:40:48:56:64:80:96:112:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=16
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
OUT_ROOT_BASE="$ROOT/campaign_runs/lean_mesh_compare_2026-08-27"
OUT_ROOT_NAIVE="$OUT_ROOT_BASE/naive"
OUT_ROOT_EQAREA="$OUT_ROOT_BASE/eqarea"
MESH_CACHE_DIR="$OUT_ROOT_BASE/mesh_cache"
mkdir -p "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA" "$MESH_CACHE_DIR"

cat > "$OUT_ROOT_BASE/manifest.json" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "lean_mesh_compare_2026-08-27",
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
  -P qfe -N s2prec_lean_mesh_compare_prewarm -j y -o "$OUT_ROOT_BASE/" \
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
  -P qfe -N s2prec_lean_mesh_compare_naive -j y -o "${OUT_ROOT_BASE}_naive_logs/" \
  -t "1-$N_TASKS" -l h_rt=12:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_NAIVE",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=700000,MESH_MODE=naive \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

qsub \
  -P qfe -N s2prec_lean_mesh_compare_eqarea -j y -o "${OUT_ROOT_BASE}_eqarea_logs/" \
  -t "1-$N_TASKS" -l h_rt=12:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_EQAREA",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=710000,MESH_MODE=equal_area,EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR" \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_lean_mesh_compare"
echo "once complete, compare with:"
echo "  scripts/symmetry_check_lean.py $OUT_ROOT_NAIVE"
echo "  scripts/symmetry_check_lean.py $OUT_ROOT_EQAREA"
