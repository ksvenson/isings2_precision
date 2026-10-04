#!/usr/bin/env bash
# equal_rp (joint circumradius+perimeter+area equalization,
# QfeLatticeS2::EqualizeCircumPerim, arXiv:2407.00459 Appendix B.1/B.2)
# sibling of submit_lean_mesh_compare_100shards.sh, submitted as a THIRD,
# separate job alongside the existing naive/equal_area 100-shard runs
# rather than replacing them -- same ladder/stats/seed-base convention
# (own output dir, own mesh cache, own seed_base range so all three jobs
# stay fully independent and mesh_mode-agnostic downstream scripts don't
# need to pool files across modes). equal_rp has only ever been run
# through ising_s2_crit.cc at low resolution/stats (push8, n_refine<=32,
# inconclusive) -- this is its first run through the lean_harmonic_stats
# driver at the same scale as the existing naive/eqarea kappa2_r plateau
# finding, targeting the still-open Objective (1) spherical-symmetry gate.
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

OUT_ROOT_BASE="$ROOT/campaign_runs/lean_mesh_compare_100shards_2026-08-27"
OUT_ROOT_EQRP="$OUT_ROOT_BASE/eqrp"
MESH_CACHE_DIR="$OUT_ROOT_BASE/mesh_cache"
mkdir -p "$OUT_ROOT_EQRP" "$MESH_CACHE_DIR" "${OUT_ROOT_BASE}_eqrp_logs"

cat > "$OUT_ROOT_BASE/manifest_eqrp.json" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "lean_mesh_compare_100shards_2026-08-27_eqrp",
  "ladder": [$(IFS=,; echo "${LADDER[*]}")],
  "n_shards": $N_SHARDS, "l_max": $L_MAX, "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "mesh_modes_compared": ["naive", "equal_area", "equal_rp"],
  "note": "third mesh mode added alongside existing naive/eqarea 100-shard runs in this stage dir; first equal_rp run through lean_harmonic_stats"
}
EOF

echo "=== step 1: prewarm equal_rp mesh cache ($N_POINTS points) ==="
PREWARM_JOB=$(qsub -terse \
  -P qfe -N s2prec_lean_mesh_compare_100sh_prewarm_eqrp -j y -o "$OUT_ROOT_BASE/" \
  -sync y -t "1-$N_POINTS" -l h_rt=2:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",L_MAX="$L_MAX",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR",MESH_MODE=equal_rp \
  "$ROOT/cluster/sge/prewarm_mesh_cache_task.sh")
echo "prewarm done: $PREWARM_JOB"
STEP_FMT=$(printf '%.3f' "$EQUAL_AREA_STEP")
for n in "${LADDER[@]}"; do
  test -f "$MESH_CACHE_DIR/q5k${n}_eqrp_step${STEP_FMT}.dat" || {
    echo "ERROR: missing cache file for n_refine=$n, aborting" >&2
    exit 1
  }
done
echo "all $N_POINTS cache files verified present"

echo "=== step 2: submit equal_rp production array ==="
N_TASKS=$((N_POINTS * N_SHARDS))

qsub \
  -P qfe -N s2prec_lean_mesh_compare_100sh_eqrp -j y -o "${OUT_ROOT_BASE}_eqrp_logs/" \
  -t "1-$N_TASKS" -l h_rt=12:00:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT_EQRP",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE=820000,MESH_MODE=equal_rp,EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR" \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_lean_mesh_compare_100sh_eqrp"
echo "once complete, compare with:"
echo "  scripts/symmetry_check_lean.py $OUT_ROOT_EQRP"
