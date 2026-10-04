#!/usr/bin/env bash
# push7: first production-scale run of --mesh_mode equal_area, to test
# whether the l>=3 kappa2_r plateau found on the naive mesh (PLAN.md step 6
# gate, journal.md 2026-08-24 "ran PLAN.md step-6 gate" entry) is a
# naive-mesh triangle-area-nonuniformity artifact. 2026-08-25 session's
# mesh-geometry diagnostics (test_gd_equal_area.cc, test_gd_symmetry.cc)
# already confirmed EqualizeFaceAreas converges cleanly across the whole
# ladder at n_iter=1000 (area std/mean 0.023/0.012/0.0058/0.0030 at
# n_refine 8/16/32/64) and preserves icosahedral symmetry to machine
# precision -- this is the first time it's run inside the actual Ising
# measurement, not just on the bare mesh.
#
# Same l_max/ladder/coupling_rule/pooled-stats target as push5+push6 (the
# naive-mesh l_max=12 baseline this is meant to compare against directly):
# - Ladder 8:16:32:64:128 (n_refine=4 dropped -- known-aliased at l_max=12,
#   push4_2026-08-22_lmax12).
# - N_TRAJ=31250/shard x N_SHARDS=32 = 1,000,000 pooled traj/point, same as
#   push6's per-point statistics target (push5's n_refine={8,16} points
#   remain the higher-stats naive-mesh baseline for those two resolutions;
#   push6 replaced n_refine={32,64,128} at the 1M-pooled target for
#   wall-clock reasons -- see submit_push6.sh). Using push6's target
#   (not push5's original 32M-pooled) here for a fair equal_area-vs-push6
#   comparison across the whole ladder.
# - JACK_BLOCK_SIZE=1250 -- native blocks/shard = 31250/1250 = 25, pooled
#   across 32 shards = 800 total blocks (same arithmetic as push6),
#   3.2x above PLAN.md's Phase 1 gate n_blocks >= 10*(2*l_max+1) = 250 at
#   l_max=12.
# - h_rt: mesh relaxation itself is negligible next to the MC run (<=30s at
#   n_refine=128 per the 2026-08-25 timing table) -- expect the same
#   ~58min/~15min/~4min (n_refine=128/64/32) and even less at 8/16 as
#   push6's estimate. h_rt kept at push6's 4:00:00 for margin.
#
# New vs. push5/push6: MESH_MODE=equal_area, EQUAL_AREA_ITERS=1000 (the
# driver's own default as of 2026-08-25, passed explicitly here for
# provenance) -- production_shard_task.sh mkdir's/invokes with the
# "_eqarea" run_id suffix the driver's mesh_mode_suffix logic requires.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

LADDER_CSV="8:16:32:64:128"
N_SHARDS=32
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_TASKS=$((N_POINTS * N_SHARDS))

L_MAX=12
COUPLING_RULE=exact_sinh
MESH_MODE=equal_area
EQUAL_AREA_ITERS=1000
N_THERM=2000
N_TRAJ=31250
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=1250
SEED_BASE=700000   # disjoint from production_2026-08-23 (90000+96), push5 (500000+160), push6 (600000+96)
H_RT=${H_RT:-4:00:00}

JOB_TAG="${1:-push7_$(date +%Y-%m-%d)_lmax12_equalarea}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_phase1_production",
  "stage": "$JOB_TAG",
  "compares_against": "push5_2026-08-24_lmax12 + push6_2026-08-24_lmax12_1Mpooled (naive mesh, same ladder/l_max/pooled stats)",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "n_shards": $N_SHARDS,
  "pooled_total_traj_per_point": $((N_SHARDS * N_TRAJ)),
  "l_max": $L_MAX,
  "coupling_rule": "$COUPLING_RULE",
  "mesh_mode": "$MESH_MODE",
  "equal_area_iters": $EQUAL_AREA_ITERS,
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "seed_base": $SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for k in "${LADDER[@]}"; do
  for s in $(seq -w 0 $((N_SHARDS - 1))); do
    mkdir -p "$OUT_ROOT/q5k${k}_eqarea/shard_${s}/q5k${k}_eqarea"
  done
done

echo "submitting $N_TASKS-task push7 array ($N_POINTS ladder points x $N_SHARDS shards), l_max=$L_MAX, mesh_mode=$MESH_MODE, 1M pooled traj/point, project qfe, unpinned queue, h_rt=$H_RT"
qsub \
  -P qfe \
  -N "s2prec_$JOB_TAG" \
  -j y \
  -o "$LOG_DIR/" \
  -t "1-$N_TASKS" \
  -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",MESH_MODE="$MESH_MODE",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/production_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
