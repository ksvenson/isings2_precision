#!/usr/bin/env bash
# push8: first production Ising run using the circumradius+perimeter mesh
# smoother (--mesh_mode equal_rp, implemented 2026-08-25 in
# QfeLatticeS2::EqualizeCircumPerim per arXiv:2407.00459 Appendix B.1/B.2,
# joint objective E=E_R+E_P+E_A per explicit user direction). Only prior
# equal_rp physics run was a login-node diagnostic
# (campaign_runs/diagnostic_eqrp_2026-08-25/, n_refine<=64, 1 seed,
# --coupling_rule duality) -- this is the first production-scale pass, and
# the first to pair equal_rp with --coupling_rule exact_sinh (the
# campaign's primary coupling rule per user direction 2026-08-21).
#
# Ladder: n_refine in {4,6,8,12,16,24,32}, per explicit user direction
# 2026-08-25. l_max=8 (matches production_2026-08-23's baseline ladder and
# lets scripts/cft_symmetry_test.py's Delta_l_pair(a) table -- computed
# for every adjacent l=1..8 pair, not just l=0/1 -- be compared directly
# against that existing dataset and push7's). n_refine=4 (162 sites) is
# below the informal "safe" threshold flagged in submit_push7_dual.sh's
# header (81 modes at l_max=8), same known/accepted risk push7 already ran
# with -- check symmetry_test.py's R_l at l=7,8 for n_refine=4 specifically
# before trusting it.
#
# Stats: N_THERM=2000/N_TRAJ=20000/N_SKIP=2/N_WOLFF=5/N_METROPOLIS=4,
# matching production_2026-08-23's per-shard budget exactly (the l_max=8
# baseline this compares against). N_SHARDS=32 (push7's shard count, not
# production_2026-08-23's 16) for tighter error bars at matched per-shard
# cost. JACK_BLOCK_SIZE=100 -> 200 native blocks/shard, pooled 6400/point,
# gate is 10*(2*8+1)=170 -- ample margin, same value production_2026-08-23
# validated at this n_traj (see directive.md's jack_block_size rule).
#
# EQUAL_AREA_ITERS/STEP reused as-is (flag is shared between equal_area and
# equal_rp, see ising_s2_crit.cc) at the driver's own defaults
# (20000/0.3) -- generous cap relative to the diagnostic run's observed
# convergence (circumradius/perimeter std/mean settled by iter ~5200/8000
# at n_refine=32 in campaign_runs/diagnostic_eqrp_2026-08-25/), so
# early-stop (rel_tol/patience) should trip well before the cap at every
# ladder point here (max n_refine=32).
#
# Mesh cache: shared mesh_cache/ dir, prewarmed serially via SGE array
# (prewarm_mesh_cache_task.sh, generalized this session to take
# MESH_MODE) before the production array is submitted -- same race
# avoidance as push7's equal_area prewarm step.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

LADDER_CSV="4:6:8:12:16:24:32"
N_SHARDS=32
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_TASKS=$((N_POINTS * N_SHARDS))

L_MAX=8
COUPLING_RULE=exact_sinh
MESH_MODE=equal_rp
EQUAL_AREA_ITERS=20000
EQUAL_AREA_STEP=0.3
N_THERM=2000
N_TRAJ=20000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=100
SEED_BASE=800000   # disjoint from production_2026-08-23 (90000+96), push5
                    # (500000+160), push6 (600000+96), push7 naive/eqarea
                    # (700000+416, 710000+416)
H_RT=${H_RT:-2:00:00}

MESH_CACHE_DIR="$ROOT/mesh_cache"
mkdir -p "$MESH_CACHE_DIR"

JOB_TAG="${1:-push8_$(date +%Y-%m-%d)_eqrp_exactsinh}"
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
  "purpose": "first production-scale equal_rp (circumradius+perimeter+area) mesh smoother run, paired with exact_sinh coupling, ladder n_refine={4,6,8,12,16,24,32}",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "mesh_mode": "$MESH_MODE",
  "n_shards": $N_SHARDS,
  "pooled_total_traj_per_point": $((N_SHARDS * N_TRAJ)),
  "l_max": $L_MAX,
  "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "equal_area_iters_cap": $EQUAL_AREA_ITERS,
  "equal_area_step": $EQUAL_AREA_STEP,
  "mesh_cache_dir": "$MESH_CACHE_DIR",
  "seed_base": $SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for k in "${LADDER[@]}"; do
  for s in $(seq -w 0 $((N_SHARDS - 1))); do
    mkdir -p "$OUT_ROOT/q5k${k}_eqrp/shard_${s}/q5k${k}_eqrp"
  done
done

echo "prewarming equal_rp mesh cache via SGE array ($N_POINTS ladder points, waits for completion)"
qsub -sync y \
  -P qfe -N "s2prec_${JOB_TAG}_prewarm" -j y -o "$LOG_DIR/" \
  -t "1-$N_POINTS" -l h_rt=00:30:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",L_MAX="$L_MAX",MESH_MODE="$MESH_MODE",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR" \
  "$ROOT/cluster/sge/prewarm_mesh_cache_task.sh"
STEP_PADDED=$(printf '%.3f' "$EQUAL_AREA_STEP")
for k in "${LADDER[@]}"; do
  test -f "$MESH_CACHE_DIR/q5k${k}_eqrp_step${STEP_PADDED}.dat" || {
    echo "ERROR: prewarm did not produce cache for n_refine=$k -- check $LOG_DIR/s2prec_${JOB_TAG}_prewarm.o*" >&2
    exit 1
  }
done
echo "prewarm complete, all $N_POINTS cache entries confirmed present"

echo "submitting equal_rp production array: $N_TASKS tasks ($N_POINTS ladder points x $N_SHARDS shards), l_max=$L_MAX, mesh_mode=$MESH_MODE, coupling_rule=$COUPLING_RULE, h_rt=$H_RT"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",MESH_MODE="$MESH_MODE",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/production_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
