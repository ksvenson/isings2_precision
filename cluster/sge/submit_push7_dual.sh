#!/usr/bin/env bash
# push7: first head-to-head naive-vs-equal_area production comparison, per
# user direction 2026-08-25. Supersedes the never-submitted
# submit_push7_equalarea.sh (a solo equal_area l_max=12 ladder extending
# push5/push6 -- kept on disk as an unused reference, not deleted, but not
# what "push7" now refers to).
#
# Design, per explicit user direction this session:
# - Ladder: n_refine in {2,3,4,6,8,12,16,24,32,48,64,96} -- powers of 2 and
#   powers of 2 times 3, up to 96 (12 points). NOT the l_max=12 aliasing-
#   safe {8..128} ladder used by push5/push6 -- this ladder deliberately
#   goes down to n_refine=2 to see the mesh-mode comparison across the
#   full refinement range, including points too coarse for l_max=12.
# - Two mesh_mode arms, both run: naive and equal_area -- the actual
#   comparison this push exists to make (does equal_area's better area
#   uniformity, verified geometrically in campaign_runs/gd_stress_2026-08-25/,
#   shrink the l>=3 kappa2_r plateau relative to naive at matched stats?).
#   A third "unprojected icosahedron" arm was discussed and then explicitly
#   dropped by the user this session -- not implemented, do not add without
#   asking again (see journal.md 2026-08-25 "push7 planning" entry for the
#   dropped-scope discussion).
# - l_max=8, per explicit user direction (overriding this script's
#   original l_max=3 choice). Not derived from a per-resolution aliasing
#   check for this ladder (PLAN.md step 3 still open) -- production_2026-
#   08-23 previously ran l_max=8 as an unconfirmed override down to
#   n_refine=4 (81 modes vs 162 sites there), itself flagged as only a
#   "plausible" risk, never separately verified. This push goes lower
#   still, to n_refine=2 (42 sites) and n_refine=3 (92 sites) -- both well
#   below the 81-mode count l_max=8 needs, so those two points in
#   particular are a real, currently-unquantified aliasing risk. Flagging
#   explicitly per user direction to proceed anyway; check
#   scripts/symmetry_test.py's R_l/kappa2_r at l=7,8 for n_refine=2,3
#   specifically before drawing any conclusion from them (an aliasing
#   blowup there would look like, and be mistaken for, an SO(3)-breaking
#   signal if not checked separately from the larger-n_refine points).
# - 10,000,000 total pooled trajectories, split evenly across all 24
#   (n_refine, mesh_mode) combinations (416,000 each) -- not per-point like
#   push5/push6's convention, per explicit user direction ("10 million
#   configs total", revised up from an initial "one million total").
# - N_SHARDS=32 per combo (768 tasks total across both arrays) for
#   "aggressive parallelism" per user direction -- same shard count
#   push5/push6 already confirmed schedules immediately with no queueing
#   on this project's allocation.
# - Mesh-position cache (--mesh_cache_dir, added this session to
#   ising_s2_crit.cc/S2.h) is pre-populated for every ladder point before
#   the equal_area array is submitted -- 32 shards per point would
#   otherwise race to fopen(...,"wb") the same cache file concurrently on
#   a cold cache. **Revised 2026-08-25 (later same day)**: the original
#   design ran this prewarm as a *serial loop on the login node* -- it
#   stalled indefinitely at n_refine=64, well past the ~11s the
#   isolated gd_stress_2026-08-25 benchmark predicted. Root cause: the
#   login node this session ran on was CPU-contended (node-wide load
#   average ~11-25 across 35 users, this session limited to 1 core), not
#   an EqualizeFaceAreas bug -- confirmed by reproducing the stall live
#   (GD_DEBUG trace never advanced past the mesh_cache-miss line) with
#   `nproc`=1 and `uptime` showing the load. Fixed by submitting the
#   prewarm itself as a small `qsub -sync y` SGE array (one task per
#   ladder point, independent cache files so no race) instead of running
#   it inline on the login node -- see prewarm_mesh_cache_task.sh.
# - Ladder extended 2026-08-25 (later same day, user direction) to add
#   n_refine=128 (13 points total, was 12/capped at 96) to match L_MAX=8's
#   already-confirmed-safe aliasing range (push4_2026-08-22_lmax12 checked
#   n_refine>=8 safe at l_max=12, so l_max=8 is safe there too) and to
#   extend the mesh-mode comparison to the same top resolution
#   production_2026-08-23/push5/push6 used.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

LADDER_CSV="2:3:4:6:8:12:16:24:32:48:64:96:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=32
N_TASKS=$((N_POINTS * N_SHARDS))

L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=13000        # 32 shards x 13000 = 416,000 pooled/combo x 24 combos = 9,984,000 (~10M)
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=100  # 130 native blocks/shard, pooled 4160/combo -- gate is 10*(2*8+1)=170
EQUAL_AREA_ITERS=20000
EQUAL_AREA_STEP=0.3
NAIVE_SEED_BASE=700000    # disjoint from production_2026-08-23 (90000+96),
EQAREA_SEED_BASE=710000   # push5 (500000+160), push6 (600000+96)
H_RT=${H_RT:-2:00:00}

MESH_CACHE_DIR="$ROOT/mesh_cache"
mkdir -p "$MESH_CACHE_DIR"

JOB_TAG="${1:-push7_$(date +%Y-%m-%d)_dual_naive_eqarea}"
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
  "purpose": "naive vs equal_area mesh_mode head-to-head at matched stats, ladder n_refine={2,3,4,6,8,12,16,24,32,48,64,96,128}",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "mesh_modes": ["naive", "equal_area"],
  "n_shards_per_combo": $N_SHARDS,
  "pooled_total_traj_per_combo": $((N_SHARDS * N_TRAJ)),
  "pooled_total_traj_all_combos": $((N_SHARDS * N_TRAJ * N_POINTS * 2)),
  "l_max": $L_MAX,
  "l_max_note": "conservative fixed choice (not a per-resolution aliasing check) -- safe for n_refine=2 (42 sites), see script header",
  "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "equal_area_iters_cap": $EQUAL_AREA_ITERS,
  "equal_area_step": $EQUAL_AREA_STEP,
  "mesh_cache_dir": "$MESH_CACHE_DIR",
  "naive_seed_base": $NAIVE_SEED_BASE,
  "eqarea_seed_base": $EQAREA_SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for mode_suffix in "" "_eqarea"; do
  for k in "${LADDER[@]}"; do
    for s in $(seq -w 0 $((N_SHARDS - 1))); do
      mkdir -p "$OUT_ROOT/q5k${k}${mode_suffix}/shard_${s}/q5k${k}${mode_suffix}"
    done
  done
done

echo "prewarming equal_area mesh cache via SGE array ($N_POINTS ladder points, waits for completion)"
qsub -sync y \
  -P qfe -N "s2prec_${JOB_TAG}_prewarm" -j y -o "$LOG_DIR/" \
  -t "1-$N_POINTS" -l h_rt=00:30:00 \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",L_MAX="$L_MAX",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR" \
  "$ROOT/cluster/sge/prewarm_mesh_cache_task.sh"
STEP_PADDED=$(printf '%.3f' "$EQUAL_AREA_STEP")
for k in "${LADDER[@]}"; do
  test -f "$MESH_CACHE_DIR/q5k${k}_step${STEP_PADDED}.dat" || {
    echo "ERROR: prewarm did not produce cache for n_refine=$k -- check $LOG_DIR/s2prec_${JOB_TAG}_prewarm.o*" >&2
    exit 1
  }
done
echo "prewarm complete, all $N_POINTS cache entries confirmed present"

echo "submitting naive array: $N_TASKS tasks ($N_POINTS ladder points x $N_SHARDS shards)"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}_naive" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",MESH_MODE="naive",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$NAIVE_SEED_BASE" \
  "$ROOT/cluster/sge/production_shard_task.sh"

echo "submitting equal_area array: $N_TASKS tasks ($N_POINTS ladder points x $N_SHARDS shards), mesh cache prewarmed"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}_eqarea" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",MESH_MODE="equal_area",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",MESH_CACHE_DIR="$MESH_CACHE_DIR",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$EQAREA_SEED_BASE" \
  "$ROOT/cluster/sge/production_shard_task.sh"

echo "submitted (detached, both arrays) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
