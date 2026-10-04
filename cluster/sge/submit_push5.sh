#!/usr/bin/env bash
# push5: l_max 8->12 production push, high stats, 32 shards/point for
# parallelism, per user direction 2026-08-24.
#
# Differences from submit_production.sh (production_2026-08-23, l_max=8):
# - l_max=12. Ladder drops n_refine=4: push4_2026-08-22_lmax12's aliasing
#   diagnostic found l_max=12 is contaminated (kappa2_r blows up) below
#   n_refine=8 -- n_refine=4 would not be trustworthy at this l_max.
# - n_traj=1,000,000/shard (was 20000), n_shards=32 (was 16), per user
#   request for high stats + max parallelism.
# - jack_block_size=4000 (was 100). Native block count = n_traj/block_size;
#   at 100 this would give 10000 blocks/file (~1GB+ per file, the same
#   multi-GB-file problem push4 hit at block_size=1) -- 4000 keeps native
#   blocks at 250, which also exactly matches PLAN.md Phase 1's statistics
#   gate n_blocks >= 10*(2*l_max+1) = 250 at l_max=12.
# - h_rt default bumped to 48:00:00 (was 12:00:00): calibrated on the login
#   node 2026-08-24 at n_refine=128, l_max=12 (2000 traj, 100 therm, 2 skip
#   = 4100 sweeps -> 222s wall = ~0.054 s/sweep, same per-sweep cost as
#   l_max=8 -- Ylm projection is not the bottleneck at this site count).
#   Extrapolated: n_refine=128 shard = (2000 + 1e6*2) sweeps * 0.054s/sweep
#   ~= 30h. 48h leaves margin above that estimate; smaller ladder points
#   finish in minutes and are not the constraint.
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
N_THERM=2000
N_TRAJ=1000000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=4000
SEED_BASE=500000   # disjoint from production_2026-08-23's 90000+96 range
H_RT=${H_RT:-48:00:00}

JOB_TAG="${1:-push5_$(date +%Y-%m-%d)_lmax12}"
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
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "n_shards": $N_SHARDS,
  "l_max": $L_MAX,
  "coupling_rule": "$COUPLING_RULE",
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
    mkdir -p "$OUT_ROOT/q5k${k}/shard_${s}/q5k${k}"
  done
done

echo "submitting $N_TASKS-task push5 array ($N_POINTS ladder points x $N_SHARDS shards), l_max=$L_MAX, project qfe, unpinned queue, h_rt=$H_RT"
qsub \
  -P qfe \
  -N "s2prec_$JOB_TAG" \
  -j y \
  -o "$LOG_DIR/" \
  -t "1-$N_TASKS" \
  -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/production_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
