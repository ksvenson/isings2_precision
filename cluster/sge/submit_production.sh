#!/usr/bin/env bash
# Production spherical-symmetry-test push, per user direction 2026-08-23:
# l_max=8 (vs PLAN.md Phase 1's aliasing-derived value, not yet checked at
# this l_max -- user override, flagged in journal.md), n_refine ladder up to
# L=128, 16 independent-seed shards per (n_refine) point for parallelism.
# One batched SGE array (N_POINTS * N_SHARDS tasks), mirroring
# Twist/cluster/sge/submit_phase2_production.sh's shape. Detached submission
# (no -sync y) per workspace convention (/projectnb/qfe/misra/CLAUDE.md) --
# use qstat/qacct to check on it later.
#
# Deliberately NOT pinned to -q y: as of 2026-08-23 every schedulable host in
# the qfe project's @y buyin group is at or near 28/28 slots (the
# apparently-idle ones are administratively disabled, state 'd', not free).
# -P qfe alone still gives this job qfe's fair-share priority but lets it
# land on any free host cluster-wide, which is how the array actually gets
# parallelism right now -- see journal.md 2026-08-23.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

LADDER_CSV="4:8:16:32:64:128"
N_SHARDS=16
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_TASKS=$((N_POINTS * N_SHARDS))

L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=20000       # per-shard; set from calib_2026-08-23 timing before submitting
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=100   # per feedback_ising_s2_jack_block_size memory: default 1 produces multi-GB jackblocks dumps
SEED_BASE=90000
H_RT=${H_RT:-12:00:00}   # per-task wall clock; generous ceiling, no calib timing available (cluster saturated 2026-08-23)

JOB_TAG="${1:-production_$(date +%Y-%m-%d)}"
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

echo "submitting $N_TASKS-task production array ($N_POINTS ladder points x $N_SHARDS shards), project qfe, unpinned queue"
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
