#!/usr/bin/env bash
# **OBSOLETE as of 2026-08-26** -- real_space_2pt_test.cc's rewrite dropped
# the z=cos(theta) binning path this script existed to compare against, so
# there is no longer a separate "orbit_avg" toggle: exact icosahedral-orbit
# averaging is now the only mode the driver has (see that file's header
# comment). Kept for historical reference (produced job 7311202, part of
# the pre-rewrite CLEAN CONFIRMATION result in journal.md) but do not reuse
# as a template -- submit_real_space_ladder.sh alone is now the current
# production-ladder submit script.
#
# Orbit-averaged twin of submit_real_space_ladder.sh (job 7311120), built
# 2026-08-25 per user direction ("build it in next") after
# test_orbit_avg.cc validated icosahedral-orbit averaging as an unbiased
# ~4-5x free variance reduction and real_space_2pt_test.cc's (pre-rewrite)
# --orbit_avg flag confirmed it agrees with the z-binned method everywhere
# except the single shortest-distance bin (already excluded by the
# existing theta_min=0.3 fit window) -- see journal.md's "built
# --orbit_avg into real_space_2pt_test.cc" entry for the full comparison.
#
# Same ladder/coupling/mesh/stats as submit_real_space_ladder.sh (matched
# design so the two pushes' Delta_inf extrapolations are directly
# comparable) with ORBIT_AVG=1 added and a disjoint seed range.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/real_space_2pt_test"
test -x "$EXE"

LADDER_CSV="8:16:24:32:48:64:96:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=24
N_TASKS=$((N_POINTS * N_SHARDS))

COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=20000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
N_BINS=100
N_REF=300
ORBIT_AVG=1
SEED_BASE=810000  # disjoint from submit_real_space_ladder.sh's 800000
H_RT=${H_RT:-2:00:00}

JOB_TAG="${1:-real_space_ladder_orbit_$(date +%Y-%m-%d)}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_real_space_delta_extraction",
  "stage": "$JOB_TAG",
  "purpose": "orbit-averaged twin of real_space_ladder_2026-08-25 (job 7311120) -- icosahedral-orbit-averaged real-space Delta(a) ladder, for direct comparison against the z-binned method",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "mesh_mode": "naive",
  "coupling_rule": "$COUPLING_RULE",
  "orbit_avg": true,
  "n_shards_per_point": $N_SHARDS,
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "n_bins": $N_BINS, "n_ref": $N_REF,
  "seed_base": $SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for k in "${LADDER[@]}"; do
  for s in $(seq -w 0 $((N_SHARDS - 1))); do
    mkdir -p "$OUT_ROOT/q5k${k}/shard_${s}"
  done
done

echo "submitting orbit-averaged real-space ladder array: $N_TASKS tasks ($N_POINTS ladder points x $N_SHARDS shards)"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",N_BINS="$N_BINS",N_REF="$N_REF",SEED_BASE="$SEED_BASE",ORBIT_AVG="$ORBIT_AVG" \
  "$ROOT/cluster/sge/real_space_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
