#!/usr/bin/env bash
# High-statistics lean_harmonic_stats ladder, per user direction 2026-08-26
# ("run it for 4,6,8,10,12,14,16,20,24,28,32" then "run at high stats with
# the new method") -- one SGE array task per n_refine point (parallel
# across compute nodes, unlike the earlier login-node sequential run this
# supersedes). Detached submission (no -sync y) per workspace convention
# (/projectnb/qfe/misra/CLAUDE.md) -- use qstat/qacct to check on it later.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/lean_harmonic_stats"
test -x "$EXE"

LADDER_CSV="4:6:8:10:12:14:16:20:24:28:32"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_TASKS=${#LADDER[@]}

L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=100000       # 5x the exploratory login-node run's 20000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=500   # n_traj/n_skip/jack_block_size = 100 native blocks
SEED=42                # same seed as the login-node run, for continuity
H_RT=${H_RT:-04:00:00}   # generous ceiling; n_refine=32 took ~94s for 10k
                          # meas on the login node, so 50k meas should be
                          # well under an hour even without SGE speedup

JOB_TAG="${1:-lean_ladder_highstats_$(date +%Y-%m-%d)}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "$JOB_TAG",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "l_max": $L_MAX,
  "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "seed": $SEED,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

echo "submitting $N_TASKS-task lean_harmonic_stats array, project qfe, unpinned queue"
qsub \
  -P qfe \
  -N "s2prec_$JOB_TAG" \
  -j y \
  -o "$LOG_DIR/" \
  -t "1-$N_TASKS" \
  -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED="$SEED" \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
