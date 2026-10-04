#!/usr/bin/env bash
# EqualizeFaceAreas performance sweep, extended to n_refine=1024, per user
# direction 2026-08-25 ("schedule the runs necessary to see the
# performance of the optimizer for L powers of 2 and powers of 2 times 3
# up to 1024. give this lots of wall time and lots of power").
#
# Ladder: 2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024
# (19 points) -- same powers-of-2 / powers-of-2x3 pattern as
# campaign_runs/gd_stress_2026-08-25/sweep_pow2_pow2x3_upto128.log,
# extended past 128. Supersedes sweep_large_256_to_2048.log, which used
# the pre-fix binary (bad early-stop bug, see journal.md 2026-08-25
# "iters-needed model" entry) and was never rerun.
#
# Extrapolating from the 2026-08-25 fixed-binary benchmark (per-iteration
# cost approx linear in n_faces=20*n_refine^2, and every point from
# n_refine=128 up already hits the 20000-iteration cap without
# early-stopping): n_refine=1024 is expected to take on the order of
# ~10-11 hours single-threaded, hence the generous h_rt/mem below rather
# than a calibrated exact figure -- this is an exploratory stress test,
# not a re-run of an already-timed point.
#
# One SGE array task per ladder point (independent, no shared state --
# unlike the mesh-cache prewarm, no serialization/`-sync y` needed here).
# Each task requests a multi-core parallel environment purely for the
# memory reservation (test_gd_stress itself is single-threaded) since SCC
# allocates memory per-core -- sized generously for the largest point's
# ~10.5M-site / ~21M-face mesh.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/test_gd_stress"
test -x "$EXE"

LADDER_CSV="2:3:4:6:8:12:16:24:32:48:64:96:128:192:256:384:512:768:1024"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}

N_ITER=20000
STEP=0.3
H_RT=${H_RT:-72:00:00}
N_CORES=${N_CORES:-16}   # memory reservation only; binary is single-threaded

JOB_TAG="gd_stress_1024_$(date +%Y-%m-%d)"
OUT_DIR="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_DIR"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

cat >"$OUT_DIR/manifest.json" <<EOF
{
  "purpose": "EqualizeFaceAreas performance sweep, powers of 2 and 2x3 up to n_refine=1024",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "n_iter_cap": $N_ITER,
  "step": $STEP,
  "h_rt": "$H_RT",
  "n_cores_reserved": $N_CORES,
  "binary": "bin/test_gd_stress",
  "simulator_checksum": "$(sha256sum "$EXE" | awk '{print $1}')"
}
EOF
echo "manifest written to $OUT_DIR/manifest.json"

echo "submitting gd_stress array: $N_POINTS tasks, h_rt=$H_RT, $N_CORES cores/task reserved for memory"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -pe omp "$N_CORES" \
  -t "1-$N_POINTS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_ITER="$N_ITER",STEP="$STEP",OUT_DIR="$OUT_DIR" \
  "$ROOT/cluster/sge/gd_stress_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
echo "results land in $OUT_DIR/gd_stress_k<n_refine>.log"
