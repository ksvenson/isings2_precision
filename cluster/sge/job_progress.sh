#!/usr/bin/env bash
# Progress of a running (or finished) lean_harmonic_stats push.
#
#   cluster/sge/job_progress.sh campaign_runs/lean_ladder_768_<date> [n_shards]
#
# Two views, because `qstat` alone answers the wrong question. qstat tells
# you what the SCHEDULER is doing; what you usually want to know is how
# much DATA has actually landed. This prints both:
#
#   1. scheduler state, your jobs only, grouped by state
#   2. shard files on disk per (n_refine, mesh_mode), vs expected
#
# Also flags short/truncated shard files -- a walltime kill can cut a
# shard mid-write, which is what corrupted 3 files in the n_refine=512
# push (journal.md 2026-09-02). analyze_lean_harmonic_stats.load() now
# drops such a block safely, but it is still worth seeing.
set -euo pipefail

OUT_ROOT="${1:?usage: job_progress.sh <campaign_runs/...> [n_shards]}"
N_SHARDS="${2:-32}"
OUT_ROOT="${OUT_ROOT%/}"

echo "=== scheduler (your jobs) ==="
if command -v qstat >/dev/null 2>&1; then
  # -g d expands an array job into its individual tasks, otherwise a
  # 352-task array shows up as a single line and you cannot see progress.
  if ! qstat -u "$USER" -g d 2>/dev/null | tail -n +3 | awk '
        {state[$5]++; total++}
        END {
          if (total == 0) { print "  (no jobs queued or running)"; exit }
          for (s in state) printf "  state %-6s %6d task(s)\n", s, state[s]
          printf "  %-12s %6d\n", "TOTAL", total
        }'; then
    echo "  (qstat failed)"
  fi
else
  echo "  (no qstat on this machine -- data view below still works)"
fi
echo "  states: qw=queued  r=running  Eqw=error  t/dr=transferring/deleting"

echo
echo "=== data on disk: $OUT_ROOT ==="
shopt -s nullglob
found_any=0
for mode_dir in "$OUT_ROOT"/*/; do
  mode="$(basename "$mode_dir")"
  files=("$mode_dir"*_lean_*.dat)
  ((${#files[@]})) || continue
  found_any=1
  echo "  $mode/"
  # group by n_refine, parsed out of the driver-generated filename
  printf '%s\n' "${files[@]}" \
    | sed -E 's#.*/q5k([0-9]+)(_eqarea|_eqrp)?_lean_.*#\1#' \
    | sort -n | uniq -c \
    | while read -r count n; do
        bar=$((count * 20 / N_SHARDS)); [[ $bar -gt 20 ]] && bar=20
        printf "    n_refine=%-5s %3d/%-3d shards  [%-20s]%s\n" \
          "$n" "$count" "$N_SHARDS" \
          "$(printf '#%.0s' $(seq 1 $bar) 2>/dev/null || true)" \
          "$([[ $count -ge $N_SHARDS ]] && echo ' complete' || true)"
      done
done
[[ $found_any -eq 1 ]] || echo "  (no *_lean_*.dat files yet -- tasks have not produced output)"

echo
echo "=== truncated-shard check (ragged final line) ==="
bad=0
for f in "$OUT_ROOT"/*/*_lean_*.dat; do
  # every data line must have the same field count; a walltime kill can
  # leave the last one short
  n=$(awk '!/^#/ && NF {print NF}' "$f" | sort -u | wc -l)
  if [[ "$n" -gt 1 ]]; then
    echo "  RAGGED: $f"
    bad=1
  fi
done
[[ $bad -eq 0 ]] && echo "  all shard files have uniform line width"

cat <<EOF

=== other things worth knowing ===
  qstat -u \$USER -g d            # your tasks, array jobs expanded
  qstat -j <jobid>                # why a job is still queued
  qacct -j <jobid> | grep -c 'exit_status  0'   # after it finishes
  tail -f ${OUT_ROOT}_*_logs/*.o<jobid>.<taskid>   # a live task's stdout
EOF
