#!/bin/bash
# MATLAB-vs-Python sweep of load_limb_sounders (see README.md).
# Usage: run_sweep.sh WORK_FOLDER
# Environment: ATMOSGW_MATLAB_CODE (required), MATLAB, PYTHON, SWEEP_JOBS,
# SWEEP_ONLY (optional grep pattern of batch files, e.g. 'batch_05[0-2]').
set -u
HERE=$(dirname "$(readlink -f "$0")")
WORK=$(readlink -f "$1")
: "${ATMOSGW_MATLAB_CODE:?set ATMOSGW_MATLAB_CODE to the folder of the MATLAB code}"
MATLAB=${MATLAB:-matlab}
PYTHON=${PYTHON:-python}
SWEEP_JOBS=${SWEEP_JOBS:-4}
export HERE WORK MATLAB PYTHON ATMOSGW_MATLAB_CODE
mkdir -p "$WORK/logs" "$WORK/out"
cd "$WORK" || exit 1
"$PYTHON" "$HERE/make_cases.py" "$WORK" || exit 1
BATCHES=$(ls batch_*.json | grep -E "${SWEEP_ONLY:-.}")
N=$(echo "$BATCHES" | wc -l)
export N
: > progress.log

one_batch() {
  b=${1%.json}
  out="$WORK/out/$b"
  mkdir -p "$out"
  # MATLAB: -batch, no display; a session that does not exit after its work
  # is stopped after 90 minutes (outputs are saved as it goes)
  SWEEP_BATCH="$WORK/$1" SWEEP_OUT="$out" env -u DISPLAY timeout 5400 \
    "$MATLAB" -batch "run('$HERE/run_batch.m')" > "$WORK/logs/${b}_ml.log" 2>&1
  status=$?
  "$PYTHON" "$HERE/compare_batch.py" "$WORK/$1" "$out" "$WORK/results_$b.jsonl" \
    > "$WORK/logs/${b}_py.log" 2>&1
  rmdir "$out" 2>/dev/null
  echo "$(date +%H:%M:%S) done $b (MATLAB exit $status)" \
    "- $(grep -c done "$WORK/progress.log") of $N done before this" >> "$WORK/progress.log"
}
export -f one_batch
echo "$BATCHES" | xargs -P "$SWEEP_JOBS" -I{} bash -c 'one_batch {}'
echo "$(date +%H:%M:%S) ALL DONE" >> progress.log
"$PYTHON" "$HERE/summarise.py" "$WORK"
