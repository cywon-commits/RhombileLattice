#!/bin/bash
# Worker pool over a job file (one argument line for annealed_long.py per line).
# Usage: nohup ./run_queue.sh <jobfile> <nworkers> > log 2>&1 &
# Re-launchable after a container restart: finished runs are skipped and
# unfinished ones resume from their checkpoints. Locks live in /tmp (cleared
# by a restart), so each job is claimed by exactly one worker per session.
JOBS=$1; NW=${2:-4}
LOCKDIR=/tmp/rl_locks; mkdir -p $LOCKDIR
worker() {
  while true; do
    got=""
    while read -r line; do
      [ -z "$line" ] && continue
      key=$(echo "$line" | tr ' ' '_')
      if mkdir "$LOCKDIR/$key" 2>/dev/null; then got="$line"; break; fi
    done < "$JOBS"
    [ -z "$got" ] && return
    python3 annealed_long.py $got | tail -n 2
  done
}
for w in $(seq 1 $NW); do worker & done
wait
echo "queue $JOBS finished"
