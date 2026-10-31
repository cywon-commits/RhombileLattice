#!/bin/bash
# Worker pool over a file of full command lines (env prefix allowed).
# Usage (from the repo root): nohup results/ising_check/run_jobs.sh <jobfile> <nworkers> > log 2>&1 &
# Re-launchable after a container restart: annealed_long.py skips finished
# runs and resumes unfinished ones from their checkpoints. Locks live in /tmp
# (cleared by a restart), so each line is claimed by one worker per session.
JOBS=$1; NW=${2:-4}
LOCKDIR=/tmp/susceptibility_locks; mkdir -p $LOCKDIR
worker() {
  while true; do
    got=""
    while read -r line; do
      [ -z "$line" ] && continue
      key=$(echo "$line" | tr ' =' '__')
      if mkdir "$LOCKDIR/$key" 2>/dev/null; then got="$line"; break; fi
    done < "$JOBS"
    [ -z "$got" ] && return
    bash -c "$got" 2>/dev/null | tail -n 2
  done
}
for w in $(seq 1 $NW); do worker & done
wait
echo "queue $JOBS finished"
