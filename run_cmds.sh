#!/bin/bash
# Throttled worker pool over a file of shell commands (one per line).
# A job starts only while fewer than MAXP simulation processes
# ("python3 annealed_*") are running machine-wide, so this pool can be
# started while another queue is still draining without oversubscribing.
# Usage: nohup ./run_cmds.sh <cmdfile> [MAXP] > log 2>&1 &
# Re-launchable: every command must itself skip finished work.
CMDS=$1; MAXP=${2:-4}
LOCKDIR=/tmp/rl_cmd_locks; mkdir -p $LOCKDIR
worker() {
  while true; do
    got=""
    exec 9>/tmp/rl_cmd_gate
    flock 9
    while [ "$(pgrep -fc '^python3 annealed_')" -ge "$MAXP" ]; do sleep 20; done
    while read -r line; do
      [ -z "$line" ] && continue
      key=$(echo "$line" | md5sum | cut -c1-16)
      if mkdir "$LOCKDIR/$key" 2>/dev/null; then got="$line"; break; fi
    done < "$CMDS"
    if [ -n "$got" ]; then bash -c "$got" & pid=$!; sleep 5; fi
    flock -u 9
    [ -z "$got" ] && return
    wait $pid
  done
}
for w in $(seq 1 $MAXP); do worker & done
wait
echo "queue $CMDS finished"
