#!/bin/bash
# Usage: experiment.sh <pid_to_pause> <rank_label>
# Samples both ranks and the point-probe file every 5 s: baseline, SIGSTOP, SIGCONT.
PAUSE=$1
LABEL=$2
PROBE=lossy_cavity.fdtd_ez_Ez_13_9_7_tm.dat
PROBE0=lossy_cavity.fdtd_ez_rank0_Ez_4_5_2_tm.dat
RANK0=$(pgrep -f "[s]emba-fdtd -i" | while read p; do tr '\0' '\n' < /proc/$p/environ | grep -q '^OMPI_COMM_WORLD_RANK=0$' && echo $p; done)
RANK1=$(pgrep -f "[s]emba-fdtd -i" | while read p; do tr '\0' '\n' < /proc/$p/environ | grep -q '^OMPI_COMM_WORLD_RANK=1$' && echo $p; done)
sample() {
  local phase=$1
  for i in 1 2 3 4 5 6; do   # 30 s
    t=$(date +%H:%M:%S)
    s0=$(awk '{print $3}' /proc/$RANK0/stat 2>/dev/null)
    s1=$(awk '{print $3}' /proc/$RANK1/stat 2>/dev/null)
    u0=$(awk '{print $14}' /proc/$RANK0/stat 2>/dev/null)
    u1=$(awk '{print $14}' /proc/$RANK1/stat 2>/dev/null)
    rows=$(wc -l < $PROBE); rows0=$(wc -l < $PROBE0)
    echo "$t phase=$phase rank0(pid=$RANK0 state=$s0 utime=$u0) rank1(pid=$RANK1 state=$s1 utime=$u1) probe_rows_owner1=$rows probe_rows_owner0=$rows0"
    sleep 5
  done
}
echo "# experiment: pause pid $PAUSE ($LABEL)"
sample baseline
kill -STOP $PAUSE; echo "$(date +%H:%M:%S) SIGSTOP sent to $PAUSE"
sample paused
kill -CONT $PAUSE; echo "$(date +%H:%M:%S) SIGCONT sent to $PAUSE"
sample resumed
