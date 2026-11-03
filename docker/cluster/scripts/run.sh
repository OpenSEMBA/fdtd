#!/bin/sh
# Runs one case on the cluster: run.sh <case.fdtd.json> [ranks]
# Copies the input file and every file of the case folder that it names (e.g. .exc
# excitations) into work/<case>/, so the outputs stay out of testData/.
set -e

case_file="$1"
ranks="${2:-3}"
[ -n "$case_file" ] || { echo "usage: $0 <case.fdtd.json> [ranks]" >&2; exit 2; }
# Resolve the case before cd, so the path is relative to the caller's directory.
case_file="$(realpath "$case_file")"
cd "$(dirname "$0")/.."

case_dir="$(dirname "$case_file")"
name="$(basename "$case_file" .fdtd.json)"
mkdir -p "work/$name"
cp "$case_file" "work/$name/"
for ref in $(grep -o '"[^"]*"' "$case_file" | tr -d '"'); do
  [ -f "$case_dir/$ref" ] && cp "$case_dir/$ref" "work/$name/"
done

docker compose exec -T node1 su mpi -c "cd /work/$name && mpirun --hostfile /cluster/hostfile -np $ranks --map-by node semba-fdtd -i $(basename "$case_file")"
