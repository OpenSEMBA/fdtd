#!/bin/sh
# Stops a running case: SIGINT to prterun on node1 (OpenMPI 5 launcher), then kills
# any leftover semba-fdtd on every node. Keeps work/ so the checkpoints can be resumed.
set -e
cd "$(dirname "$0")/.."
docker compose exec -T node1 pkill -INT -f 'prterun|mpirun' || true
sleep 3
for n in node1 node2 node3; do
  docker compose exec -T "$n" pkill -9 -f semba-fdtd || true
done
echo "job stopped, work/ kept"
