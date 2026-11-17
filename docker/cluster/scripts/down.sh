#!/bin/sh
# Stops the cluster and removes the outputs in work/ (keeps SSH keys and .gitkeep).
set -e
cd "$(dirname "$0")/.."
docker compose down -v --remove-orphans
find work -mindepth 1 ! -name .gitkeep -delete
echo "cluster down, work/ emptied"
