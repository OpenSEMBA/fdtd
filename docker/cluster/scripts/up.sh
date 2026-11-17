#!/bin/sh
# Creates the SSH keys on first run, builds the solver image and starts the cluster.
set -e
cd "$(dirname "$0")/.."

if [ ! -f ssh/id_ed25519 ]; then
    ssh-keygen -t ed25519 -N '' -f ssh/id_ed25519 -q
    cp ssh/id_ed25519.pub ssh/authorized_keys
fi

# The commit is recorded in the binary (see GIT_COMMIT in Dockerfile).
export GIT_COMMIT="$(git -C ../.. rev-parse --short HEAD 2>/dev/null || echo UNKNOWN)"
# The three nodes share one image: build it once, then start them.
docker compose build node1
docker compose up -d

for n in node1 node2 node3; do
    i=0
    until docker compose exec -T "$n" sh -c 'pgrep -x sshd >/dev/null' 2>/dev/null; do
        i=$((i + 1))
        [ "$i" -gt 30 ] && { echo "sshd did not start on $n" >&2; exit 1; }
        sleep 1
    done
done
echo "cluster up: node1 node2 node3 (solver commit $GIT_COMMIT)"
