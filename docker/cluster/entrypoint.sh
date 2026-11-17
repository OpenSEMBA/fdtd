#!/bin/sh
# Starts sshd so mpirun can launch ranks on this node.
# Host keys are generated per container; user keys come from ./ssh (mounted read-only).
set -e
ssh-keygen -A >/dev/null
mkdir -p /run/sshd /home/mpi/.ssh
cp /cluster-ssh/id_ed25519 /cluster-ssh/id_ed25519.pub /cluster-ssh/authorized_keys /home/mpi/.ssh/
printf 'Host *\n  StrictHostKeyChecking no\n  UserKnownHostsFile /dev/null\n' > /home/mpi/.ssh/config
chown -R 1000:1000 /home/mpi/.ssh
chmod 700 /home/mpi/.ssh
chmod 600 /home/mpi/.ssh/id_ed25519 /home/mpi/.ssh/authorized_keys /home/mpi/.ssh/config
exec /usr/sbin/sshd -D -e
