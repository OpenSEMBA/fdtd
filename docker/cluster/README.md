# Simulated MPI cluster for semba-fdtd

Three Docker nodes (`node1`..`node3`) running the solver over Open MPI and SSH.
Each node is capped at 1 GiB of RAM with no swap and one CPU.

The solver is compiled from the repository sources with the `rls-mpi` preset
during the image build. The final image contains only the runtime libraries.

## Usage

```sh
./scripts/up.sh      # create SSH keys on first run, build the image, start the nodes
./scripts/down.sh    # stop the nodes and empty work/ (SSH keys are kept)
```

The solver commit is recorded in the binary and printed by `up.sh`.

Run a case with `scripts/run.sh`, which copies the input file to `work/<case>/`
and runs `semba-fdtd` under `mpirun` on the three nodes:

```sh
./scripts/run.sh ../../testData/cases/mpi/bundles_for_mpi.fdtd.json 3
```

## Notes

- `/work` is a bind mount of `docker/cluster/work/`, shared by all nodes. Outputs
  are written there.
- The solver writes its outputs next to the input file, so each case runs in its
  own directory under `work/`.
- The memory limit is enforced by the cgroup. `free` inside a container reports
  host memory, so use `docker stats` to see per-node usage.
- The disk is not capped: `work/` lives on the host filesystem.
