# Docker

This document explains how to use Docker to build, test, and run semba-fdtd without installing any dependencies on your machine.

## Prerequisite: initialize submodules

Submodules must be initialized on the host before building the image, as the build depends on them:

```bash
git submodule update --init --recursive
```

## Included files

| File | Description |
|---|---|
| `Dockerfile` | Multi-stage build: `builder` (compilation + tests) and `runtime` (binary only) |
| `docker-compose.yml` | Services `solver` (run simulations) and `test` (build and test) |
| `.dockerignore` | Excludes unnecessary files from the build context |

---

## Building the images

```bash
docker compose build test    # image for tests
docker compose build solver  # runtime image for simulations
```

`build` only constructs and saves the image to disk — it does not start any container. You only need to re-run it when the code changes.

### Build arguments

The build mode and optional features can be configured via `--build-arg`:

| Argument | Values | Default |
|---|---|---|
| `BUILD_TYPE` | `Release`, `Debug` | `Release` |
| `ENABLE_MPI` | `ON`, `OFF` | `OFF` |
| `ENABLE_MTLN` | `ON`, `OFF` | `ON` |

```bash
# Debug build
docker compose build --build-arg BUILD_TYPE=Debug test

# Combining arguments
docker compose build \
  --build-arg BUILD_TYPE=Debug \
  --build-arg ENABLE_MPI=ON \
  --build-arg ENABLE_MTLN=OFF \
  test
```

**Release** (`-Ofast`): optimized for speed, no debug information.  
**Debug** (`-g -O0 -fcheck=all -fbacktrace`): no optimization, with runtime checks and backtraces on error — useful for diagnosing crashes.

### Base image digest

The `Dockerfile` pins the base image using a SHA256 digest instead of just the tag:

```dockerfile
FROM ubuntu:22.04@sha256:eb29ed27... AS builder
```

`ubuntu:22.04` is a mutable tag — Canonical can update it at any time. The digest identifies an exact, immutable image, so builds are fully reproducible regardless of when or where they run.

The downside is that OS security patches are not picked up automatically. To update the digest:

```bash
docker pull ubuntu:22.04
docker inspect ubuntu:22.04 --format='{{index .RepoDigests 0}}'
```

Then replace both occurrences of the digest in the `Dockerfile` (builder and runtime stages).

---

## Running the tests

```bash
docker compose run --rm test
```

This runs in sequence:
1. `build/bin/fdtd_tests` — unit tests (GoogleTest)
2. `python3 -m pytest test/ --durations=20` — Python integration tests

To run only part of the test suite:

```bash
# Unit tests only
docker compose run --rm test build/bin/fdtd_tests

# pytest with a specific marker
docker compose run --rm test python3 -m pytest test/ -m mtln
docker compose run --rm test python3 -m pytest test/ -m hdf
```

To open an interactive shell inside the container:

```bash
docker compose run --rm --entrypoint bash test
```

---

## Debugging a simulation

This uses `gdbserver` inside the container and connects VSCode to it via the C/C++ extension.

### Prerequisites

- VSCode extension: **C/C++** (`ms-vscode.cpptools`)
- `gdb` installed on the host:
  ```bash
  sudo apt install gdb
  ```

### Step 1 — Build the debug image

```bash
docker compose build debug
```

### Step 2 — Extract the binary for local symbol loading

VSCode's GDB client needs a local copy of the binary to load debug symbols. Extract it from the image once after each build:

```bash
docker create --name tmp-debug fdtd-debug
docker cp tmp-debug:/src/build/bin/semba-fdtd ./build/bin/semba-fdtd
docker rm tmp-debug
```

### Step 3 — Start gdbserver

Place your input files in `simulations/` and run:

```bash
docker compose run --rm -p 2345:2345 debug case.fdtd.json
```

The container starts and waits, printing something like:

```
Process /src/build/bin/semba-fdtd created; pid = 7
Listening on port 2345
```

### Step 4 — Connect from VSCode

`.vscode/launch.json` is tracked in the repository and already contains the **Docker: attach gdbserver** configuration.

Open the **Run and Debug** panel (`Ctrl+Shift+D`), select **Docker: attach gdbserver**, and press `F5`. VSCode connects, the simulation starts, and breakpoints work normally.

```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Docker: attach gdbserver",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/build/bin/semba-fdtd",
            "miDebuggerServerAddress": "localhost:2345",
            "miDebuggerPath": "gdb",
            "MIMode": "gdb",
            "cwd": "${workspaceFolder}/simulations",
            "sourceFileMap": {
                "/src": "${workspaceFolder}"
            },
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        }
    ]
}
```

> The `sourceFileMap` maps `/src` (container path baked into debug info) to `${workspaceFolder}` on the host, so source files display correctly.

---

## Running a simulation

Place your input files in `simulations/` (created at the repo root) and run:

```bash
docker compose run --rm solver case.fdtd.json
```

The `simulations/` directory is mounted at `/work` inside the container, which is the solver's working directory.

---

## Managing images and containers

### Where are images stored?

Docker manages them internally (under `/var/lib/docker/` on Linux), not in a project folder. To inspect them:

```bash
docker images          # list all saved images
docker image prune     # remove unused images
```

### Viewing active containers

```bash
docker ps      # currently running containers
docker ps -a   # running + stopped containers
```

### Stopping and removing containers

```bash
docker stop <ID_or_name>   # stops the container (does not remove it)
docker rm <ID_or_name>     # removes it
docker rm -f <ID_or_name>  # stop and remove in one step
docker container prune     # remove all stopped containers
```

The full ID is not required — the first 3–4 characters are enough:

```bash
docker stop a1b2
```

> With `--rm` (used in all `docker compose run` commands) the container is removed automatically when it finishes, so no manual cleanup is needed.

### Does closing the terminal stop the container?

- **Without `-d`** (normal mode, what we use): yes, closing the terminal stops the container because it is attached to it.
- **With `-d`** (detached mode): no, it keeps running in the background even after the terminal is closed.

```bash
# Run in the background
docker compose run -d --rm test
```

For tests it is better to omit `-d` so you can see the output in real time.
