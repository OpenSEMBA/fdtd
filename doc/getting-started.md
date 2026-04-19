# Getting started

This page is the shortest path from clone to first simulation. For platform-specific details and debugging workflows, continue to [](development.md).

## Option 1: use a prebuilt release

Prebuilt binaries are published on the project's GitHub releases page:

- <https://github.com/OpenSEMBA/fdtd/releases>

On Windows, the runtime still depends on Intel oneAPI runtime libraries, as described in [](development.md).

## Option 2: build from source

Before configuring CMake, initialize the submodules from the repository root:

```bash
git submodule update --init --recursive
```

Then configure and build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The main executable is written to `build/bin/semba-fdtd`.

## Run a case

The primary workflow uses an input file in the `.fdtd.json` format:

```bash
build/bin/semba-fdtd -i CASE_NAME.fdtd.json
```

The detailed input reference lives in [](fdtdjson.md). Example inputs are available in `testData/input_examples/`.

## Recommended next steps

- Read [](running-simulations.md) for the normal execution workflow.
- Use [](tutorials/index.md) for worked examples.
- Use [](docker.md) if you prefer containerized build and test workflows.
