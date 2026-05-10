---
name: hpc-engineer
description: 1D Z-axis MPI domain decomposition, OpenMP collapse(2) on YEE kernels, compiler optimization flags (GCC/Intel), profiling with gprof and NVTX placeholders, performance bottleneck analysis, and GPU roadmap with OpenACC placeholders
---

## When to use

- Optimizing simulation performance for large 3D domains
- Configuring MPI domain decomposition for cluster runs
- Tuning OpenMP parallelization for multi-core nodes
- Selecting compiler optimization flags for target hardware
- Profiling and identifying performance bottlenecks
- Evaluating GPU acceleration strategies (OpenACC roadmap)
- Diagnosing MPI communication overhead or load imbalance

## Key Files to Read

### MPI Parallelization
- `src_main_pub/mpicomm.F90` (1801 lines) — Full MPI infrastructure: domain decomposition, field exchange, wire communication
- `src_main_pub/timestepping.F90` — Time-stepping loop with MPI flush calls

### OpenMP Kernels
- `src_main_pub/timestepping.F90` — Core YEE update kernels with `!$OMP PARALLEL DO COLLAPSE(2)`
- `src_wires_pub/wires.F90` — Thin-wire model parallel loops
- `src_main_pub/planewaves.F90` — Plane wave source update loops (100+ OpenMP references)

### Build Configuration
- `CMakeLists.txt` — Compiler flags, optional features, optimization options
- `.github/workflows/ubuntu.yml` — CI build matrix with compiler variants
- `set_precompiled_libraries.cmake` — Prebuilt library configuration

### Profiling
- `testData/cases/paul/gprof_output.txt` — Example gprof profiling output
- `testData/cases/holland/gprof_paul.txt` — Example gprof profiling output
- `testData/cases/multilines_opamp/gprof_output.txt` — Example gprof profiling output

## MPI Parallelization

### Domain Decomposition Strategy

The project uses **1D domain decomposition along the Z-axis** (z-slicing). The computational domain is divided into horizontal slices, each assigned to one MPI process.

**Key functions in `src_main_pub/mpicomm.F90`:**

| Function | Line | Purpose |
|----------|------|---------|
| `InitGeneralMPI()` | 71-81 | Initialize MPI, get communicator size and rank |
| `MPIdivide()` | 86-304 | Divide domain into Z-slices, assign sub-domains with PML handling |
| `InitMPI()` / `InitMPI_Cray()` | 310, 1238-1423 | Initialize communication buffers for field exchange |
| `FlushMPI_H()` / `FlushMPI_E()` | 434-587 | Exchange H/E field ghost cells (non-blocking I/O) |
| `FlushMPI_H_Cray()` / `FlushMPI_E_Cray()` | 1425-1593 | Cray-optimized flush variants |
| `MPIupdateMin()` | 361-367 | MPI_AllReduce with MPI_MIN for time-step synchronization |
| `MPIupdateBloques()` | 372-380 | MPI_AllReduce with MPI_SUM for field observation |
| `newInitWiresMPI()` / `newFlushWiresMPI()` | 600, 869 | Wire current/charge node data exchange |

### MPI Communication Pattern

Each time step performs **synchronous boundary exchange** between adjacent slices:

    ! In solver_run / step (timestepping.F90):
    call advanceE_fields(...)
    call FlushMPI_E_Cray(...)    ! Exchange E-field ghost cells

    call advanceH_fields(...)
    call FlushMPI_H_Cray(...)    ! Exchange H-field ghost cells

**Communication details:**
- Uses **non-blocking I/O** (`MPI_ISEND`/`MPI_IRECV`) with `MPI_WAITALL` for overlap
- Exchanges ghost cells between **adjacent ranks only** (neighbors in Z direction)
- Wire models require additional MPI communication for current/charge node data at slice boundaries
- Anisotropic materials, SGBC, and multiport BCs need extra flush buffers (`InitExtraFlushMPI`)

### MPI Command-Line Usage

```bash
# Basic MPI run
mpirun -n 4 semba-fdtd -i my_case.fdtd.json

# With hostfile for cluster
mpirun -n 64 -hostfile hosts semba-fdtd -i large_case.fdtd.json

# Via Python pyWrapper
solver = FDTD(
    input_filename='my_case.fdtd.json',
    path_to_exe='/path/to/semba-fdtd',
    mpi_command='mpirun -n 8'
)
solver.run()
```

## OpenMP Parallelization

### Always Enabled

OpenMP is **unconditionally compiled** in all builds. The `CompileWithOpenMP` define is always set.

### Compiler Flags

| Compiler | OpenMP Flag |
|----------|-------------|
| GNU (gfortran) | `-fopenmp` |
| IntelLLVM (ifx) | `-qopenmp` |

### Core YEE Update Kernels

All six field update kernels in `timestepping.F90` use `COLLAPSE(2)` on the inner j,k loops:

| Kernel | Line | OpenMP Directive |
|--------|------|-----------------|
| `advanceEx` | 2206-2225 | `!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)` |
| `advanceEy` | 2249-2266 | `!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)` |
| `advanceEz` | 2295-2312 | `!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)` |
| `advanceHx` | 2358-2376 | `!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)` |
| `advanceHy` | 2400-2417 | `!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)` |
| `advanceHz` | 2439-2455 | `!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2)` |

**Parallelization strategy:** `COLLAPSE(2)` on the j,k loops leaves the k-loop as the outer OpenMP dimension. This is well-suited for 3D FDTD where Y and Z dimensions are typically larger than X.

### Other Parallelized Loops

- **Wire charge node updates** (`wires.F90`): `!$OMP PARALLEL DO` on wire charge nodes
- **Plane wave sources** (`planewaves.F90`): Many loops with `!$OMP PARALLEL DO`
- **Observation FFT** (`observation.F90`): Frequency-domain probe computation
- **Tag computation** (`timestepping.F90` `fillMtag()`): Three separate `!$OMP PARALLEL DO` regions

### OpenMP Configuration

```bash
# Set number of threads per MPI rank
export OMP_NUM_THREADS=4
mpirun -n 8 ./build/bin/semba-fdtd -i case.fdtd.json
# Total: 8 MPI ranks x 4 OpenMP threads = 32 cores

# Hybrid parallelism (recommended for multi-node clusters)
# Node 1: 2 MPI ranks x 12 OpenMP threads = 24 threads
# Node 2: 2 MPI ranks x 12 OpenMP threads = 24 threads
export OMP_NUM_THREADS=12
mpirun -n 4 ./build/bin/semba-fdtd -i case.fdtd.json
```

## Compiler Optimization Flags

### GNU Compiler (gcc/gfortran)

| Mode | Flags |
|------|-------|
| Release | `-Ofast` |
| Debug | `-g -O0 -fno-inline -fcheck=all -fbacktrace` |
| Common | `-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch` |

### IntelLLVM Compiler (ifx)

| Mode | Flags |
|------|-------|
| Release | `-O3 -fp-model fast=2w` |
| Debug | `-check all,nouninit -debug full -traceback` |
| Common | `-qopenmp -fpp -static-intel` |

### Optional Intel Optimizations

| CMake Option | Flag | Purpose |
|--------------|------|---------|
| `SEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION` | `-xHost` | CPU-specific auto-vectorization |
| `SEMBA_FDTD_ENABLE_INTEL_IPO` | `-ipo` | Interprocedural optimization |

### Build Commands for Performance

```bash
# Standard release build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# With Intel CPU-specific optimization
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DSEMBA_FDTD_ENABLE_INTEL_XHOST_OPTIMIZATION=ON \
  -DSEMBA_FDTD_ENABLE_INTEL_IPO=ON
cmake --build build -j

# Double precision (may be slower but more accurate)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION=ON
cmake --build build -j
```

## Profiling Infrastructure

### gprof (Historical)

The project has been profiled with gprof in the past. Reference outputs exist in:
- `testData/cases/paul/gprof_paul.txt`
- `testData/cases/holland/gprof_output.txt`
- `testData/cases/multilines_opamp/gprof_output.txt`

To profile with gprof, rebuild with `-pg` flag and run the simulation:

```bash
# Add -pg to Fortran flags in CMakeLists.txt temporarily
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run simulation
./build/bin/semba-fdtd -i case.fdtd.json

# Analyze
gprof ./build/bin/semba-fdtd gmon.out > profile.txt
```

### NVTX (NVIDIA Tools Extension — Placeholder)

The code contains NVTX range markers under `#ifdef CompileWithProfiling` in `timestepping.F90` (15 locations):

    #ifdef CompileWithProfiling
    call nvtxStartRange("Antes del bucle N")       ! Entire time loop
    call nvtxStartRange("Antes del bucle EX/EY/EZ") ! E-field components
    call nvtxStartRange("Antes del bucle HX/HY/HZ") ! H-field components
    call nvtxEndRange
    #endif

However, `CompileWithProfiling` is never defined in the build system, and the `nvtx` module does not exist in the repository. This is **forward-looking infrastructure** for GPU profiling with NVIDIA tools.

## Performance Bottlenecks

### Primary Bottlenecks

1. **MPI communication** — Synchronous boundary exchange at every time step between adjacent Z-slices. This is the dominant bottleneck for large MPI runs.
2. **OpenMP thread scaling** — `COLLAPSE(2)` works well for large Y,Z dimensions but may underutilize threads for thin domains.
3. **Wire model overhead** — Thin-wire models add sequential dependencies (charge node updates) that limit parallelism.
4. **Memory bandwidth** — FDTD is memory-bound; field arrays are large and accessed repeatedly each time step.

### Optimization Strategies

| Strategy | Expected Impact | Effort |
|----------|----------------|--------|
| Increase MPI ranks (more Z-slices) | High for large domains | Low |
| Tune `OMP_NUM_THREADS` per rank | Medium | Low |
| Enable `-xHost` (Intel) | Medium (vectorization) | Low |
| Enable IPO (`-ipo`) | Medium (inlining) | Low |
| Use double precision only when needed | High (memory bandwidth) | Low |
| 2D/3D MPI decomposition | High (less communication) | High |
| GPU acceleration (OpenACC) | High | Very high |
| Loop unrolling / SIMD hints | Low-Medium | Medium |

## GPU Roadmap

### Current State

**No active GPU acceleration.** The codebase has:
- No CUDA, OpenCL, ROCm, or SYCL code
- No NVIDIA HPC compiler in the active CI matrix
- `CompileWithACC` placeholders in `timestepping.F90` (6 locations) suggesting **OpenACC planning**

### OpenACC Placeholders (Inactive)

In `timestepping.F90`, the YEE kernels have commented-out OpenACC directives:

    #ifdef CompileWithACC
    !$ACC parallel loop DEFAULT(present) collapse(2) \
        private(...) copyin(...) copyout(...)
    #endif

These exist for all six field update kernels (lines 2209, 2252, 2298, 2361, 2403, 2442).

### Recommended GPU Strategy

1. **Start with OpenACC** — The existing placeholders map directly to the YEE kernels. Enable with NVHPC or NVIDIA HPC SDK compiler.
2. **Target the E/H field updates** — These are the most parallelizable sections with `COLLAPSE(2)` loop structure.
3. **Handle MPI separately** — GPU-aware MPI (UCX, NCCL) for inter-node communication.
4. **Wire models are harder** — Sequential charge node updates may need algorithmic changes for GPU.

### NVHPC Compiler Support (Stub)

NVHPC is partially supported in `CMakeLists.txt` but marked as `TODO: Tune flags for NVHPC`. The compiler detection exists but optimization flags are not configured.

## Performance Tuning Checklist

### Before Running Large Simulations

- [ ] Verify CFL condition is satisfied (check log for adjusted dt)
- [ ] Set `OMP_NUM_THREADS` appropriate for node core count
- [ ] Use MPI rank count matching available nodes/cores
- [ ] Enable compiler optimizations (`-Ofast` for GCC, `-O3` for Intel)
- [ ] Consider `-xHost` for Intel CPUs (if geometry is large enough)
- [ ] Use single precision unless accuracy requires double precision
- [ ] Disable unused features (MTLN, HDF, SMBJSON) to reduce memory

### For MPI Scaling Tests

- [ ] Run strong scaling test: fixed problem size, increasing MPI ranks
- [ ] Run weak scaling test: fixed work per rank, increasing total size
- [ ] Monitor MPI communication time in logs
- [ ] Check for load imbalance (some ranks finishing much earlier)
- [ ] Verify Z-slice decomposition matches aspect ratio of domain

### Debugging Performance Issues

- [ ] Check `SEMBA_FDTD_temp.log` for timing information per time step
- [ ] Profile with gprof to identify hot functions
- [ ] Verify OpenMP thread count with `OMP_NUM_THREADS`
- [ ] Check for unnecessary data copying or reallocation in hot loops
- [ ] Ensure contiguous array allocation (check `contiguous` attribute on field pointers)
- [ ] Review pointer aliasing patterns for potential false dependencies
