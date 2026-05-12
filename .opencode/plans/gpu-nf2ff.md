# GPU Acceleration Plan: Near-to-Far-Field (NF2FF)

## Algorithm Overview

NF2FF computes far-field radiation patterns from near-field E/H data sampled on a Huygens box surface. Two phases:

### Phase 1: Time-Domain DFT Accumulation (`UpdateFarField`)
- Called every `NDecim` timesteps during simulation
- For each of 6 Huygens box faces × 2 field components:
  - Read E/H field values at cell locations on the face
  - Multiply by DFT phase factor `auxExp_E/H(ii)` for each frequency
  - Accumulate into complex buffer arrays `buffer(j,k,ii)`
- **Loop**: `faces × components × Nj × Nk × NumFreqs` complex multiply-adds
- **Data**: Reads from YEE grid fields, writes to 12 complex buffer arrays (3D)
- **Memory**: `12 × Nj × Nk × NumFreqs × 16 bytes` (double-precision complex)

### Phase 2: Far-Field Pattern Computation (`FlushFarfield`)
- Called once at end of simulation (or at flush intervals)
- For each frequency × each (theta, phi) angle pair:
  - For each of 3 face-pairs (Tr/Fr, Iz/De, Ab/Ar):
    - For each face in pair (2 faces):
      - For each cell on Huygens surface:
        - Compute equivalent currents M (magnetic), J (electric)
        - Compute 6 complex phase factors `exp(j·k·r·n̂)`
        - Accumulate L_θ, L_φ, N_θ, N_φ (Huygens integral)
        - Apply PEC/PMC symmetry clones (up to 12 per cell)
  - MPI_AllReduce 4 complex values (across MPI ranks)
  - Compute |E_θ|, |E_φ|, RCS
  - Write output file
- **Loop**: `NumFreqs × Ntheta × Nphi × 3 face-pairs × 2 faces × Ncells × (1 + Nclones)`
- **Core kernel**: `update_LN` — 6 complex exps + 12 complex muls + 8 complex adds per cell

## Current GPU Status

**NONE.** The existing GPU infrastructure (`gpu_state_t` in `gpu_core_m.F90`) only manages YEE fields, CPML buffers, and MUR buffers. No NF2FF state on device.

## GPU Acceleration Strategy

### Priority 1: `FlushFarfield` — Far-Field Pattern Kernel (HIGH IMPACT)

This is the compute-heavy phase with massive parallelism across (freq, theta, phi).

**Kernel Design:**
```
farfield_pattern_kernel<<<grid, block>>>(
    // Device pointers to DFT buffers (read-only)
    ExIz_d, ExDe_d, ExAb_d, ExAr_d, EyFr_d, EyTr_d, EyAb_d, EyAr_d,
    EzIz_d, EzDe_d, EzFr_d, EzTr_d,
    HxIz_d, HxDe_d, HxAb_d, HxAr_d, HyFr_d, HyTr_d, HyAb_d, HyAr_d,
    HzIz_d, HzDe_d, HzFr_d, HzTr_d,
    HxIz2_d, ..., HzTr2_d,  // Schneider averaging buffers

    // Device pointers to geometry (read-only)
    phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d, ...  // 18 coordinate arrays
    dyh_d, dze_d, dye_d, dzh_d,  // cell dimensions

    // Output: far-field results per (freq, theta, phi)
    Etheta_d, Ephi_d, RCS_d,

    // Configuration
    NumFreqs, Ntheta, Nphi,
    thetaStart, thetaStep, phiStart, phiStep,
    freq, comun, cluz, z0,
    // Symmetry flags (bitmask or boolean arrays)
    sym_flags_d
)
```

**Launch Configuration:**
- **Grid**: `(NumFreqs, Ntheta, Nphi)` — each block handles one (freq, theta, phi)
- **Block**: `ceil(Ncells / 2)` threads — each thread handles 1-2 Huygens cells
- **Shared memory**: Precomputed direction cosines, cell positions for current face

**Data Dependencies:**
- Each (freq, theta, phi) is **independent** — no synchronization needed
- MPI reduction: GPU computes per-rank sums, CPU performs `MPI_AllReduce` on 4 complex values per angle (tiny: 128 bytes per angle)
- Two passes (aritmetica/geometrica): sequential within same (freq, theta, phi), handled by kernel

**Symmetry Clones Handling:**
- Current code: `if (flag) call cloneTrFr(...)` — branch-heavy
- GPU strategy: Use **predicated execution** — each thread evaluates all 12 symmetry flags with `if` guards
- Alternative: Precompute symmetry transformation table on host, launch separate kernel passes per symmetry group

**Memory Layout:**
- DFT buffers: Keep as-is `buffer(j, k, ii)` — 2D face access is coalesced along j dimension
- Geometry arrays: Coalesce along j,k dimensions on each face
- Output: `Etheta(freq, theta, phi)` — coalesced along phi dimension

**Estimated Speedup:** 50-200x for typical cases (1000 freq × 100×100 angles, 100×100 Huygens box)

### Priority 2: `UpdateFarField` — DFT Accumulation (MODERATE IMPACT)

This runs every timestep and is more memory-bound. The strided write pattern is the challenge.

**Kernel Design:**
```
dft_accumulate_kernel<<<grid, block>>>(
    // Device pointers to YEE fields (read-only)
    Ex_d, Ey_d, Ez_d, Hx_d, Hy_d, Hz_d,

    // Device pointers to DFT buffers (read-write, atomic-friendly)
    ExIz_d, ExDe_d, ..., HzTr2_d,  // 18 arrays

    // DFT phase factors (read-only, broadcast)
    auxExp_E_d, auxExp_H_d,

    // Face geometry and indices
    face_indices_d,  // (field, j_start, j_end, k_start, k_end) per face/component
    dyh_d, dye_d, dze_d, dzh_d,

    NumFreqs
)
```

**Launch Configuration:**
- **Grid**: `(Nj, Nk)` — each block handles one face cell
- **Block**: `NumFreqs` threads (or use loop if NumFreqs > block size)
- **Alternative**: Grid `(Nj, Nk, NumFreqs)` — each thread handles one (j, k, ii) triple

**Memory Coalescing Challenge:**
- Reads from YEE fields: coalesced along j dimension ✅
- Writes to `buffer(j, k, ii)` with ii innermost: **strided, not coalesced** ❌
- **Mitigation**: Transpose buffer layout to `buffer(ii, j, k)` — makes writes coalesced when threads differ in j/k

**Estimated Speedup:** 5-20x (memory-bound, strided writes)

## Implementation Steps

### Step 1: Extend `gpu_state_t` with NF2FF buffers
```fortran
type gpu_state_t
    ! ... existing fields ...
    
    ! NF2FF DFT buffers (device)
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: ExIz_d, ExDe_d, ExAb_d, ExAr_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: EyFr_d, EyTr_d, EyAb_d, EyAr_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: EzIz_d, EzDe_d, EzFr_d, EzTr_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: HxIz_d, HxDe_d, HxAb_d, HxAr_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: HyFr_d, HyTr_d, HyAb_d, HyAr_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: HzIz_d, HzDe_d, HzFr_d, HzTr_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: HxIz2_d, HxDe2_d, HxAb2_d, HxAr2_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: HyFr2_d, HyTr2_d, HyAb2_d, HyAr2_d
    complex(kind=rkind), pointer, device, dimension(:,:,:) :: HzIz2_d, HzDe2_d, HzFr2_d, HzTr2_d
    
    ! NF2FF frequency arrays (device)
    complex(kind=rkind), pointer, device, dimension(:) :: expIwdt_d, auxExp_E_d, auxExp_H_d
    
    ! NF2FF geometry (device)
    real(kind=rkind), pointer, device, dimension(:) :: phys_x_Mx_d, phys_y_Mx_d, phys_z_Mx_d
    real(kind=rkind), pointer, device, dimension(:) :: phys_x_My_d, phys_y_My_d, phys_z_My_d
    real(kind=rkind), pointer, device, dimension(:) :: phys_x_Mz_d, phys_y_Mz_d, phys_z_Mz_d
    real(kind=rkind), pointer, device, dimension(:) :: phys_x_Jx_d, phys_y_Jx_d, phys_z_Jx_d
    real(kind=rkind), pointer, device, dimension(:) :: phys_x_Jy_d, phys_y_Jy_d, phys_z_Jy_d
    real(kind=rkind), pointer, device, dimension(:) :: phys_x_Jz_d, phys_y_Jz_d, phys_z_Jz_d
    
    ! NF2FF configuration
    logical :: nf2ff_initialized = .false.
    integer(kind=4) :: nf2ff_num_cells
    integer(kind=4) :: nf2ff_num_freqs
end type gpu_state_t
```

### Step 2: Create `gpu_nf2ff_m.F90` with GPU kernels

**Files to create:**
- `src_main_pub/gpu_nf2ff_m.F90` — GPU NF2FF kernels module

**Kernels to implement:**
1. `gpu_init_nf2ff_buffers()` — allocate and copy DFT buffers + geometry to device
2. `gpu_update_nf2ff()` — DFT accumulation kernel (Priority 2)
3. `gpu_flush_nf2ff()` — far-field pattern kernel (Priority 1)
4. `gpu_destroy_nf2ff_buffers()` — deallocate device buffers

### Step 3: Wire into `timestepping.F90`

**At initialization** (after `InitFarField`):
```fortran
call gpu_init_nf2ff_buffers(this%gpu, FF)
```

**During timestep** (in `UpdateFarField` path):
```fortran
if (this%gpu_initialized .and. this%gpu%nf2ff_initialized) then
    call gpu_update_nf2ff(this%gpu, this%bounds)
else
    call UpdateFarField(...)  // CPU fallback
endif
```

**At end of simulation** (in `FlushFarfield` path):
```fortran
if (this%gpu_initialized .and. this%gpu%nf2ff_initialized) then
    call gpu_flush_nf2ff(this%gpu, this%bounds, ...)
    ! D2H transfer results
    call cudaMemcpy(host_results, gpu_results, ...)
else
    call FlushFarfield(...)  // CPU fallback
endif
```

### Step 4: MPI Handling

**Option A (GPU-aware MPI):** Use `MPI_Allreduce` with device pointers if available (NVHPC + CUDA-aware MPI).

**Option B (CPU reduction):** GPU computes per-rank sums, then CPU performs `MPI_AllReduce` on 4 complex values per angle. This is the simpler approach:
- GPU kernel accumulates L_θ, L_φ, N_θ, N_φ per (freq, theta, phi) per rank
- After all angles: copy 4 × NumFreqs × Ntheta × Nphi × 16 bytes to host
- CPU performs `MPI_AllReduce` on these small arrays
- Host computes final |E_θ|, |E_φ|, RCS and writes output

**Recommended:** Option B — the reduction data is tiny (128 bytes per angle), so CPU-side MPI is negligible.

## Data Transfer Analysis

### `UpdateFarField` (during simulation)
- **H2D**: None (fields already on device, buffers stay on device)
- **D2H**: None
- **Bottleneck**: Strided writes to buffer arrays (mitigated by transpose or coalesced access pattern)

### `FlushFarfield` (post-simulation)
- **H2D**: `18 × Nj × Nk × NumFreqs × 16 bytes` — one-time transfer at start
  - Example: 18 × 100 × 100 × 1000 × 16 = **288 MB** (significant but one-time)
- **D2H**: `4 × NumFreqs × Ntheta × Nphi × 16 bytes` — transfer per-rank sums for MPI
  - Example: 4 × 1000 × 100 × 100 × 16 = **64 MB** (one-time at end)
- **Mitigation**: Use persistent device buffers — transfer once at simulation start, not per-timestep

## Performance Estimates

### Sphere case (5 freq, 3 theta × 5 phi = 15 angles, 80×80×80 grid)
- Huygens box: ~80×80 cells per face × 6 faces = ~38,400 cells
- Total work: 5 × 15 × 3 × 2 × 38,400 × 200 FLOPs = **~1.4e10 FLOPs**
- GPU estimate: ~0.5s (vs ~2s current total — NF2FF is ~25% of runtime)

### Conformal sphere (200 freq, 3 theta × 3 phi = 9 angles, larger grid)
- Huygens box: ~200×200 cells per face × 6 faces = ~240,000 cells
- Total work: 200 × 9 × 3 × 2 × 240,000 × 200 FLOPs = **~5.2e12 FLOPs**
- GPU estimate: ~5-10s (vs ~30-60s current — NF2FF is ~50-80% of runtime)

### Production case (1000+ freq, 100×100 angles, large Huygens box)
- Total work: 1000 × 10,000 × 3 × 2 × 100,000 × 200 FLOPs = **~1.2e16 FLOPs**
- GPU estimate: ~500-1000s (vs ~10,000-20,000s current — NF2FF is ~90%+ of runtime)

## Risks & Mitigations

1. **Memory capacity**: 18 complex buffer arrays for large grids may exceed GPU memory
   - Mitigation: Check `cudaMemGetInfo` before allocating; fall back to CPU if insufficient
   
2. **Strided memory access**: DFT buffer writes may not be coalesced
   - Mitigation: Transpose buffer layout to `buffer(ii, j, k)` or use shared memory tiling

3. **MPI complexity**: GPU-aware MPI may not be available
   - Mitigation: Use CPU-side reduction (Option B above)

4. **Symmetry clone branches**: 12 conditional branches per cell may cause warp divergence
   - Mitigation: Predicated execution (NVHPC handles this well) or precompute symmetry groups

5. **Two-pass algorithm**: Geometric + arithmetic passes are sequential
   - Mitigation: Handle both passes in same kernel with a `pasadas` loop variable

## Implementation Order

1. **Step 1**: Extend `gpu_state_t` with NF2FF buffers (1 hour)
2. **Step 2**: Implement `gpu_flush_nf2ff` kernel (Priority 1) — 4-6 hours
3. **Step 3**: Wire `gpu_flush_nf2ff` into `timestepping.F90` (1 hour)
4. **Step 4**: Test with sphere case (1 hour)
5. **Step 5**: Implement `gpu_update_nf2ff` kernel (Priority 2) — 3-4 hours
6. **Step 6**: Wire `gpu_update_nf2ff` into timestep loop (1 hour)
7. **Step 7**: Test with conformal sphere case (2 hours)
8. **Step 8**: Benchmark and optimize (2 hours)

**Total estimated time: 15-18 hours**

## Key Files

| File | Lines | Changes Needed |
|------|-------|----------------|
| `src_main_pub/gpu_nf2ff_m.F90` | 0 → ~800 | **NEW** — GPU NF2FF kernels |
| `src_main_pub/gpu_core_m.F90` | 1,561 | +100 lines — extend `gpu_state_t` |
| `src_main_pub/farfield.F90` | 3,524 | No changes (CPU fallback remains) |
| `src_main_pub/timestepping.F90` | 3,317 | +30 lines — wire GPU NF2FF calls |
| `src_main_pub/observation.F90` | 5,463 | No changes (calls `FlushFarfield`/`UpdateFarField`) |
