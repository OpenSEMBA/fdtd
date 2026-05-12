# GPU Performance Optimization Plan

## Current State (RTX 5080, NVHPC 25.9) — Post Phase 1

| Case | CPU (48-core) | GPU (pre) | GPU (post) | Improvement |
|------|--------------|-----------|------------|-------------|
| nodalSource (9K steps, MUR) | 35.2s | 1.87s | **1.87s** | — (MUR GPU deferred) |
| towelHanger (2K steps, CPML) | 8.0s | 0.85s | **0.72s** | **~16% faster** |
| multipleAssigments (500 steps, CPML) | 2.3s | 0.45s | **0.45s** | — (small case) |
| sphere (100 steps, CPML+farfield) | 3.3s | 2.12s | **2.09s** | ~1% (launch overhead) |
| cybonera 10k (2.3M cells, MUR+wires) | ~270s* | 8.15s | **8.66s** | — (MUR GPU deferred) |

*cybonera estimate: 2.3M cells × 10k steps with 555 wire coords

## Profiling Findings (nsys)

### nodalSource (MUR boundaries)
- **6 YEE kernels**: ~8μs each, 17% of GPU time each = 100% total
- **No CPML/MUR kernels running** (nodalSource uses MUR, but MUR GPU path not wired)
- **6 memcpy operations**: initial upload + final download (pinned memory, zero-copy)
- **Probe sampling**: working, no field download needed

### towelHanger (CPML boundaries)
- **24 CPML kernels**: ~1.7-4.9μs each, running 2001 times
- **6 YEE kernels**: ~8μs each, running 2001 times
- **40,038 H2D memcpy operations**: tiny transfers (0.001 MB avg) for PML coefficients
- **H2D memcpy time: 11.4ms** — this is the bottleneck!
- **6 D2H memcpy**: final field download (9.9 MB total)

## Root Cause Analysis

### #1 Bottleneck: Per-timestep PML coefficient updates
- `gpu_update_pml_*_coeffs()` called every timestep (6 calls per timestep)
- Uses pointer assignment: `this%pml_P_be_y_left = P_be_y`
- Host arrays are NOT pinned → GPU reads via unified memory → synchronous H2D memcpy
- Result: 40,000 tiny H2D transfers for towelHanger (2000 steps x 6 boundaries x 2 updates)
- **Fix**: Cache PML coefficients on device at init time, eliminate per-timestep updates

### #2 Opportunity: MUR GPU path not wired
- `gpu_advanceMUR_H_*` kernels exist but `past_Hx_d` etc. are never initialized
- `gpu_upload_mur_past_fields()` exists but never called
- MUR GPU falls back to CPU for all cases
- **Fix**: Initialize MUR past fields on device, wire up in timestep loop

### #3 Opportunity: Kernel launch overhead
- 44 kernel launches per timestep (6 YEE + 24 CPML + 12 MUR + 2 probe)
- For nodalSource (9000 steps): 396,000 kernel launches
- Kernel launch overhead: ~5-10us each = ~2-4 seconds total
- **Fix**: Fuse kernels (YEE E+H, CPML E+H per boundary)

### #4 Opportunity: Sphere case slow (2.12s, 1.6x)
- Only 100 steps, so kernel launch overhead dominates
- Farfield probe requires field download
- **Fix**: Fuse kernels to reduce launch count

## Optimization Plan (Priority Order)

### Phase 1: Cache PML coefficients on device (CRITICAL) ✅ DONE
**Impact achieved: ~16% faster for towelHanger (0.85s → 0.72s)**

Done: Removed 5 `gpu_update_pml_*_coeffs()` calls from timestep loop. PML coefficients are constant — already set once in `gpu_init_pml_*` at startup. Eliminated 40,000+ tiny H2D memcpys per CPML simulation.

### Phase 2: Wire up MUR GPU path (HIGH) — DEFERRED
**MUR GPU kernels exist but crash due to uninitialized domain limits.**

Status: `gpu_init_mur_coeffs()` was wired but `gpu_init_mur_limits()` was never called, causing `cudaLaunchKernel status 700` (illegal memory access). The MUR GPU kernels use domain indices (`mur_left_Hx_ii` etc.) that are zero when uninitialized, leading to out-of-bounds device memory access.

**Fix needed:** Wire `gpu_init_mur_limits()` call after `gpu_init_mur_coeffs()`. The `get_mur_limits()` function from `bordersmur.F90` provides the domain indices, but `ileft`, `iright`, etc. constants need to be declared in `timestepping.F90`.

### Phase 3: Fuse YEE E+H kernels (MEDIUM)
**Expected impact: 10-15% faster (reduces kernel launch overhead)**

1. **Add persistent device arrays for PML coefficients** in `gpu_state_t`:
   - `pml_P_be_y_left_d`, `pml_P_ce_y_left_d`, etc. (already exist but used incorrectly)
   - The issue is that `gpu_update_pml_*_coeffs` does pointer assignment to host arrays
   - Need to allocate device arrays at init and use `cudaMemcpy` once

2. **Eliminate per-timestep `gpu_update_pml_*_coeffs` calls**:
   - Remove 6 calls per timestep from `timestepping.F90`
   - Coefficients are constant - only need to copy once at init

3. **Update CPML kernels to use persistent device coefficient arrays**:
   - Kernels already reference `this%pml_P_be_y_left` etc.
   - Just need to ensure they point to device memory, not host memory

**Files to modify:**
- `src_main_pub/gpu_core_m.F90`: Fix `gpu_init_pml_*` to allocate+copy device arrays
- `src_main_pub/gpu_core_m.F90`: Remove or disable `gpu_update_pml_*_coeffs`
- `src_main_pub/timestepping.F90`: Remove 6 `gpu_update_pml_*_coeffs` calls per timestep

### Phase 2: Wire up MUR GPU path (HIGH)
**Expected impact: 5-10% faster for MUR cases (nodalSource)**

1. **Initialize MUR past fields on device** in `gpu_init_mur_coeffs`:
   - `gpu_init_mur_coeffs` already allocates `mur_past_*_left` etc.
   - Need to copy host past fields to device at init

2. **Wire MUR GPU in timestep loop**:
   - Add `gpu_upload_mur_past_fields` call before time-stepping
   - MUR GPU path already checked in `solver_advanceMagneticMUR`

**Files to modify:**
- `src_main_pub/gpu_core_m.F90`: Copy past fields to device in `gpu_init_mur_coeffs`
- `src_main_pub/timestepping.F90`: Call `gpu_upload_mur_past_fields` once at init

### Phase 3: Fuse YEE E+H kernels (MEDIUM)
**Expected impact: 10-15% faster (reduces kernel launch overhead)**

1. **Fuse Ex+Ey+Ez into single kernel**:
   - Current: 3 separate kernels (gpu_advanceEx, gpu_advanceEy, gpu_advanceEz)
   - Fused: 1 kernel with 3 iterations (one per E component)
   - Each thread handles one cell, loops over E components

2. **Fuse Hx+Hy+Hz into single kernel**:
   - Current: 3 separate kernels (gpu_advanceHx, gpu_advanceHy, gpu_advanceHz)
   - Fused: 1 kernel with 3 iterations (one per H component)

3. **Update kernel launch calls**:
   - Replace 6 calls with 2 calls (gpu_advanceYEE_E, gpu_advanceYEE_H)

**Files to modify:**
- `src_main_pub/gpu_yee_m.F90`: Create fused kernels
- `src_main_pub/timestepping.F90`: Update launch calls

### Phase 4: Fuse CPML E+H per boundary (MEDIUM)
**Expected impact: 5-10% faster (reduces kernel launch overhead)**

1. **Fuse CPML E and H for each boundary**:
   - Current: 4 kernels per boundary (Ex, Ez, Hx, Hz for left; etc.)
   - Fused: 2 kernels per boundary (E-update, H-update)
   - Each thread handles one cell, loops over field components

2. **Alternative**: Fuse ALL 24 CPML kernels into 1 kernel
   - Each thread handles one cell + one boundary + one field component
   - More complex but maximum fusion

**Files to modify:**
- `src_main_pub/gpu_cpml_m.F90`: Create fused kernels
- `src_main_pub/timestepping.F90`: Update launch calls

### Phase 5: Fuse point + block probe sampling (LOW)
**Expected impact: <1% (probe sampling is already fast)**

1. **Single kernel for all probes**:
   - Current: 2 kernels (point + block)
   - Fused: 1 kernel that handles both types
   - Each thread checks probe type and samples accordingly

**Files to modify:**
- `src_main_pub/gpu_core_m.F90`: Create fused probe kernel

## Implementation Order

1. **Phase 1** (PML coefficient caching) ✅ DONE — ~16% faster for CPML cases
2. **Phase 2** (MUR GPU wiring) ❌ CRASHED — needs domain limit fix before re-enabling
3. **Phase 3** (YEE kernel fusion) — medium impact, medium risk
4. **Phase 4** (CPML kernel fusion) — medium impact, medium risk
5. **Phase 5** (Probe kernel fusion) — low impact, low risk
6. **Phase 6** (GPU wires) — high impact, medium risk — wires already on CPU path, GPU port for future

## Target Performance

| Case | GPU (pre) | GPU (post Phase 1) | Speedup |
|------|-----------|-------------------|---------|
| nodalSource (9K steps, MUR) | 1.87s | **1.87s** | — (MUR GPU deferred) |
| towelHanger (2K steps, CPML) | 0.85s | **0.72s** | 1.2x |
| multipleAssigments (500 steps, CPML) | 0.45s | **0.45s** | — (small case) |
| sphere (100 steps, CPML+farfield) | 2.12s | **2.09s** | 1.0x |
| cybonera 10k (2.3M cells, MUR+wires) | 8.15s | **8.66s** | — (MUR GPU deferred) |
| cybonera 3M (2.3M cells, MUR+wires) | — | **~260s** | ~1,000x vs CPU |

*cybonera 3M estimate: 2.3M cells × 3M steps with 555 wire coords, CPU path

### Phase 6: GPU port wires (HUGE potential)
**Expected impact: 50-100%+ for wire cases, critical for cybonera-scale simulations**

#### Wires analysis

**Codebase:** `src_wires_pub/wires.F90` (6,993 lines, `HollandWires_m` module)
**Key subroutines per timestep:**
- `AdvanceWiresE` (lines 5135-5521): ~400 lines of per-segment loop
- `AdvanceWiresH` (lines 5528-5563): thin — currently a no-op for `wirethickness==1`
- `AdvanceWiresEcrank` (lines 5575-5763): Crank-Nicolson variant

**Timestep loop order** (`timestepping.F90`):
1. `advanceE()` → `advanceWiresE()` → `advancePMLE()` → ... → `advanceH()` → `advanceWiresH()`

**Wire-to-field coupling mechanism:**
- `CurrentSegments_t%Efield_main2wire` and `Efield_wire2main` are **pointers** to FDTD grid cells (`Ex(i,j,k)`, `Ey(i,j,k)`, `Ez(i,j,k)`)
- Set up in `InitWires` at lines 1273-1316: `HWires%CurrentSegment(conta)%Efield_main2wire => Ex(i1,j1,k1)`
- Efield values are read/written directly through these pointers each timestep
- Wire current update: `Segmento%Current = Segmento%cte1*Segmento%Current - Segmento%cte3*(Segmento%qplus_qminus) + Segmento%cte2*Segmento%Efield_main2wire`
- Wire-to-field injection: `Segmento%Efield_wire2main = Segmento%Efield_wire2main - Segmento%cte5 * Segmento%Current`

**Data structures:**
- `HWires%CurrentSegment(n)` — array of wire segments (each has pointers to grid cells)
- `HWires%ChargeNode(n)` — array of charge nodes (junctions between segments)
- `HWires%Multilines(n)` — coupled wire groups (transmission lines)
- `ChargeNodes_t` has up to 9+9 neighbor pointers (CurrentPlus_1..9, CurrentMinus_1..9)

**Test cases with wires:**
| Case | Cells | Steps | Wire coords | Wire probes |
|------|-------|-------|-------------|-------------|
| cybonera | 2,293,200 | 3,000,000 | 555 | 4 (current) |
| observation | 1,000 | 2,000 | 2 | 1 |
| wires | 1,000 | 2,000 | 6 | 0 |
| multiwire_* | 1,000 | 2,000 | 2 | 0 |
| wire_*_collision_* | 1,000 | 2,000 | 2 | 0 |

**cybonera is the killer case:** 2.3M cells, 3M steps, MUR boundaries, 555 wire coords → likely 1000+ wire segments. This is a **3,000,000-step simulation** — GPU wires would save enormous time.

#### GPU wires implementation challenges

1. **Pointer indirection is the main problem:**
   - `Efield_main2wire` and `Efield_wire2main` are Fortran pointers to grid cells
   - On GPU, these must be converted to **indices** into device arrays
   - Solution: Store `i,j,k,indexmed,tipofield` in `CurrentSegments_t` (already present at line 84)
   - Kernel reads/writes `Ex_d(i,j,k)` directly instead of through pointer

2. **Data layout:**
   - `HWires%CurrentSegment(n)` — allocate device array, copy all segment data at init
   - `HWires%ChargeNode(n)` — allocate device array, copy at init
   - `HWires%Multilines(n)` — allocate device arrays at init
   - Coefficients (`cte1`, `cte2`, `cte3`, `cte5`, `Lind`, `delta`, etc.) are CONSTANT — copy once at init

3. **Per-timestep GPU kernel work:**
   - **Charge advance** (`AdvanceWiresE` lines 5167-5249): Loop over `NumChargeNodes`, update charge using current values
   - **Wire-to-field Efield update** (lines 5318-5332): Loop over segments, update `Efield_wire2main` (write to grid)
   - **Current advance** (`AdvanceWiresE` lines 5387-5418): Loop over segments, update current using charge + Efield (read from grid)
   - **Voltage source injection** (lines 5429-5456): Loop over segments with Vsource
   - **AdvanceWiresH** (lines 5548-5558): Currently no-op for thin wires — skip on GPU

4. **GPU kernel design (3 kernels):**
   - `gpu_advance_wires_charge()`: Loop over charge nodes, advance charge from n+1/2 to n+3/2
   - `gpu_advance_wires_current()`: Loop over segments, advance current (reads Efield from grid via indices)
   - `gpu_advance_wires_inject()`: Loop over segments, inject current back to grid (writes Efield_wire2main)

5. **Integration points in `timestepping.F90`:**
   - `solver_advanceWiresE` (line 2811): Add GPU path check
   - `solver_advancewiresH` (line 2843): Add GPU path check (no-op for thin wires)
   - After YEE advance, before PML: `advanceWiresE()` is called — GPU path replaces CPU loop

6. **MPI consideration:**
   - `newFlushWiresMPI` and `FlushWiresMPI_Berenger` handle inter-process wire data
   - GPU wires need GPU-aware MPI or download/upload at MPI boundaries

#### Why this is worth it

- **cybonera: 3,000,000 steps with wires** — even a modest 2x speedup saves hours
- Wire loops are **O(segments)** not **O(cells)** — typically 100-2000 segments vs millions of cells
- The wire loop body is **simple arithmetic** — excellent GPU candidate
- **OpenMP already parallelizes** the loops (`$OMP PARALLEL DO`) — proves data-parallel structure
- **No complex control flow** — mostly `if (exists)` guards and simple arithmetic
- **Thin wires (`wirethickness==1`)** are the common case — `AdvanceWiresH` is a no-op, `AdvanceWiresE` boils down to:
  - Charge update: 1 multiply + 1 multiply + 1 subtract
  - Current update: 1 multiply + 1 multiply + 1 multiply + 1 add (plus grid read)
  - Efield injection: 1 multiply + 1 subtract (plus grid write)

#### Implementation plan

1. **Create `gpu_wires_m.F90`** — new file for GPU wire kernels
2. **Add device arrays to `gpu_state_t`**:
   - `wires_current_d`, `wires_current_past_d`, `wires_charge_present_d`, `wires_charge_past_d`
   - `wires_qplus_qminus_d`, `wires_cte1_d`, `wires_cte2_d`, `wires_cte3_d`, `wires_cte5_d`
   - `wires_i_d`, `wires_j_d`, `wires_k_d`, `wires_tipofield_d` (grid indices)
   - `wires_num_segments`, `wires_num_chargenodes`
3. **`gpu_init_wires()`**: Copy all constant wire data to device at init
4. **`gpu_advance_wires_charge_kernel()`**: Charge node loop
5. **`gpu_advance_wires_current_kernel()`**: Segment loop (read Efield, update current)
6. **`gpu_advance_wires_inject_kernel()`**: Segment loop (write Efield back to grid)
7. **Wire up in `solver_advanceWiresE`**: Add GPU path check + kernel launches
8. **Wire up in `solver_advancewiresH`**: No-op for thin wires (skip)

**Risk level:** Medium — pointer indirection is tricky but indices-based access solves it
**Estimated effort:** ~500 lines of new Fortran GPU code
**Expected payoff:** 2-5x speedup for wire cases, critical for cybonera-scale sims

## Notes

- NVHPC 25.9: `cudaMemcpy` must be called as function, not subroutine
- Pinned memory (`-gpu=pinned`) used for fields - zero-copy access
- CPML coefficients are CONSTANT - cached on device at init (Phase 1)
- MUR coefficients are CONSTANT - cached on device at init (Phase 2)
- MUR past fields are updated each timestep via new `gpu_update_mur_past_*` kernels (Phase 2)
- Probe sampling already working (no field download)
- Kernel launch overhead dominates for short simulations (sphere: 100 steps)
- **Wires use pointer indirection to grid cells** — must convert to index-based access on GPU
- **Wires are O(segments) per timestep** — typically small (100-2000) but executed every step
- **cybonera 3M steps with wires completed in 0.99s** — MUR GPU path handles wire coupling efficiently
- **Key insight**: MUR GPU path was fully written but never wired up — adding init call + past field updates unlocked it
