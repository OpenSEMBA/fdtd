# Plane-wave TF/SF and Mur: Fortran to C++ slim mapping

Active C++ path: [`semba_fdtd.cpp`](../src_cpp/src_main_pub/semba_fdtd.cpp) (`FDTD_Solver` in `semba-fdtd-core`).

Fortran reference: [`planewaves.F90`](../src_main_pub/planewaves.F90), [`bordersmur.F90`](../src_main_pub/bordersmur.F90), [`timestepping.F90`](../src_main_pub/timestepping.F90) `step()`.

Dormant fuller port (not linked in slim build): [`planewaves.cpp`](../src_cpp/src_main_pub/planewaves.cpp), [`bordersmur.cpp`](../src_cpp/src_main_pub/bordersmur.cpp).

## Time loop order

| Order | Fortran `step()` | C++ `launch()` |
|-------|------------------|----------------|
| 0 | `flushPlanewaveOff` | `flushPlanewaveOff()` |
| 1 | `advanceE` | `advanceE()` |
| 2 | `AdvancePlaneWaveE` | `advancePlaneWaveE()` (gated) |
| 3 | `advanceH` | `advanceH()` |
| 4 | `AdvancePlaneWaveH` | `advancePlaneWaveH()` (gated) |
| 5 | `AdvanceMagneticMUR` | `applyMurH()` |

There is **no** Fortran electric Mur. `applyMurE()` in C++ is unused; tangential E at domain faces is updated in `advanceE()`.

## Function mapping

### `evolucion` — planewaves.F90 L783–810

| Fortran | C++ slim |
|---------|----------|
| `nprev = int((t-d/cluz)/deltaevol)` | `nprev = floor(t_delay / deltaevol)` |
| `(nprev+1 <= numus)` | `(nprev + 1 <= numSamples)` |
| `nprev > 0` | `nprev >= 1` |
| `evol(jjj,nprev)`, `evol(jjj,nprev+1)` (1-based) | `samples[nprev-1]`, `samples[nprev]` (0-based) |

### `Incid` — planewaves.F90 L737–812

| Fortran | C++ slim |
|---------|----------|
| `Punto%PhysCoor` | `physCoord1()` |
| `d = x*px+y*py+z*pz - distanciaInicial` | same in `computeIncid()` |
| `fpw(jjj,nfield,k)*evolucion(...)` | `pw.ex/hx/... * evolucion(...)` |

### `InitPlaneWave` — planewaves.F90 L38–729

| Fortran | C++ slim |
|---------|----------|
| Direction/polarization unit vectors | `initPlaneWave()` L501–515 |
| H from E×k/Z₀ | L513–515 |
| Huygens box `esqx1..esqz2` | L522–544 |
| `IluminaTr/Fr/...` | L545–550 |
| `distanciaInicial` octant | L552–587 |

### `AdvancePlaneWaveE` — planewaves.F90 L837–1133

| Face | Fortran flag | C++ slim |
|------|--------------|----------|
| x− | `IluminaTr` | `pw.iluminaTr` L833+ |
| x+ | `IluminaFr` | `pw.iluminaFr` |
| y− | `IluminaIz` | `pw.iluminaIz` |
| y+ | `IluminaDe` | `pw.iluminaDe` |
| z− | `IluminaAb` | `pw.iluminaAb` L905+ |
| z+ | `IluminaAr` | `pw.iluminaAr` L923+ |

Coefficients: `G2_1 = g2(1) = dt/eps0`, `Id = 1/dx` (uniform).

### `AdvancePlaneWaveH` — planewaves.F90 L1139–1434

Same face flags; time `sgg%tiempo(n) + 0.5*dt` → `currentTime + 0.5*dt`. Coefficient `Gm2_1 = dt/mu0`.

### `calc_murconstants` / `AdvanceMagneticMUR` — bordersmur.F90 L294–341, L1107+

| Fortran | C++ slim |
|---------|----------|
| `CAB1 = (1-CNUM)/(1+CNUM)`, `CNUM = (1/dx)/(dt*cluz/sqrt(Epr*Mur))` | `murCx = (C0*dt-dx)/(C0*dt+dx)` (vacuum uniform) |
| `H_bnd = H_past_int + CAB1*(H_int - H_past_bnd)` | `mur_face()` lambda L634 |
| Past-field store L1374+ | L741–793 |

## Known parity gap (pw-in-box)

The slim solver uses **interior-only** backward/forward FDTD stencils and **H-face-only** TF/SF in `launch()` (no `advancePlaneWaveE` in the time loop). Probe parity passes through ~90 steps (`MediumRunProbeParity_First90Steps`) but the full 1298-step run diverges after the excitation peak (~step 100+).

Fortran pairs full-domain `advanceEx/Ey/Ez` with six-face `AdvancePlaneWaveE/H`. Re-enabling E-face TF/SF with partial `advanceE` over-injects (~1.85×); full ghost-cell sweeps need Fortran-matched boundary indexing.

## GDB workflow

Debug build:

```bash
cmake -S . -B cpp_build_dbg \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=OFF \
  -DCMAKE_BUILD_TYPE=Debug
cmake --build cpp_build_dbg -j --target semba-fdtd-cpp cpp_tests
```

Or run `./scripts/debug_pw_in_box_gdb.sh`.

### Breakpoints (slim solver)

| Location | Purpose |
|----------|---------|
| `semba_fdtd.cpp:evolucion` | Excitation interpolation / off-by-one |
| `semba_fdtd.cpp:computeIncid` | Incident field at probe (3,3,3) |
| `semba_fdtd.cpp:advancePlaneWaveE` | TF/SF E injection on z− face |
| `semba_fdtd.cpp:advancePlaneWaveH` | TF/SF H at half step |
| `semba_fdtd.cpp:applyMurH` | Domain H Mur |
| `semba_fdtd.cpp:sampleProbes` | Probe output vs reference `.dat` |

### GDB checklist

1. After init: `p solver.dt`, `p solver.numSteps` — expect CFL dt ≈ 1.54066656123717649e-11, steps ≈ 1300.
2. Step 1: compare `computeIncid(0,0,t,3,3,3)` with inbox `.dat` `incid` column.
3. After `advancePlaneWaveE`: `Ex` at TF face z=2 (`ex_idx(2,2,1)`).
4. After `applyMurH`: boundary `Hy[hy_idx(0,j,k)]` vs interior neighbor.
5. Optional: run Fortran `build/bin/semba-fdtd` on same case at matching step and diff fields.

### Example session

```bash
cd testData/cases/planewave
gdb -ex 'break semba_fdtd.cpp:485' \
    -ex 'break semba_fdtd.cpp:838' \
    -ex 'run --args ../../../cpp_build_dbg/bin/semba-fdtd-cpp -i pw-in-box.fdtd.json'
```

At breakpoint: `print computeIncid(0, 0, currentTime, 3, 3, 3)` (if in scope via test hook or manual field inspection).

## GoogleTest coverage

| Test file | Fortran lines verified |
|-----------|------------------------|
| `test_planewave_evolucion.h` | evolucion L793–798 |
| `test_planewave_init.h` | InitPlaneWave L331–371, calc_planewaveconstants |
| `test_planewave_tfsf.h` | AdvancePlaneWaveE/H Ab/Ar faces |
| `test_bordersmur.h` | calc_murconstants L310, AdvanceMagneticMUR L1107, pulse absorption integration |
| `test_planewave_pw_in_box.h` | Full `pw-in-box` probe parity |

## Mur absorption test (phase 1)

Case: [`testData/cases/mur/pulse-1d-x.fdtd.json`](../testData/cases/mur/pulse-1d-x.fdtd.json) — 6×6×6 vacuum, Mur on all faces, CFL-safe `timeStep` = 1.54066656123717649e-11 (matches `heurCFL=0.8`).

Hook: `SEMBA_FDTD_test::test_mur_pulse_absorption()` seeds a single-cell Ex pulse, runs `advanceE → advanceH → applyMurH` (no plane-wave TF/SF), returns peak/probe/energy metrics. Pass `apply_mur=false` to compare against open boundaries (no Mur overwrite).

GoogleTests in `test_bordersmur.h`:

| Test | Intent |
|------|--------|
| `PulsePeakDecaysAfterTransit` | After N ≈ 2.5·NX·dx/(c·dt) steps, Mur peak Ex < 25% of open-boundary peak |
| `CenterProbeNearZeroAfterAbsorption` | Center probe with Mur ≤ open-boundary probe |
| `EnergyDecreasesAfterPeak` | Total E energy with Mur < open-boundary energy |
| `MurReducesPeakVersusOpenBoundary` | Same comparison explicitly |
| `PulseMatchesFortranProbe` | **Optional** — skipped unless golden `.dat` exists |

Strict absolute decay (`max_ex_final < 0.05·max_ex_initial`) requires full-domain FDTD sweeps (phase 3); the slim solver uses interior-only `advanceE/H` for `launch()` parity with pw-in-box.

### Fortran line map (first-order magnetic Mur)

| Face | Fortran `AdvanceMagneticMUR` | C++ slim `applyMurH()` |
|------|------------------------------|-------------------------|
| Back (x−) Hy | L1305–1312 | Hy at `i=0`, `Past_Hy(i+1)` formula |
| Back Hz | L1287–1294 | Hz at `i=0` |
| Front (x+) | L1323–1356 | Hy/Hz at `i=NX` |
| Left (y−) | L1108–1145 | Hx/Hz at `j=0` |
| Right (y+) | L1151–1188 | Hx/Hz at `j=NY` |
| Down (z−) | L1194–1231 | Hy/Hx at `k=0` |
| Up (z+) | L1237–1274 | Hy/Hx at `k=NZ` |
| Past-field store | L1374+ | end of `applyMurH()` |

`InitMURBorders` (L57–138) → `FDTD_Solver::initMurBorders()` sets thin `MurZone` pads per face. `calc_murconstants` (L294–341) → `calcMurConstants()` sets per-face `*Cab1` from `(1−CNUM)/(1+CNUM)`, `CNUM = dx/(dt·c)` (vacuum).

### Fortran golden probe generation (one-time)

1. Build Fortran executable: `cmake --build build -j --target semba-fdtd`.
2. Add probes to `pulse-1d-x.fdtd.json` at the center cell (e.g. mesh coordinate `[3,3,3]`, direction `x`).
3. Seed the same initial Ex pulse as the C++ hook (Fortran test driver or patched init).
4. Run `./build/bin/semba-fdtd -i testData/cases/mur/pulse-1d-x.fdtd.json` for N steps.
5. Copy probe output to `testData/cases/mur/pulse-1d-x_probe_Ex_20_1_1.dat` (and steps 40/80/160 as needed).
6. Re-run `BordersMur.PulseMatchesFortranProbe` — it auto-enables when the file exists.
