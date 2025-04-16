![{63B38107-5D31-417E-8F6D-494075945788}](https://github.com/user-attachments/assets/32788197-5bb0-430c-a70a-5d4137bf4f55)[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

[![ubuntu](https://github.com/OpenSEMBA/fdtd/actions/workflows/ubuntu.yml/badge.svg?branch=main)](https://github.com/OpenSEMBA/fdtd/actions/workflows/ubuntu.yml)
[![windows](https://github.com/OpenSEMBA/fdtd/actions/workflows/windows.yml/badge.svg?branch=main)](https://github.com/OpenSEMBA/fdtd/actions/workflows/windows.yml)

# semba-fdtd

In a nutshell, semba-fdtd capabilities are

+ Cluster working capabilites through MPI.
+ Multiple threads per processor through OpenMP.
+ Closed/symmetric problems by means of PEC and PMC conditions.
+ Open problems by means of PML boundary conditions (CPML formulation) or by Mur ABCs.
+ Non-uniformly meshed domains by means of rectilinear (or graded) meshes.
+ Bulk lossless and lossy dielectrics.
+ Materials with frequency dependent relative permittivity and/or permeability, with an arbitrary number of complex-conjugate pole-residue pairs.
+ Bulk anisotropic lossless and lossy dielectrics.
+ Equivalent models of multilayered skin-depth materials.
+ Branched multiwires.
+ Multiconductor transmission lines networks embedded within 3D FDTD solvers.
+ Coupling with SPICE solvers (ngspice).
+ Junctions of wires of different radii.

  1. Junctions of multiwires.
  2. Wire bundles.
  3. Loaded with p.u.l resistance and inductance wires.
  4. Grounding through lumped elements.

+ Plane-wave illumination with arbitrary time variation.
+ Multiple planewaves illumination for reverberation chamber modeling.
+ Hertzian dipole sources.
+ Equivalent Huygens surfaces.
+ [Low frequency thin composites and lossy surfaces.](https://doi.org/10.1109/TMTT.2016.2637348)
+ [Thin slots.](https://doi.org/10.1109/TAP.2024.3484673)
+ Time, frequency and transfer function probes.
+ Near-to-far field transformation.

# Usage

Compilation and debugging instructions can be found [here](doc/development.md).

The main binary is `semba-fdtd` which uses [the SMBJSON format](doc/smbjson.md) as input files.
It can be run with

```shell
  semba-fdtd -i CASE_NAME.fdtd.json
```

Tests must be run from the root folder. `python` wrapper test assumes that `semba-fdtd` has been compiled successfully and is located in folder `build/bin/`. For intel compilation it also assumes that the intel runtime libraries are accessible.

# License

This code is licensed under the terms of the [MIT License](LICENSE). All rights reserved by the University of Granada (Spain)

## Funding

- Spanish Ministry of Science and Innovation (MICIN/AEI) (Grant Number: PID2022-137495OB-C31)
- European Union, HECATE project. (HE-HORIZON-JU-Clean-Aviation-2022-01)
- iSense Project. In-Situ Monitoring of Electromagnetic Interference. (HE-HORIZON-MSCA-2023-DN-01)
