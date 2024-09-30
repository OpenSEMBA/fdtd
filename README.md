[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

[![ubuntu-gnu](https://github.com/OpenSEMBA/fdtd/actions/workflows/ubuntu-gnu.yml/badge.svg)](https://github.com/OpenSEMBA/fdtd/actions/workflows/ubuntu-gnu.yml)
[![ubuntu-intelLLVM](https://github.com/OpenSEMBA/fdtd/actions/workflows/ubuntu-intelLLVM.yml/badge.svg)](https://github.com/OpenSEMBA/fdtd/actions/workflows/ubuntu-intelLLVM.yml)
[![windows-intelLLVM](https://github.com/OpenSEMBA/fdtd/actions/workflows/windows-intelLLVM.yml/badge.svg)](https://github.com/OpenSEMBA/fdtd/actions/workflows/windows-intelLLVM.yml)

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
+ Low frequency thin composites and lossy surfaces.
+ Thin slots.
+ Time, frequency and transfer function probes.
+ Near-to-far field transformation.

# Usage

## Running from prebuilt binary releases

Prebuilt binares are available at [releases](https://github.com/OpenSEMBA/fdtd/releases).

In windows, you need to install [intel oneapi runtime libraries](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html).

## GNU/Linux Compilation

It is important to point out the repository has dependencies which are available as submodules. It is necessary to run `git submodule init` and `git submodule update` from the root folder before running any `cmake` or `build` commands.

If you use intel oneapi compiler, make sure to run

```shell
  source /opt/intel/oneapi/compiler/latest/env/vars.sh
  export FC=ifx
  export CC=icx
  export CXX=icpx
```

### Compilation options

#### HDF5 Libraries

HDF5 precompiled libraries for ubuntu are used by default (see [precompiled libraries cmake script](set_precompiled_libraries.cmake)).  

You can compile HDF5 for your specific platform downloading the latest sources from this [link](https://www.hdfgroup.org/downloads/hdf5/source-code/).
Extract to a folder and build and install with the following commands

```shell
  cmake -S . -B build-ifx --fresh -DHDF5_BUILD_FORTRAN=YES
  cmake --build build-ifx -j 
  cmake --install build-ifx --prefix ~/hdf5-installed 
```

A specific HDF5 library can be set with the option `-DHDF5_ROOT=<path-to-library>`, e.g.

```shell
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=~/hdf5-installed -DHDF5_USE_STATIC_LIBRARIES=TRUE --fresh 
  cmake --build build -j
```

#### MTLN and ngspice

MTLN depends on `lapack` and `ngspice`. Precompiled versions are included for windows (intelLLVM) and ubuntu (intelLLVM and GNU).
For other platform/compilers these will need to be compiled.

##### Compiling ngspice

In linux, when using some of the provided scripts you may find problems with carriage returns. These can be fixed with:

```shell
sed -i -e 's/\r$//' compile_linux.sh
sed -i -e 's/\r$//' autogen.sh
find . -name \*.m4|xargs dos2unix\nfind . -name \*.ac|xargs dos2unix\nfind . -name \*.am|xargs dos2unix
```

the `ngspice` static library can be compiled doing the following:

1. Edit `configure.ac`, to `AC_SUBST([STATIC], [-static])`
2. Edit `compile_linux_shared.sh`, to 
```
libngspice_la_CFLAGS = -static

libngspice_la_LDFLAGS = -static -version-info @LIB_VERSION@
```

#### MPI

If you use intel oneapi, make sure to load the mpi environment variables:

```shell
  source /opt/intel/oneapi/mpi/latest/env/vars.sh
```

## Windows (intelLLVM) Compilation

Clone this repository:

```shell
  git clone https://github.com/OpenSEMBA/fdtd.git
```

or, if using SSH keys:

```shell
  git clone git@github.com:OpenSEMBA/fdtd.git
```

navigate to the `/fdtd/` folder that has been created, this folder will be referred to as `root` for any future purposes.

### Prerequisites

This compilation process will use the already available precompiled libraries included with the project, thus it's not required to build them manually.
This repository has dependencies that are available as submodules. It is necessary to run `git submodule init` and `git submodule update` from the root folder before running any `cmake` or `build` commands.

Additionally, this software requires [Windows BaseKit](https://registrationcenter-download.intel.com/akdlm/IRC_NAS/62641e01-1e8d-4ace-91d6-ae03f7f8a71f/w_BaseKit_p_2024.0.0.49563_offline.exe) and [Windows HPCKit](https://registrationcenter-download.intel.com/akdlm/IRC_NAS/5b36181e-4974-4733-91c7-0c10c54900a5/w_HPCKit_p_2024.0.0.49588_offline.exe). Install these packages with all their features selected.

If not done already, install [CMake](https://cmake.org/download/) and [Ninja](https://github.com/seanmiddleditch/gha-setup-ninja), follow their respective installation steps.

### Compilation process

Open a command prompt with OneAPI variables initialised, to do this open a new command prompt and type:

```shell
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
```

This will load the OneAPI environment for x64.

Navigate to the fdtd root folder, choose between "Debug"/"Release" for `-DCMAKE_BUILD_TYPE`, and "Yes"/"No" for `-DSEMBA_FDTD_ENABLE_MPI`, for example, a Release version with MPI Support would be:

```shell
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release -DSEMBA_FDTD_ENABLE_MPI=Yes
```

Then,

```shell
cmake --build build -j
```

We should now find the compiled executables in `\build\bin\`.

### Usage

In order to use semba-fdtd, the executable must have access to the dynamic libraries it has dependencies on. Either move the libraries to the same folder as the executable, or run the executable through a console with the OneAPI environment loaded:

```shell
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
```

Once the environment is loaded, follow the steps in the next section.

## Running cases

`semba-fdtd` uses the following [format](doc/smbjson.md) as input files.
These can be run with

```shell
  semba-fdtd -i CASE_NAME.fdtd.json
```

## Testing

Tests must be run from the root folder. `python` wrapper test assumes that `semba-fdtd` has been compiled successfully and is located in folder `build/bin/`. For intel compilation it also assumes that the intel runtime libraries are accessible.

# License

This code is licensed under the terms of the [MIT License](LICENSE). All rights reserved by the University of Granada (Spain)

# References

+ Miguel Ruiz Cabello, Maksims Abalenkovs, Luis Diaz Angulo, Clemente Cobos Sanchez, Franco Moglie, Salvador Gonzalez Garcia, *Performance of parallel FDTD method for shared- and distributed-memory architectures: Application to bioelectromagnetics*. PLOS ONE. 2020. https://doi.org/10.1371/journal.pone.0238115
+ Luis Diaz Angulo, Miguel Ruiz Cabello, Jesus Alvarez, Amelia Rubio Bretones, Salvador Gonzalez Garcia, *From Microscopic to Macroscopic Description of Composite Thin Panels: A Road Map for Their Simulation in Time Domain*. IEEE Transactions on Microwave Theory and Techniques. 2018. https://doi.org/10.1109/TMTT.2017.2786263.
+ Miguel Ruiz Cabello, Luis Diaz Angulo, Jesus Alvarez, Ian Flintoft, Samuel Bourke, John Dawson, *A Hybrid Crank–Nicolson FDTD Subgridding Boundary Condition for Lossy Thin-Layer Modeling*. IEEE Transactions on Microwave Theory and Techniques. 2017. https://doi.org/10.1109/TMTT.2016.2637348.
+ Miguel Ruiz Cabello, Luis Diaz Angulo, Amelia Rubio Bretones, Rafael Gomez Martin, Salvador Gonzalez Garcia and Jesus Alvarez, *A novel subgriding scheme for arbitrarily dispersive thin-layer modeling*, 2017 IEEE MTT-S International Conference on Numerical Electromagnetic and Multiphysics Modeling and Optimization for RF, Microwave, and Terahertz Applications (NEMO), Seville, Spain, 2017.
https://doi.org/10.1109/NEMO.2017.7964255.
+ Guadalupe Gutierrez Gutierrez, Daniel Mateos Romero, Miguel Ruiz Cabello, Enrique Pascual-Gil, Luis Diaz Angulo, David Garcia Gomez, Salvador Gonzalez Garcia, 
*On the Design of Aircraft Electrical Structure Networks*, 
IEEE Transactions on Electromagnetic Compatibility. 2016. https://doi.org/10.1109/TEMC.2016.2514379.
