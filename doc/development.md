# Compilation and debugging

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
  cmake -S . -B build -DHDF5_BUILD_FORTRAN=YES -DHDF5_ENABLE_Z_LIB_SUPPORT=NO --fresh
  cmake --build build -j 
  cmake --install build --prefix ~/hdf5-installed 
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
In the .gitmodules file, the submodules use the SSH remote URL by default. If not using a SSH-key in the computer where the following process will be performed, the remote addresses for each submodule must be individually changed to their HTTPS alternative.

This software requires [Windows BaseKit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html) and [Windows HPCKit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html). Install these packages with all their features selected.

Additionally, if not done already, install [CMake](https://cmake.org/download/) and [Ninja](https://github.com/ninja-build/ninja), follow their respective installation steps.

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

### Debugging with Visual Studio

1. **Install necessary tools**:
   - Install **Intel Base Kit** and **Intel HPC Kit**.
   - Install Visual Studio 2019. This must be done **after** the intel Intel compilers. Visual Studio 2022 is not supported (as of Oct. 2024).
   - Ensure **CMake** is installed. It's recommended to use **CMake GUI** to simplify configuration on Windows.

2. **Generate the Project with CMake GUI**:
   - In **CMake GUI**, select the option to create a project for **Visual Studio 2019** from the CMake files.
   - Set `-T fortran=ifx` to use intelLLVM compiler.
   - Specify the output folder where the `.sln` file for the project will be generated.

3. **Open the Project in Visual Studio**:
   - Open the generated `.sln` file in Visual Studio 2019.
   - You will see multiple projects in the solution, one for each library and executable.

4. **Compilation Configuration**:
   - You may encounter specific errors in the **fdtd-tests** projects during compilation, particularly related to the `/Mtd` or `/MD` options.
   - To resolve this issue:
     - Right-click on the problematic project (`fdtd_tests`) and select **Properties**.
     - Go to **Code Generation** and change the corresponding option to `/Mtd`.

5. **Set the Main Project**:
   - Select the `semba-fdtd` project and set it as the main project.

6. **Debugger Configuration**:
   - In the project properties, go to the **Debugger** section.
   - Configure the following parameters:
     - **Working Directory**: Set the working directory based on the case you want to debug.
     - **Command Arguments**: Enter `-i filename` (where `filename` is the required input file).

7. **Debug the Project**:
   - Launch the project by pressing **F5** to begin debugging.

Following these steps, the project should be set up and ready to compile and debug in Visual Studio 2019 with Intel tools.
