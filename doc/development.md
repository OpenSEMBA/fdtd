# Compilation and debugging

## Contents

- [Prebuilt binary releases](#running-from-prebuilt-binary-releases)
- [GNU/Linux compilation](#gnulinux-compilation)
  - [Compilation options](#compilation-options)
  - [HDF5 libraries](#hdf5-libraries)
  - [MTLN and ngspice](#mtln-and-ngspice)
  - [MPI](#mpi)
- [Windows (intelLLVM) compilation](#windows-intelllvm-compilation)
  - [Prerequisites](#prerequisites)
  - [Compilation process](#compilation-process)
  - [Visual Studio debugging](#debugging-with-visual-studio)
- [WSL2 and Visual Studio Code setup](#wsl2--visual-studio-code--gfortran-setup-guide)
- [Debugging the project](#debugging-the-project)
  - [MPI debugging](#debugging-with-mpi)
  - [Troubleshooting](#troubleshooting)

## Running from prebuilt binary releases

Prebuilt binares are available at [releases](https://github.com/OpenSEMBA/fdtd/releases).

In windows, you need to install [intel oneapi runtime libraries](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html).

## GNU/Linux Compilation

The repository has dependencies available as submodules. Before running CMake, initialise them from the repository root:

```shell
git submodule update --init --recursive
```

CMake presets use separate build directories, such as `build-rls/` and
`build-dbg-mpi/`, so their configurations can coexist.
Run CMake with `--fresh` when changing the options of an existing build
directory so its previous cache is discarded.

If you use intel oneapi compiler, make sure to run

```shell
  source /opt/intel/oneapi/compiler/latest/env/vars.sh
  export FC=ifx
  export CC=icx
  export CXX=icpx
```

### Compilation options

#### HDF5 Libraries

GNU builds use the system HDF5 installation by default.
Intel builds use the bundled serial HDF5 installation unless a different
installation is selected with `HDF5_ROOT` or `HDF5_DIR`.

You can compile HDF5 for your platform by downloading the latest sources from
the [HDF5 website](https://www.hdfgroup.org/downloads/hdf5/source-code/).
Extract the archive, then build and install a serial version with:

```shell
cmake -S . -B build \
  -DHDF5_BUILD_FORTRAN=ON \
  -DHDF5_ENABLE_Z_LIB_SUPPORT=NO \
  --fresh
cmake --build build -j
cmake --install build --prefix ~/hdf5-installed
```

A specific HDF5 installation can be selected with
`-DHDF5_ROOT=<path-to-library>`:

```shell
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DHDF5_ROOT=~/hdf5-installed \
  -DHDF5_USE_STATIC_LIBRARIES=TRUE \
  --fresh
cmake --build build -j
```

Parallel movie output requires HDF5 built with both MPI and Fortran support.
Configure HDF5 with MPI compiler wrappers and `HDF5_ENABLE_PARALLEL=ON`:

```shell
CC=mpicc FC=mpifort cmake -S . -B build-hdf5-parallel \
  -DHDF5_BUILD_FORTRAN=ON \
  -DHDF5_ENABLE_PARALLEL=ON \
  -DHDF5_ENABLE_Z_LIB_SUPPORT=NO
cmake --build build-hdf5-parallel -j
cmake --install build-hdf5-parallel --prefix ~/hdf5-parallel
```

MPI builds automatically prefer parallel HDF5 when it is available.
Select the parallel HDF5 installation while enabling the normal MPI option:

```shell
cmake -S . -B build-parallel \
  -DSEMBA_FDTD_ENABLE_MPI=ON \
  -DSEMBA_FDTD_ENABLE_HDF=ON \
  -DSEMBA_FDTD_ENABLE_OUTPUT_MODULE=ON \
  -DHDF5_ROOT=~/hdf5-parallel
```

If HDF5 reports parallel support, configuration verifies that its Fortran MPIO
interfaces compile and link with the selected MPI library.
MPI builds using serial HDF5 remain valid, but the parallel HDF5 movie backend
is not available in those builds.

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
This repository has dependencies available as submodules. Initialise them from the root folder before running CMake:

```shell
git submodule update --init --recursive
```

The default submodule URLs use HTTPS.

This software requires [Windows BaseKit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html) and [Windows HPCKit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html). Install these packages with all their features selected.

Additionally, if not done already, install [CMake](https://cmake.org/download/) and [Ninja](https://github.com/ninja-build/ninja), follow their respective installation steps.

### Compilation process

Open a command prompt with OneAPI variables initialised, to do this open a new command prompt and type:

```shell
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
```

This will load the OneAPI environment for x64.

Navigate to the fdtd root folder, choose between "Debug"/"Release" for `-DCMAKE_BUILD_TYPE`, and "ON"/"OFF" for `-DSEMBA_FDTD_ENABLE_MPI`, for example, a Release version with MPI Support would be:

```shell
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release -DSEMBA_FDTD_ENABLE_MPI=ON --fresh
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
   - Install Visual Studio 2022. This must be done **after** the intel Intel compilers.
   - Ensure **CMake** is installed. 

3. **Open a terminal with the Intel One API variables loaded**:
   - Launch Visual Studio 2022 with
     ```
     "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\devenv.exe"  
     ```
   - Open the opensemba/fdtd cloned repo as a folder.


4. **Ensure cmake has nofpp option**:
   - The following command 
     ```
     add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:/nofpp>)
     ```
     must be present in CMakeLists.txt due to an issue https://gitlab.kitware.com/cmake/cmake/-/issues/21816

6. **Debugger Configuration**:
   - Select your current build as start up item
   - In Debug options, open `Debug and launch settings` for your start up item
   - Configure the following parameters:
     - **Working Directory**: Set the working directory based on the case you want to debug.
     - **Command Arguments**: Enter `-i filename` (where `filename` is the required input file).
     You should end up with something like this.
     ```
      {
        "configurations": [
                {
                    "type": "default",
                    "project": "CMakeLists.txt",
                    "projectTarget": "semba-fdtd.exe (bin\\semba-fdtd.exe)",
                    "name": "semba-fdtd.exe (bin\\semba-fdtd.exe)",
                    "currentDir": "<root-folder>\\tmp_cases\\coated_antenna",
                    "args": [
                        "-i",
                        "coated_antenna.fdtd.json",
                        "-mtlnwires"
                        ]
                }
        ],
        "defaults": {},
        "version": "0.2.1"
      }
     ```

7. **Debug the Project**:
   - Launch the project by pressing **F5** to begin debugging.

## WSL2 + Visual Studio Code + GFortran Setup Guide

### Prerequisites

- **Windows 10** (version 1903 or higher) or **Windows 11**
- **WSL2** installed
- **Visual Studio Code (VSCode)** installed on Windows

### Step 1: Install WSL2

#### 1.1 Install WSL (Windows Subsystem for Linux)

1. Open **PowerShell** as Administrator and run the following command to enable the required features:

```bash
  wsl --install
  wsl --set-default-version 2
```

2. After installation is complete, you can choose a Linux distribution from the Microsoft Store (e.g., Ubuntu). for example if you want to install Ubuntu:

```bash
  wsl --install -d Ubuntu-24.04
```


3. Launch the installed distribution from the Start menu, and it will complete the installation by setting up a user.

### Step 2: Install Visual Studio Code

#### 2.1 Download Visual Studio Code
Go to the official Visual Studio Code website and download the installer:
https://code.visualstudio.com/

Once downloaded, double-click the installer and follow the on-screen instructions to complete the installation.

#### 2.2 Launch Visual Studio Code

After the installation is complete, you can open Visual Studio Code either by:
- Searching for "Visual Studio Code" in the Start Menu
- Or by launching the VSCode application from the shortcut created during installation.

#### 2.3 Installing VSCode Extensions

To work with these project is mandatory to install the next extensions

- Modern Fortran (Fortran development)
- Python (Python development)
- C/C++ (C/C++ development)
- CMake (CMake integration)
- C++ TestMate (C++ testing)
- Remote - WSL

#### 2.4 Connect to WSL2

After installing the Remote - WSL extension, follow these steps to open VSCode in the WSL2 environment:
1. Press Ctrl + Shift + P to open the Command Palette.
2. Type `Remote-WSL: New Window` and press Enter.
3. VSCode will open a new window and connect to your default WSL2 distribution.

Now you're set up to work with your WSL2 environment directly from VSCode

Now we are ready to clone the repo

#### Step 1: Clone the Repository
```bash
git clone <repository_url>
cd <repository_name>
```

This project has submodule dependencies remember to initiate an update the

```bash
git submodule update --init --recursive
```

#### Step 2: Install Python Requirements

If the project has Python dependencies listed in a requirements.txt file, you can install them using the following command:

1. Make sure you have Python and pip installed on your system. It is recomended to use a venv.
2. Run the following command to install the required dependencies:

```bash
python3 -m pip install -r requirements.txt
```

#### Step 3: Install hdf5 dependencies
It is needed to install manually the hdf5 dependencies. On the terminal type the next command:
```
sudo apt install libhdf5-dev libopenmpi-dev
```
#### Step 4: Build the Project with CMake

Now, you need to build the project using CMake.

1. Open the Command Palette in Visual Studio Code by pressing Ctrl + Shift + P.
2. Type CMake: Build and select it from the list of commands.

   This will trigger the build process using the current CMake configuration.


#### Step 5: Change CMake Settings (Optional)

If you want to modify the CMake settings (such as build options or configurations), follow these steps:

1. Open the Command Palette in Visual Studio Code by pressing Ctrl + Shift + P.
2. Type CMake: Configure and select it from the list of commands.

   This will open the CMake configuration interface where you can modify build settings.

3. After making any changes to the settings, rebuild the project by running the CMake: Clean Rebuild command again from the Command Palette.


#### Step 6: Configure and Run Python Tests

To run Python unit tests, you need to configure the testing framework.

1. Open the Command Palette in Visual Studio Code by pressing Ctrl + Shift + P.
2. Type Python: Configure Tests and select it from the list of commands.
3. When prompted, select unittest as the testing framework.
4. Select the folder containing your test files. Typically, this folder is named tests or something similar. Once selected, Visual Studio Code will automatically configure and discover the tests.

To run the tests:

1. Open the Command Palette again by pressing Ctrl + Shift + P.
2. Type Python: Run All Tests to run the unit tests in your project.

## Debugging the project

For a correct debugging experience configuring a launch.json file is needed. This file usually is created by vscode automatically. In case it does not exist. You can create your own on .vscode folder.

An example of launch.json filke is given. This will use a file as argument when calling to semba-fdtd.
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Fortran Launch (GDB)",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceRoot}/build-dbg/bin/semba-fdtd",
            "miDebuggerPath": "gdb",
            "args": ["-i", "shieldingEffectiveness.fdtd.json"],
            "stopAtEntry": false,
            "cwd": "${workspaceRoot}/tmp_cases/sgbcShieldingEffectiveness/"
        }
    ]
}
```

Now you are ready to work with the project.

### Debugging with MPI

#### Overview

GDB controls one process per debug session.
To debug an MPI job, each MPI rank is therefore started under its own
`gdbserver`, and VS Code creates one `cppdbg` session for each rank.

The checked-in configuration supports a two-rank solver job:

```text
VS Code: MPI: debug all ranks (2 ranks)
  |-- GDB session: rank 0 -> localhost:20000 -> gdbserver -> MPI rank 0
  `-- GDB session: rank 1 -> localhost:20001 -> gdbserver -> MPI rank 1
```

Both sessions are shown separately in the VS Code **Call Stack** panel.
Breakpoints are sent to both sessions, although rank-specific control flow can
mean that only one rank reaches a particular breakpoint.

#### Configuration files

| File | Responsibility |
|---|---|
| `.vscode/launch.dev.json` | Version-controlled template for debug configurations. |
| `.vscode/launch.json` | Active local configuration; ignored by Git. |
| `.vscode/settings.json` | Local input file, working directory, and test filter values. |
| `.vscode/tasks.json` | Build tasks and the MPI-safe `fdtd_tests` preparation task. |
| `scripts/debug-mpi-gdbserver.sh` | Starts MPI ranks under `gdbserver`, waits for ports, and cleans stale jobs. |

Copy the version-controlled files template when initially configuring the
workspace, or when the template changes:

```shell
cp .vscode/launch.dev.json .vscode/launch.json
```

Add the following project-specific values to the local
`.vscode/settings.json` file:

```json
{
  "semba-fdtd.debug.inputFile": "pw-in-box.fdtd.json",
  "semba-fdtd.debug.inputCwd": "testData/cases/planewave",
  "semba-fdtd.debug.mpiGtestFilter": "conformal.geometry_coord_position"
}
```

`inputFile` is relative to `inputCwd`.
Set `inputCwd` to the directory containing the JSON input and all files that
the JSON references with relative paths, such as excitation files.

#### Build requirements

Debug launches do not configure or rebuild the project automatically.
Build an MPI-enabled Debug executable before starting VS Code debugging:

```shell
cmake --fresh --preset dbg-mpi
cmake --build --preset dbg-mpi -j
```

The `dbg-mpi` preset builds in `build-dbg-mpi/`.
Its configuration can coexist with other presets, so a Release or non-MPI
preset build does not replace the MPI Debug executable.
If the options of `dbg-mpi` itself change, rerun the commands above to refresh
that build directory.

#### Starting all ranks

1. Open the VS Code **Run and Debug** view.
2. Select `MPI: debug all ranks (2 ranks)`.
3. Press F5.
4. Wait for both `MPI all ranks: rank 0` and
   `MPI all ranks: rank 1` to appear in **Call Stack**.
5. Continue each session once after the initial entry stop.

The two ranks must both be allowed to continue.
If one remains stopped before `MPI_Init` or another collective operation, the
other rank can appear blocked while it waits for that rank.

`MPI: debug solver rank 0 (2 ranks)` is a simpler alternative.
It runs a two-rank MPI job but attaches GDB only to rank 0;
rank 1 runs normally.

#### Startup sequence

The all-rank configuration is a VS Code compound containing two hidden launch
configurations.
No `preLaunchTask` or problem matcher is used for the solver.
Instead, the C/C++ extension directly owns the helper processes through
`debugServerPath` and waits for readiness through `serverStarted`.

The startup sequence is:

1. The rank 0 launch configuration invokes
   `scripts/debug-mpi-gdbserver.sh` with `--foreground-all 2`.
2. The script validates `--workdir`, changes to that directory, and executes
   one `mpirun -np 2` job.
3. Each MPI process calculates its debugger port as
   `20000 + OMPI_COMM_WORLD_RANK` and then executes `gdbserver`.
4. The rank 0 adapter waits for `Listening on port 20000` and connects its GDB.
5. The rank 1 launch configuration invokes the same script with
   `--wait-for-port 20001` instead of starting a second MPI job.
6. The waiter checks `/proc/net/tcp` and `/proc/net/tcp6` without opening a
   debugger connection.
7. When port 20001 is listening, the rank 1 adapter connects its own GDB.
8. The compound's `stopAll` option stops both sessions when either session is
   terminated.

It is important that rank 1 only waits for its port.
Starting `mpirun` from both hidden configurations would create two unrelated
MPI jobs rather than two debugger views of the same job.

#### Port and rank mapping

| MPI rank | `gdbserver` address | VS Code session |
|---:|---|---|
| 0 | `localhost:20000` | `MPI all ranks: rank 0` |
| 1 | `localhost:20001` | `MPI all ranks: rank 1` |

The shell script supports more ranks, but the checked-in compound explicitly
defines two GDB sessions.
Supporting additional ranks requires another hidden launch configuration and
port waiter for each extra rank.

#### Working directory

The target process is launched by `gdbserver`, not directly by `cppdbg`.
Consequently, the `cwd` property alone does not reliably set the inferior's
working directory.

The launch configuration passes the directory explicitly:

```text
--workdir ${workspaceFolder}/${config:semba-fdtd.debug.inputCwd}
```

The script verifies that the directory exists and is writable, then changes to
it before starting `mpirun`.
This is required for relative JSON resources, output files, and solver control
files such as `running`, `pause`, `relaunch`, and `forcestop`.

#### Script modes

The helper script has the following modes:

| Mode | Purpose |
|---|---|
| `--foreground-all <ranks> <program> ...` | Run every MPI rank under a separate `gdbserver`. |
| `--foreground-debug-rank <rank> <ranks> <program> ...` | Debug one rank and run the remaining ranks normally. |
| `--wait-for-port <port>` | Wait for another launch configuration's `gdbserver`. |
| `--debug-rank <rank> <ranks> <program> ...` | Detached preparation mode used by task-based workflows. |
| `--stop` | Stop a detached MPI debug job recorded by the script. |

The `--foreground-*` modes are preferred for solver debugging because
`OpenDebugAD7` owns their lifetime directly.
This avoids races in which a background task exits before GDB connects.

#### Debugging MPI unit tests

The full `fdtd_tests` suite should normally be debugged as one process, even
when linked against an MPI-enabled build.
Several tests write fixed file names and are not safe to execute concurrently
on every rank.

The `MPI: debug fdtd_tests (2 ranks)` compound is intended only for an
MPI-safe filtered test.
Set `semba-fdtd.debug.mpiGtestFilter` in `.vscode/settings.json` before using
that compound.

#### Troubleshooting

**GDB connection timeout**

Confirm that the selected configuration is
`MPI: debug all ranks (2 ranks)` and reload the VS Code window after changing
`launch.json`.
The solver configuration must use `debugServerPath`, `serverStarted`, and the
foreground script modes; it must not depend on a background `preLaunchTask`.

Check for stale MPI or `gdbserver` processes before retrying:

```shell
pgrep -af 'gdbserver|prterun|mpirun'
```

**Cannot create `running` or another relative file**

Verify `semba-fdtd.debug.inputCwd` and confirm that the directory is writable.
The debug output prints `MPI working directory: ...` before `mpirun` starts.

**A breakpoint is not reached**

Confirm that the correct rank executes that code path and that the breakpoint
was installed before the one-time initialization code ran.
Also confirm that the active executable is an MPI-enabled Debug build.
A valid source breakpoint cannot force execution through a false runtime
condition.

**Both sessions connect but the program does not advance**

Select each rank in **Call Stack** and continue it.
One stopped rank can hold the other rank inside an MPI collective operation.

**Warnings about unavailable system-library debug information**

Messages about missing separate debug information for MPI or system libraries
are non-fatal when debugging project sources.
Install the corresponding system debug packages only when stepping inside those
libraries is required.

#### Manual attach fallback

The native compound is preferred, but GDB can also attach manually to an
already running process.
Start the MPI job in a terminal:

```shell
mpirun -np 2 build-dbg-mpi/bin/semba-fdtd -i input_file.fdtd.json
```

Then use the `Attach to process` configuration and select one `semba-fdtd`
process.
This method provides one attached rank per debug session and does not perform
the automatic port coordination described above.

If Linux blocks manual attachment because of `ptrace_scope`, temporarily relax
the restriction only on a trusted development machine:

```shell
sudo sysctl kernel.yama.ptrace_scope=0
```

Restore the normal restriction after debugging:

```shell
sudo sysctl kernel.yama.ptrace_scope=1
```

See the [MIEngine troubleshooting guide][miengine-troubleshooting]
for more information.

[miengine-troubleshooting]: https://github.com/Microsoft/MIEngine/wiki/Troubleshoot-attaching-to-processes-using-GDB


