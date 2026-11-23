# Compilation and debugging

## Contents

- [Prebuilt binary releases](#running-from-prebuilt-binary-releases)
- [GNU/Linux compilation](#gnulinux-compilation)
  - [Compilation options](#compilation-options)
  - [HDF5 libraries](#hdf5-libraries)
  - [MTLN and ngspice](#mtln-and-ngspice)
  - [MPI](#mpi)
- [Windows (Intel oneAPI) compilation](#windows-intel-oneapi-compilation)
  - [Prerequisites](#prerequisites)
  - [Build environment](#build-environment)
  - [Compilation process](#compilation-process)
  - [MPI on Windows](#mpi-1)
  - [Usage](#usage)
  - [Visual Studio Code (native Windows)](#visual-studio-code-native-windows)
  - [Visual Studio debugging](#debugging-with-visual-studio)
  - [Windows troubleshooting](#windows-troubleshooting)
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
Select the parallel HDF5 installation while enabling MPI:

```shell
cmake -S . -B build-parallel \
  -DSEMBA_FDTD_ENABLE_MPI=ON \
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

## Windows (Intel oneAPI) Compilation

Native Windows builds use the Intel oneAPI compilers (`ifx` for Fortran, `icx` for C and C++)
together with the MSVC linker and Windows SDK, CMake and Ninja.
HDF5 and LAPACK are already provided in `precompiled_libraries/` and are selected automatically, so they do not need to be built.

The steps below were validated on Windows 11 with the versions listed in the prerequisites.
The Windows CI workflow (`.github/workflows/windows.yml`) uses the same Intel compiler version.

### Prerequisites

#### Minimal installation

| Component | Validated version | Why it is needed |
|---|---|---|
| [Git for Windows](https://git-scm.com/download/win) | any | Clone the repository and its submodules. |
| [Visual Studio 2022 Build Tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022), workload *Desktop development with C++* (`Microsoft.VisualStudio.Workload.VCTools`) with its recommended components | 17.14 (MSVC 14.44) | `ifx` uses the MSVC `link.exe` and the Windows SDK. The recommended component *C++ CMake tools for Windows* also provides CMake 3.31 and Ninja. A full Visual Studio 2022 installation with the same workload works too. |
| Intel Fortran compiler, component `intel.oneapi.win.ifort-compiler` | 2025.2.1 | Provides `ifx`. |
| Intel C/C++ compiler, component `intel.oneapi.win.cpp-dpcpp-common` | 2025.2.1 | Provides `icx`. ngspice needs a C99 compiler with complex number support, which MSVC `cl` does not provide. |
| Intel MPI Library, component `intel.oneapi.win.mpi.devel` | 2021.16.1 | **Only for MPI builds.** Installing the compilers brings only the MPI runtime, not the headers, Fortran modules and libraries. |
| Python 3 | 3.x | **Optional.** Python tests (`requirements.txt`) and the `fortls` language server for VS Code. |

> **CMake version:** CMake 4.x cannot detect the `ifx` Fortran preprocessor on Windows.
> Use CMake 3.x. The one bundled with Visual Studio 2022 17.14 is 3.31.6.
> If you install CMake separately, use 3.31.x.

Both installers need administrator rights.

Install the Visual Studio Build Tools, for example with `winget`:

```powershell
winget install --id Microsoft.VisualStudio.2022.BuildTools -e --override "--quiet --wait --norestart --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
```

Install only the required Intel components from the offline HPC Toolkit installer.
This is the same installer used by the CI through `fortran-lang/setup-fortran`. The download is about 2.8 GB:

```powershell
curl.exe -L -C - --retry 5 -o hpckit.exe https://registrationcenter-download.intel.com/akdlm/IRC_NAS/e63ac2b4-8a9a-4768-979a-399a8b6299de/intel-oneapi-hpc-toolkit-2025.2.1.46_offline.exe
.\hpckit.exe -s -x -f hpckit
.\hpckit\bootstrapper.exe -s --action install --eula=accept --components=intel.oneapi.win.ifort-compiler:intel.oneapi.win.cpp-dpcpp-common:intel.oneapi.win.mpi.devel -p=NEED_VS2017_INTEGRATION=0 -p=NEED_VS2019_INTEGRATION=0 -p=NEED_VS2022_INTEGRATION=0
```

- The installer runs silently and shows no window. It is done when the command returns.
- Leave out `intel.oneapi.win.mpi.devel` if you do not need MPI.
- To add components later, use `--action modify` and list **all** the components, including the ones already installed.
  Running `--action install` again fails with *"It is already installed"*.
- Use `-p=NEED_VS2022_INTEGRATION=1` if you also want to debug from the Visual Studio IDE.

#### Full installation (alternative)

Instead of the individual components, you can install the
[Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html) and the
[Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)
with all their features selected. This includes everything listed above plus many unused libraries, so it needs much more disk space.
You still need the Visual Studio 2022 C++ workload and a CMake 3.x.

#### Repository

```shell
git clone --recurse-submodules https://github.com/OpenSEMBA/fdtd.git
```

or, in an existing clone, `git submodule update --init --recursive`.

### Build environment

All commands must run in a `cmd` console with the MSVC and oneAPI environments loaded.
`scripts/windows/env.bat` does this:

```bat
call scripts\windows\env.bat
```

It performs these steps:

1. Clears `NoDefaultCurrentDirectoryInExePath`.
   If set, `setvars.bat` fails because it calls `vars.bat` through a relative path.
2. Finds Visual Studio with `vswhere.exe` and calls `vcvars64.bat`.
   This provides `link.exe`, the Windows SDK, CMake and Ninja.
   `setvars.bat` does not load MSVC by itself.
3. Calls the oneAPI `setvars.bat intel64`. This provides `ifx`, `icx` and, if installed, Intel MPI.
4. Sets `CC=icx`, `CXX=icx` and `FC=ifx`. Otherwise CMake would select MSVC `cl` for C.
5. Warns if CMake is 4.x, or if the first `link.exe` on `PATH` is not the MSVC linker.
   Git for Windows ships a `link.exe` in `usr\bin`.

Optional variables:
- `ONEAPI_ROOT`: oneAPI installation folder. Defaults to `C:\Program Files (x86)\Intel\oneAPI`.
- `SEMBA_TOOLS_PATH`: folders prepended to `PATH`, for example a standalone CMake 3.x and Ninja.

The equivalent manual steps are:

```bat
set NoDefaultCurrentDirectoryInExePath=
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
set CC=icx
set CXX=icx
set FC=ifx
```

### Compilation process

Use the CMake presets. `scripts/windows/build.bat` loads the environment, then configures and builds a preset:

```bat
scripts\windows\build.bat rls
```

The equivalent commands in a console with the environment loaded are:

```bat
cmake --preset rls
cmake --build --preset rls
```

- The executables are written to `build-<preset>\bin\`, for example `build-rls\bin\semba-fdtd.exe`.
- The presets that work on Windows are `rls`, `dbg`, `rls-nomtln`, `dbg-nomtln`, `rls-mpi` and `dbg-mpi`.
  The `intel-*` presets are Linux only: they use `/opt/intel` paths and the `mpiifx` wrappers.
- Add `--fresh` to the configure step when you change options in an existing build folder.

### MPI

MPI builds need the Intel MPI development component (see [prerequisites](#minimal-installation)).
Build with the `rls-mpi` or `dbg-mpi` preset:

```bat
scripts\windows\build.bat rls-mpi
```

- The bundled Windows HDF5 is serial, so the configuration reports `HDF5 parallel: OFF`.
  Everything works except the parallel HDF5 movie backend.
- Windows builds are compiled without OpenMP. MPI builds link the Intel OpenMP runtime (`libiomp5md`) only because
  `timestepping.F90` calls `omp_set_num_threads`.

### Usage

The executable needs the oneAPI runtime DLLs. Run it with `scripts/windows/run.bat` from the case folder,
or from a console with the environment loaded:

```bat
cd testData\cases\planewave
..\..\..\scripts\windows\run.bat -i pw-in-box.fdtd.json
```

`run.bat` uses `build-rls` by default. Set `SEMBA_BUILD_DIR` to use another build, for example `set SEMBA_BUILD_DIR=build-dbg`.

Run MPI builds with `mpiexec` through `scripts/windows/run-mpi.bat`, which uses `build-rls-mpi` by default:

```bat
..\..\..\scripts\windows\run-mpi.bat 2 -i pw-in-box.fdtd.json
```

- Local `mpiexec` runs need neither the `hydra_service` nor any credentials.
- An MPI executable also runs without `mpiexec`, as a single process.

The native tests run with CTest from a console with the environment loaded:

```bat
ctest --test-dir build-rls --output-on-failure
```

In MPI build folders, CTest also runs the MPI tests through `mpiexec`.

> Known issue: `sync_timeout_aborts_on_stalled_neighbour` times out with Intel MPI on Windows.
> The test stalls rank 1 with the Unix command `sleep 120`, which runs if Git's `usr\bin` is on `PATH`.
> Rank 0 then calls `MPI_Abort(SUBCOMM_MPI, ...)`, but the test never initialises `SUBCOMM_MPI`.
> Intel MPI rejects that communicator, so the stalled rank is not killed.

### Visual Studio Code (native Windows)

#### Extensions

The minimal set of extensions is:

| Extension | Why |
|---|---|
| Modern Fortran (`fortran-lang.linter-gfortran`) | Fortran syntax, hover, go to definition and symbols, through `fortls`. |
| CMake Tools (`ms-vscode.cmake-tools`) | Select a preset, build, and run CTest. |
| C/C++ (`ms-vscode.cpptools`) | Provides the `cppvsdbg` debugger, which reads the `.pdb` files produced by `ifx`. `gdb` does not work with these binaries. |

The other extensions in `.vscode/extensions.json` are optional on Windows.

Install `fortls` so that Modern Fortran can navigate the code:

```shell
python -m pip install --user fortls
```

If `fortls.exe` is not on `PATH`, set `fortran.fortls.path` to its location.
With `--user` it is installed in `%APPDATA%\Python\Python3XX\Scripts`.

#### Opening the project

VS Code, CMake Tools, the integrated terminal and the debugger inherit the environment of the process that starts VS Code.
Close every VS Code window, then open the project with:

```bat
scripts\windows\open-vscode.bat
```

> If another VS Code window is already open, the new one joins that process and does **not** get the environment.
> CMake then fails with *"The Fortran compiler identification is unknown"* or *"could not find 'link'"*.

Then select a configure preset, for example `dbg`, in the CMake Tools status bar, and build with **F7**.
The *CMake: configure/build* tasks in `.vscode/tasks.json` work as well.
The *Test:* and MPI debugging tasks in that file are Linux only.

#### Settings and launch configurations

`.vscode/settings.json` and `.vscode/launch.json` are not tracked, so create them locally.
A minimal `settings.json`:

```json
{
  "cmake.useCMakePresets": "always",
  "cmake.configureOnOpen": false,
  "fortran.fortls.disableAutoupdate": true,
  "fortran.linter.compiler": "Disabled",
  "fortran.formatting.formatter": "Disabled",
  "C_Cpp.default.configurationProvider": "ms-vscode.cmake-tools",
  "terminal.integrated.defaultProfile.windows": "Command Prompt"
}
```

- The linter is disabled because its default compiler is `gfortran`.
- `Command Prompt` makes the integrated terminal inherit the build environment.

A `launch.json` for the `dbg` and `dbg-mpi` presets:

```json
{
  "version": "0.2.0",
  "configurations": [
    {
      "name": "Win: debug semba-fdtd (dbg)",
      "type": "cppvsdbg",
      "request": "launch",
      "program": "${workspaceFolder}/build-dbg/bin/semba-fdtd.exe",
      "args": ["-i", "pw-in-box.fdtd.json"],
      "cwd": "${workspaceFolder}/testData/cases/planewave",
      "console": "integratedTerminal"
    },
    {
      "name": "Win: debug fdtd_tests (dbg)",
      "type": "cppvsdbg",
      "request": "launch",
      "program": "${workspaceFolder}/build-dbg/bin/fdtd_tests.exe",
      "args": ["--gtest_filter=*"],
      "cwd": "${workspaceFolder}",
      "console": "integratedTerminal"
    },
    {
      "name": "Win MPI: debug semba-fdtd, 1 process (dbg-mpi)",
      "type": "cppvsdbg",
      "request": "launch",
      "program": "${workspaceFolder}/build-dbg-mpi/bin/semba-fdtd.exe",
      "args": ["-i", "pw-in-box.fdtd.json"],
      "cwd": "${workspaceFolder}/testData/cases/planewave",
      "console": "integratedTerminal"
    },
    {
      "name": "Win MPI: attach to a semba-fdtd process",
      "type": "cppvsdbg",
      "request": "attach",
      "processId": "${command:pickProcess}"
    }
  ]
}
```

Running a case writes its outputs next to the input file. Copy the case to a scratch folder if you do not want them in `testData/`.

To debug an MPI run with several processes:

1. Start the run from the integrated terminal with the debug build:

   ```bat
   set SEMBA_BUILD_DIR=build-dbg-mpi
   scripts\windows\run-mpi.bat 2 -i <case>.fdtd.json
   ```

2. Use *Win MPI: attach to a semba-fdtd process* and pick the rank's `semba-fdtd.exe`.

### Debugging with Visual Studio

1. Install the [prerequisites](#prerequisites) with the Visual Studio 2022 IDE (for example Community) instead of the Build Tools.
   Enable the Intel Visual Studio integration with `-p=NEED_VS2022_INTEGRATION=1`, or select it in the graphical installer.
2. Open a `cmd` console and load the environment with `call scripts\windows\env.bat`.
3. Start Visual Studio from that console with `devenv .`, so that it opens the repository folder and inherits the environment.
4. Select a configure preset, for example `dbg`, and build.
5. If you hit [CMake issue 21816](https://gitlab.kitware.com/cmake/cmake/-/issues/21816) with the Visual Studio integration, add
   `add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:/nofpp>)` to `CMakeLists.txt`.
   The Ninja command-line and VS Code builds described above do not need it.
6. **Debugger configuration:**
   - Select your build as the startup item.
   - In the debug options, open `Debug and launch settings` for your startup item.
   - Set **Working Directory** to the folder of the case you want to debug.
   - Set **Command Arguments** to `-i filename`, where `filename` is the input file.

   You should end up with something like this:

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

7. Press **F5** to start debugging.

### Windows troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `"vars.bat" is not recognized as an internal or external command` while loading oneAPI | `NoDefaultCurrentDirectoryInExePath` is set, so `cmd` will not run batch files from the current folder | Clear it (`set NoDefaultCurrentDirectoryInExePath=`). `env.bat` does this. |
| CMake: *"The Fortran compiler identification is unknown"*, with `link: unknown option -- s` in `CMakeConfigureLog.yaml` | Git for Windows `usr\bin\link.exe` is found before the MSVC linker | Load `vcvars64.bat` (`env.bat`), which puts the MSVC linker first, or remove `Git\usr\bin` from `PATH`. |
| `ifx: error #10037: could not find 'link'` | The MSVC environment is not loaded. `setvars.bat` alone does not load it | Call `vcvars64.bat` before `setvars.bat`, or use `env.bat`. |
| CMake: *"Missing variable is: CMAKE_Fortran_PREPROCESS_SOURCE"* | Usually one of the two errors above. With CMake 4.x it also happens with a correct environment | Fix the environment, and use CMake 3.x. |
| `ninja: error: ...F90-pp.f90.d: no outputs declared` | An interrupted build left truncated dependency files | Delete the affected `.d` files or the build folder, and build again. |
| Intel installer: *"Cannot install ... It is already installed"* | The product is already installed and you want to add a component | Use `--action modify` with the full list of components. |
| `semba-fdtd.exe` fails to start because of missing DLLs (`libifcoremd.dll`, `impi.dll`, ...) | oneAPI runtime not on `PATH` | Run it through `run.bat`/`run-mpi.bat`, or load `env.bat` first. |

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
