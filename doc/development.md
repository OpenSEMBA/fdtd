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
git submodule init
git submodule update
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
            "program": "${workspaceRoot}/build/bin/semba-fdtd",
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

Instructions following [this link](https://github.com/microsoft/vscode-cpptools/issues/1723#issuecomment-1983158512).

1. First install gdbserver:

```bash
sudo apt install gdbserver
```

2.  In a terminal in vscode execute this command to have gdbserver set a series of debugging servers, that will hold listening until a process is attached to them:
```bash
mpiexec -n 1 gdbserver :20000 ./build/bin/semba-fdtd : -n 1 gdbserver :20001 .build/bin/semba-fdtd
```
This example assumes we are `${workspaceRoot}` folder, hence the path to the semba-fdtd executable

3. Modify the file launch.json to connect to the debugging servers when launched:

```json
{
    "version": "0.2.0",
    "compounds": [
        {
          "name": "connect F (NProcs=2)",
          "configurations": [
            "_rank0",
            "_rank1"
          ],
          "stopAll": true
        }
      ],
      "configurations": [
        {
          "name": "_rank0",
          "type": "cppdbg",
          "request": "launch",
          "program": "${workspaceRoot}/build/bin/semba-fdtd",
          "programArgs" :  "-i holland1981.fdtd.json -mtlnwires",
          "miDebuggerServerAddress": "localhost:20000",
          "MIMode": "gdb",
          "miDebuggerPath": "/usr/bin/gdb",
          "cwd": "${workspaceRoot}/testData/cases/holland_mpi"
        },
        {
          "name": "_rank1",
          "type": "cppdbg",
          "request": "launch",
          "program": "${workspaceRoot}/build/bin/semba-fdtd",
          "programArgs" :  "-i holland1981.fdtd.json -mtlnwires",
          "miDebuggerServerAddress": "localhost:20001",
          "MIMode": "gdb",
          "miDebuggerPath": "/usr/bin/gdb",
          "cwd": "${workspaceRoot}/testData/cases/holland_mpi"
        }
      ]
}
```