@echo off
REM Load the semba-fdtd Windows build environment into the current cmd session:
REM   MSVC (link.exe, Windows SDK, bundled CMake and Ninja) + Intel oneAPI (ifx, icx, Intel MPI).
REM
REM Usage: call scripts\windows\env.bat
REM
REM Optional variables:
REM   ONEAPI_ROOT       oneAPI installation folder (default: %ProgramFiles(x86)%\Intel\oneAPI)
REM   SEMBA_TOOLS_PATH  folders prepended to PATH, e.g. a standalone CMake 3.x and Ninja

REM Intel's setvars.bat calls vars.bat through a relative path, which fails when this is set.
set NoDefaultCurrentDirectoryInExePath=

REM --- MSVC: ifx needs link.exe and the Windows SDK. setvars.bat does not load them.
if defined VSCMD_VER goto msvc_done
set "SEMBA_VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%SEMBA_VSWHERE%" (
  echo [env] ERROR: vswhere.exe not found. Install Visual Studio 2022 or its Build Tools with the C++ workload.
  exit /b 1
)
set "SEMBA_VSINSTALL="
for /f "usebackq delims=" %%i in (`"%SEMBA_VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "SEMBA_VSINSTALL=%%i"
if not defined SEMBA_VSINSTALL (
  echo [env] ERROR: no Visual Studio installation with the MSVC x64 tools was found.
  exit /b 1
)
call "%SEMBA_VSINSTALL%\VC\Auxiliary\Build\vcvars64.bat" >nul
if not defined VSCMD_VER (
  echo [env] ERROR: vcvars64.bat failed in "%SEMBA_VSINSTALL%".
  exit /b 1
)
:msvc_done

REM --- Intel oneAPI: ifx, icx and, if installed, Intel MPI.
if defined SETVARS_COMPLETED goto oneapi_done
if not defined ONEAPI_ROOT set "ONEAPI_ROOT=%ProgramFiles(x86)%\Intel\oneAPI"
if not exist "%ONEAPI_ROOT%\setvars.bat" (
  echo [env] ERROR: "%ONEAPI_ROOT%\setvars.bat" not found. Install the Intel oneAPI Fortran compiler or set ONEAPI_ROOT.
  exit /b 1
)
call "%ONEAPI_ROOT%\setvars.bat" intel64 >nul
if not defined SETVARS_COMPLETED (
  echo [env] ERROR: "%ONEAPI_ROOT%\setvars.bat" failed.
  exit /b 1
)
:oneapi_done

if defined SEMBA_TOOLS_PATH set "PATH=%SEMBA_TOOLS_PATH%;%PATH%"

set CC=icx
set CXX=icx
set FC=ifx

REM --- Sanity checks.
where ifx >nul 2>&1 || (echo [env] ERROR: ifx not found in the oneAPI environment. & exit /b 1)
where cmake >nul 2>&1 || (echo [env] ERROR: cmake not found. Install the "C++ CMake tools for Windows" VS component or set SEMBA_TOOLS_PATH. & exit /b 1)
where ninja >nul 2>&1 || (echo [env] ERROR: ninja not found. Install the "C++ CMake tools for Windows" VS component or set SEMBA_TOOLS_PATH. & exit /b 1)
set "SEMBA_CMAKE_VERSION="
for /f "tokens=3" %%v in ('cmake --version ^| findstr /b /c:"cmake version"') do set "SEMBA_CMAKE_VERSION=%%v"
if "%SEMBA_CMAKE_VERSION:~0,1%"=="4" echo [env] WARNING: CMake %SEMBA_CMAKE_VERSION% cannot detect the ifx Fortran preprocessor on Windows. Use CMake 3.x.
REM Git for Windows ships a coreutils link.exe that breaks ifx if it precedes the MSVC linker.
set "SEMBA_LINK="
for /f "delims=" %%l in ('where link 2^>nul') do if not defined SEMBA_LINK set "SEMBA_LINK=%%l"
if "%SEMBA_LINK:\VC\Tools\MSVC\=%"=="%SEMBA_LINK%" echo [env] WARNING: first link.exe on PATH is "%SEMBA_LINK%", not the MSVC linker.
exit /b 0
