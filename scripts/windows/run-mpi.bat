@echo off
REM Run semba-fdtd with Intel MPI from the current folder (usually the case folder).
REM
REM Usage: scripts\windows\run-mpi.bat NUM_PROCESSES -i CASE.fdtd.json [other semba-fdtd arguments]
REM   SEMBA_BUILD_DIR selects the build folder (default: build-rls-mpi).
setlocal
REM Keep these before "shift", which also shifts %0.
set "SEMBA_SCRIPTS=%~dp0"
for %%I in ("%~dp0..\..") do set "SEMBA_ROOT=%%~fI"
set "NP=%~1"
if "%NP%"=="" (
  echo Usage: run-mpi.bat NUM_PROCESSES -i CASE.fdtd.json
  exit /b 1
)
shift
set ARGS=
:collect
if "%~1"=="" goto launch
set ARGS=%ARGS% %1
shift
goto collect
:launch
if "%SEMBA_BUILD_DIR%"=="" set SEMBA_BUILD_DIR=build-rls-mpi

call "%SEMBA_SCRIPTS%env.bat" || exit /b 1

mpiexec -n %NP% "%SEMBA_ROOT%\%SEMBA_BUILD_DIR%\bin\semba-fdtd.exe"%ARGS%
