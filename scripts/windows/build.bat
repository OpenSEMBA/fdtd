@echo off
REM Configure and build semba-fdtd with one of the CMake presets (default: rls).
REM
REM Usage: scripts\windows\build.bat [PRESET]
REM   e.g. build.bat rls | build.bat dbg | build.bat rls-mpi | build.bat dbg-mpi
REM Output: build-PRESET\bin (see CMakePresets.json).
setlocal
for %%I in ("%~dp0..\..") do set "SEMBA_ROOT=%%~fI"
set "PRESET=%~1"
if "%PRESET%"=="" set PRESET=rls

call "%~dp0env.bat" || exit /b 1

cd /d "%SEMBA_ROOT%"
cmake --preset %PRESET% || exit /b 1
cmake --build --preset %PRESET% || exit /b 1
echo [build] preset %PRESET% built successfully.
