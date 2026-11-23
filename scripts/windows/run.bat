@echo off
REM Run semba-fdtd from the current folder (usually the case folder) with the oneAPI runtime loaded.
REM
REM Usage: scripts\windows\run.bat -i CASE.fdtd.json [other semba-fdtd arguments]
REM   SEMBA_BUILD_DIR selects the build folder (default: build-rls).
setlocal
for %%I in ("%~dp0..\..") do set "SEMBA_ROOT=%%~fI"
if "%SEMBA_BUILD_DIR%"=="" set SEMBA_BUILD_DIR=build-rls

call "%~dp0env.bat" || exit /b 1

"%SEMBA_ROOT%\%SEMBA_BUILD_DIR%\bin\semba-fdtd.exe" %*
