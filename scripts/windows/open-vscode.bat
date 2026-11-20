@echo off
REM Open Visual Studio Code on the repository with the build environment loaded, so that
REM CMake Tools, the integrated terminal and the debugger inherit MSVC, oneAPI and Intel MPI.
REM
REM Close every VS Code window first: a new window joins the running VS Code process
REM and does NOT inherit this environment.
for %%I in ("%~dp0..\..") do set "SEMBA_ROOT=%%~fI"
call "%~dp0env.bat" || (pause & exit /b 1)

set "SEMBA_CODE=%LOCALAPPDATA%\Programs\Microsoft VS Code\Code.exe"
if not exist "%SEMBA_CODE%" set "SEMBA_CODE=%ProgramFiles%\Microsoft VS Code\Code.exe"
if exist "%SEMBA_CODE%" (
  start "" "%SEMBA_CODE%" "%SEMBA_ROOT%"
) else (
  start "" code "%SEMBA_ROOT%"
)
