@echo off
title Sphere claudecorr - MPI 1 x OMP 20 - 6000 steps
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul 2>&1
if errorlevel 1 (
  echo Intel runtime initialization failed.
  exit /b 1
)
cd /d "%~dp0"
set "PATH=C:\Users\salva\work\elemwave\fdtd\build-rls-mpi-omp-fast-relocated\bin;%PATH%"
set "OMP_NUM_THREADS=20"
set "OMP_DYNAMIC=FALSE"
set "OMP_PLACES=cores"
set "OMP_PROC_BIND=spread"
set "OMP_DISPLAY_ENV=FALSE"
set "MKL_NUM_THREADS=1"
set "I_MPI_PIN=0"
echo Running fast-math sphere: 1 MPI rank, 20 OpenMP threads, 6000 steps.
echo Output directory: %CD%
"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin\mpiexec.exe" -n 1 "C:\Users\salva\work\elemwave\fdtd\build-rls-mpi-omp-fast-relocated\bin\semba-fdtd_omp_mpi_fast.exe" -i sphere.fdtd_claudecorr.json
set "SOLVER_EXIT=%ERRORLEVEL%"
echo %SOLVER_EXIT%>solver_exit_code.txt
echo Solver finished with exit code %SOLVER_EXIT%.
exit /b %SOLVER_EXIT%
