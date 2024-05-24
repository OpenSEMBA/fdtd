REM SPDX-FileCopyrightText: 2022 Intel Corporation
REM
REM SPDX-License-Identifier: MIT

set LANGUAGE=%1
set VS_VER=%2
set SAMPLES_TAG=%3

IF "%VS_VER%"=="2017_build_tools" (
    @call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
) ELSE (
    IF "%VS_VER%"=="2019_build_tools" (
        @call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    ) ELSE (
          @call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%
    )
)

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=Release -DCompileWithMPI=NO FC=ifx CC=icx
cmake --build build -j
