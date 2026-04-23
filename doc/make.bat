@ECHO OFF

set SOURCEDIR=.
set BUILDDIR=_build

if "%SPHINXBUILD%"=="" (
    set SPHINXBUILD=python -m sphinx
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR%
