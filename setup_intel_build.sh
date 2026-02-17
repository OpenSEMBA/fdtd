#!/bin/bash
# ========================================================
# Script para limpiar CMake y configurar Intel 2025.3
# ========================================================

# --------------------------
# Configuración del proyecto
# --------------------------
SRC_DIR="/home/username/elemwave/fdtd"
BUILD_DIR="$SRC_DIR/build-dbg"

# --------------------------
# Carga Intel oneAPI 2025.3
# --------------------------
echo "Cargando Intel oneAPI 2025.1..."
source /opt/intel/oneapi/setvars.sh intel64 --force --compiler-version=2025.1

# --------------------------
# Definir compiladores Intel
# --------------------------
export CC=/opt/intel/oneapi/compiler/2025.1/bin/icx
export CXX=/opt/intel/oneapi/compiler/2025.1/bin/icpx
export FC=/opt/intel/oneapi/compiler/2025.1/bin/ifx

echo "Compiladores configurados:"
echo "  CC = $CC"
echo "  CXX = $CXX"
echo "  FC = $FC"

# --------------------------
# Limpiar build anterior
# --------------------------
if [ -d "$BUILD_DIR" ]; then
    echo "Eliminando cache de CMake y archivos previos..."
    rm -rf "$BUILD_DIR/CMakeCache.txt" "$BUILD_DIR/CMakeFiles"
else
    mkdir -p "$BUILD_DIR"
fi

# --------------------------
# Ejecutar CMake con Ninja
# --------------------------
echo "Ejecutando CMake con Intel..."
cmake -G Ninja \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_Fortran_COMPILER=$FC \
      -S "$SRC_DIR" \
      -B "$BUILD_DIR"

echo "¡Configuración completa! Ahora puedes compilar con:"
echo "  ninja -C $BUILD_DIR"
