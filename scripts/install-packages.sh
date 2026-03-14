#!/bin/bash
# Install script for ugrfdtd development environment
# Generated for Ubuntu 24.04 LTS (Noble Numbat) x86_64

set -e

echo "=== Updating package lists ==="
sudo apt update

echo "=== Installing build tools ==="
sudo apt install -y \
    build-essential \
    cmake \
    ninja-build \
    git \
    gdb \
    pkg-config

echo "=== Installing Fortran compiler ==="
sudo apt install -y \
    gfortran

echo "=== Installing MPI ==="
sudo apt install -y \
    libopenmpi-dev \
    openmpi-bin \
    openmpi-common

echo "=== Installing HDF5 ==="
sudo apt install -y \
    libhdf5-dev \
    hdf5-helpers

echo "=== Installing Python development tools ==="
sudo apt install -y \
    python3-dev \
    python3-pip


echo "=== Installing Python packages ==="
pip3 install --break-system-packages fortls

COMPILER_VERSION="2025.1"

echo "=== Adding Intel oneAPI apt repository ==="
KEY="GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB"
curl -fsSL "https://apt.repos.intel.com/intel-gpg-keys/$KEY" \
    | sudo apt-key add -
echo "deb https://apt.repos.intel.com/oneapi all main" \
    | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt-get update -q

echo "=== Installing Intel Fortran + C/C++ compilers (${COMPILER_VERSION}) ==="
sudo apt-get install -y \
    intel-oneapi-compiler-fortran-${COMPILER_VERSION} \
    intel-oneapi-compiler-dpcpp-cpp-${COMPILER_VERSION}

echo "=== Installing Intel MPI ==="
sudo apt-get install -y \
    intel-oneapi-mpi \
    intel-oneapi-mpi-devel

echo ""
echo "=== All packages installed successfully ==="
echo "Source the Intel environment before building:"
echo "  source /opt/intel/oneapi/setvars.sh"
