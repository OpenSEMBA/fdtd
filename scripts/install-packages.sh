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

echo "=== Installing math/numeric libraries ==="
sudo apt install -y \
    libgmp-dev \
    libmpfr-dev \
    libexpat1-dev \
    libreadline-dev \
    libncurses-dev \
    libcurl4-openssl-dev

echo "=== Installing SLURM / job scheduler ==="
sudo apt install -y \
    slurm-client \
    slurmd \
    munge \
    libmunge-dev

echo "=== Installing system utilities ==="
sudo apt install -y \
    vim \
    htop \
    iotop \
    screen \
    nload \
    net-tools \
    ncat \
    traceroute \
    wget \
    curl \
    rsync \
    zsh \
    byobu \
    openssh-server \
    jq

echo "=== Installing Python packages ==="
pip3 install --break-system-packages fortls

echo ""
echo "=== Intel oneAPI (MPI + Fortran compiler) ==="
echo "Intel oneAPI is not installed via apt in this script."
echo "To install it, follow the official guide:"
echo "  https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html"
echo ""
echo "Packages needed:"
echo "  intel-oneapi-compiler-fortran"
echo "  intel-oneapi-compiler-dpcpp-cpp"
echo "  intel-oneapi-mpi"
echo "  intel-oneapi-mpi-devel"
echo ""
echo "After installing, source the environment with:"
echo "  source /opt/intel/oneapi/setvars.sh"

echo ""
echo "=== All packages installed successfully ==="
