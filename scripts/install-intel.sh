#!/usr/bin/env bash
# install-intel.sh
# Installs Intel oneAPI Fortran compiler 2025.1 and Intel MPI,
# mirroring what fortran-lang/setup-fortran@v1 and mpi4py/setup-mpi@v1 do.

set -euo pipefail

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

echo "=== Done ==="
echo "Source the environment before building:"
echo "  source /opt/intel/oneapi/setvars.sh"
