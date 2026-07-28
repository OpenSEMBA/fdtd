#!/usr/bin/env bash
set -e

# Initialise compiler and Intel MPI paths for every interactive container.
source /opt/intel/oneapi/setvars.sh --force >/dev/null

export HDF5_ROOT="${HDF5_ROOT:-/home/developer/workspaces/fdtd/precompiled_libraries/linux-intel/hdf5}"
export LD_LIBRARY_PATH="${HDF5_ROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

exec "$@"
