# ─── Stage 1: Builder ─────────────────────────────────────────────────────────
FROM ubuntu:22.04@sha256:eb29ed27b0821dca09c2e28b39135e185fc1302036427d5f4d70a41ce8fd7659 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    gfortran \
    g++ \
    cmake \
    make \
    libhdf5-dev \
    libopenmpi-dev \
    python3 \
    python3-pip \
    gdb \
    gdbserver \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# Install Python test/wrapper dependencies
RUN python3 -m pip install --no-cache-dir -r requirements.txt

# Build (MPI off by default; override at build time with --build-arg ENABLE_MPI=ON)
ARG ENABLE_MPI=OFF
ARG ENABLE_MTLN=ON
ARG BUILD_TYPE=Release

RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DSEMBA_FDTD_ENABLE_MPI=${ENABLE_MPI} \
        -DSEMBA_FDTD_ENABLE_HDF=ON \
        -DSEMBA_FDTD_ENABLE_MTLN=${ENABLE_MTLN} \
        -DSEMBA_FDTD_ENABLE_SMBJSON=ON \
        -DSEMBA_FDTD_ENABLE_TEST=ON \
    && cmake --build build -j$(nproc)

# ─── Stage 2: Runtime ─────────────────────────────────────────────────────────
# Minimal image with only the shared libraries the binary needs at runtime.
# HDF5, LAPACK, BLAS, and ngspice are all statically linked on Linux,
# so only the Fortran/OpenMP runtimes and HDF5 transitive deps are required.
FROM ubuntu:22.04@sha256:eb29ed27b0821dca09c2e28b39135e185fc1302036427d5f4d70a41ce8fd7659 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    libgfortran5 \
    libgomp1 \
    zlib1g \
    libaec2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/build/bin/semba-fdtd /usr/local/bin/semba-fdtd

WORKDIR /work

ENTRYPOINT ["semba-fdtd", "-i"]
