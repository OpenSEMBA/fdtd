# Build context is the repository root (see compose.yml), so the solver is
# compiled from source with the same preset used for local MPI builds.

# Build stage: MPI toolchain, HDF5 and the solver sources.
FROM ubuntu:26.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    gfortran \
    cmake \
    ninja-build \
    libhdf5-dev \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# The build reads the commit from git_info.txt when .git is not in the context.
ARG GIT_COMMIT=UNKNOWN
RUN echo "GIT_COMMIT=\"${GIT_COMMIT}\"" > git_info.txt \
    && cmake --fresh --preset rls-mpi -DSEMBA_FDTD_ENABLE_TEST=OFF \
    && cmake --build --preset rls-mpi --target semba-fdtd

# Runtime stage: only the libraries the solver links against, plus the MPI and SSH tools.
FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    openssh-server \
    openmpi-bin \
    libhdf5-310 \
    libhdf5-fortran-310 \
    libgfortran5 \
    procps \
    iputils-ping \
    && rm -rf /var/lib/apt/lists/*

# Unprivileged user that runs mpirun and the solver. Ubuntu already ships uid 1000
# as "ubuntu"; rename it so the uid matches the owner of the host work/ directory.
RUN usermod -l mpi -d /home/mpi -m ubuntu \
    && mkdir -p /run/sshd

COPY --from=build /src/build-rls-mpi/bin/semba-fdtd /usr/local/bin/semba-fdtd
COPY docker/cluster/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod 0755 /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
