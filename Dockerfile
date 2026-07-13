# syntax=docker/dockerfile:1.7

# ─── Base ─────────────────────────────────────────────────────────────────────
FROM ubuntu:26.04@sha256:b7f48194d4d8b763a478a621cdc81c27be222ba2206ca3ca6bc42b49685f3d9e AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    locales \
    gfortran \
    g++ \
    cmake \
    make \
    ninja-build \
    libhdf5-dev \
    libopenmpi-dev \
    python3 \
    python3-pip \
    python3-venv \
    gdb \
    gdbserver \
    && locale-gen en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*

ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

# ─── Stage 1: Builder ─────────────────────────────────────────────────────────
FROM base AS builder

WORKDIR /src
COPY . .

# Build (MPI off by default; override at build time with --build-arg ENABLE_MPI=ON)
ARG ENABLE_MPI=OFF
ARG ENABLE_MTLN=ON
ARG BUILD_TYPE=Release
ARG ENABLE_TEST=OFF

RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DSEMBA_FDTD_ENABLE_MPI=${ENABLE_MPI} \
        -DSEMBA_FDTD_ENABLE_HDF=ON \
        -DSEMBA_FDTD_ENABLE_MTLN=${ENABLE_MTLN} \
        -DSEMBA_FDTD_ENABLE_SMBJSON=ON \
        -DSEMBA_FDTD_ENABLE_TEST=${ENABLE_TEST} \
    && cmake --build build -j$(nproc) \
    && cp build/bin/semba-fdtd /tmp/semba-fdtd \
    && if [ "${BUILD_TYPE}" = "Release" ]; then strip /tmp/semba-fdtd; fi

# ─── Stage 1b: Test Builder ───────────────────────────────────────────────────
FROM builder AS test-builder

RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DSEMBA_FDTD_ENABLE_MPI=${ENABLE_MPI} \
        -DSEMBA_FDTD_ENABLE_HDF=ON \
        -DSEMBA_FDTD_ENABLE_MTLN=${ENABLE_MTLN} \
        -DSEMBA_FDTD_ENABLE_SMBJSON=ON \
        -DSEMBA_FDTD_ENABLE_TEST=ON \
    && cmake --build build -j$(nproc)

# Install Python test/wrapper dependencies only for images that run tests.
RUN python3 -m venv /opt/fdtd-venv \
    && /opt/fdtd-venv/bin/python -m pip install --upgrade pip \
    && /opt/fdtd-venv/bin/python -m pip install --no-cache-dir -r requirements.txt

ENV PATH=/opt/fdtd-venv/bin:${PATH}

# ─── Stage 2: Development Tools ────────────────────────────────────────────────
FROM base AS dev-tools

USER root

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    git \
    gh \
    openssh-client \
    nodejs \
    npm \
    paraview \
    sudo \
    && rm -rf /var/lib/apt/lists/*

RUN --mount=type=cache,target=/root/.npm \
    npm install -g opencode-ai

RUN --mount=type=cache,target=/root/.cache/pip \
    python3 -m pip install --break-system-packages fortls fprettify

# ─── Stage 2b: Development User ────────────────────────────────────────────────
FROM dev-tools AS dev

ARG USERNAME=developer
ARG USER_UID=1000
ARG USER_GID=1000

ENV USERNAME=${USERNAME} \
    USER_UID=${USER_UID} \
    USER_GID=${USER_GID}

RUN userdel --remove ubuntu 2>/dev/null || true \
    && groupdel ubuntu 2>/dev/null || true \
    && groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME} \
    && mkdir -p /home/${USERNAME}/.config \
    && mkdir -p /home/${USERNAME}/.ssh \
    && mkdir -p /home/${USERNAME}/.local/share/opencode \
    && chmod 700 /home/${USERNAME}/.ssh \
    && chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.config \
    && chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.ssh \
    && chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.local \
    && usermod -aG sudo ${USERNAME} \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" >/etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME}

USER ${USERNAME}
WORKDIR /home/${USERNAME}/workspaces/fdtd

# ─── Stage 3: Runtime ─────────────────────────────────────────────────────────
# Minimal image with only the shared libraries the binary needs at runtime.
# HDF5, LAPACK, BLAS, and ngspice are all statically linked on Linux.
# Keep OpenMPI runtime libraries available for MPI-enabled image builds.
FROM ubuntu:26.04@sha256:b7f48194d4d8b763a478a621cdc81c27be222ba2206ca3ca6bc42b49685f3d9e AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    libgfortran5 \
    libgomp1 \
    libcurl4t64 \
    zlib1g \
    libaec0 \
    libsz2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /tmp/semba-fdtd /usr/local/bin/semba-fdtd

WORKDIR /work

ENTRYPOINT ["semba-fdtd", "-i"]
