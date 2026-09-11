# syntax=docker/dockerfile:1.7

# Keep shared build dependencies first so both environments reuse this layer.
FROM ubuntu:26.04@sha256:b7f48194d4d8b763a478a621cdc81c27be222ba2206ca3ca6bc42b49685f3d9e AS quality-base

# Prevent package installation from prompting during image builds.
ENV DEBIAN_FRONTEND=noninteractive

# Cache apt metadata and packages between BuildKit builds without retaining them
# in the final image layer.
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
    paraview \
    && locale-gen en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*

# Use a UTF-8 locale consistently for compiler, test, and terminal output.
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

# These can be overridden to match the host user and avoid bind-mount ownership
# mismatches when running through Docker Compose.
ARG USERNAME=developer
ARG USER_UID=1000
ARG USER_GID=1000

# Reuse a pre-existing UID/GID where Ubuntu provides one, otherwise create the
# requested account.
RUN set -eux; \
    existing_group="$(getent group ${USER_GID} | cut -d: -f1)"; \
    if [ -n "${existing_group}" ] && [ "${existing_group}" != "${USERNAME}" ] && ! getent group ${USERNAME} >/dev/null; then \
        groupmod -n ${USERNAME} ${existing_group}; \
    elif [ -z "${existing_group}" ]; then \
        groupadd --gid ${USER_GID} ${USERNAME}; \
    fi; \
    existing_user="$(getent passwd ${USER_UID} | cut -d: -f1)"; \
    if [ -n "${existing_user}" ] && [ "${existing_user}" != "${USERNAME}" ]; then \
        usermod -l ${USERNAME} -d /home/${USERNAME} -m ${existing_user}; \
    elif [ -z "${existing_user}" ]; then \
        useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME}; \
    fi; \
    mkdir -p /home/${USERNAME}/workspaces/fdtd; \
    chown -R ${USER_UID}:${USER_GID} /home/${USERNAME}/workspaces

# Minimal environment for compiling and running the project.
FROM quality-base AS quality

ENV HOME=/home/${USERNAME}

USER ${USERNAME}
WORKDIR /home/${USERNAME}/workspaces/fdtd

# Full development environment. This target deliberately follows quality so its
# build reuses every compiler and ParaView layer.
FROM quality AS dev

USER root

# Development-only tools: version control, debugging, and the OpenCode CLI.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    git \
    gh \
    openssh-client \
    nodejs \
    npm \
    gdb \
    gdbserver \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Preserve package-manager caches across BuildKit builds.
RUN --mount=type=cache,target=/root/.npm \
    npm install -g opencode-ai

RUN --mount=type=cache,target=/root/.cache/pip \
    python3 -m pip install --break-system-packages fortls fprettify

# Prepare mount points and grant passwordless sudo for interactive development.
RUN mkdir -p /home/${USERNAME}/.config \
    /home/${USERNAME}/.ssh \
    /home/${USERNAME}/.local/share/opencode \
    && chmod 700 /home/${USERNAME}/.ssh \
    && chown -R ${USER_UID}:${USER_GID} /home/${USERNAME}/.config \
    /home/${USERNAME}/.ssh \
    /home/${USERNAME}/.local \
    && usermod -aG sudo ${USERNAME} \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" >/etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME}

USER ${USERNAME}
WORKDIR /home/${USERNAME}/workspaces/fdtd

# Intel oneAPI 2025.0 ships the compilers, Intel MPI, and matching runtime
# libraries as a tested toolkit rather than relying on Intel's APT repository.
FROM intel/oneapi-hpckit:2025.0.0-0-devel-ubuntu24.04 AS intel-quality

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    locales \
    cmake \
    make \
    ninja-build \
    python3 \
    python3-pip \
    python3-venv \
    libegl1 \
    libgl1 \
    && locale-gen en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*

ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

ARG USERNAME=developer
ARG USER_UID=1000
ARG USER_GID=1000

RUN set -eux; \
    existing_group="$(getent group ${USER_GID} | cut -d: -f1)"; \
    if [ -n "${existing_group}" ] && [ "${existing_group}" != "${USERNAME}" ] && ! getent group ${USERNAME} >/dev/null; then \
        groupmod -n ${USERNAME} ${existing_group}; \
    elif [ -z "${existing_group}" ]; then \
        groupadd --gid ${USER_GID} ${USERNAME}; \
    fi; \
    existing_user="$(getent passwd ${USER_UID} | cut -d: -f1)"; \
    if [ -n "${existing_user}" ] && [ "${existing_user}" != "${USERNAME}" ]; then \
        usermod -l ${USERNAME} -d /home/${USERNAME} -m ${existing_user}; \
    elif [ -z "${existing_user}" ]; then \
        useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME}; \
    fi; \
    mkdir -p /home/${USERNAME}/workspaces/fdtd;

COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /home/${USERNAME}/workspaces/fdtd/.venv \
    && /home/${USERNAME}/workspaces/fdtd/.venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt \
    && rm /tmp/requirements.txt

COPY docker/intel-quality-entrypoint.sh /usr/local/bin/intel-quality-entrypoint
RUN chmod 0755 /usr/local/bin/intel-quality-entrypoint
COPY . /home/${USERNAME}/workspaces/fdtd
RUN chown -R ${USER_UID}:${USER_GID} /home/${USERNAME}/workspaces
USER ${USERNAME}
WORKDIR /home/${USERNAME}/workspaces/fdtd
ENV VIRTUAL_ENV=/home/${USERNAME}/workspaces/fdtd/.venv
ENV PATH=${VIRTUAL_ENV}/bin:${PATH}

ENTRYPOINT ["/usr/local/bin/intel-quality-entrypoint"]
CMD ["bash", "-l"]

# Intel development environment kept separate from the GNU-based dev target.
# The inherited entrypoint enables oneAPI and Intel MPI for every shell.
FROM intel-quality AS intel-dev

USER root

RUN apt-get update && apt-get install -y \
    git \
    gh \
    openssh-client \
    nodejs \
    npm \
    gdb \
    gdbserver \
    paraview \
    sudo \
    && rm -rf /var/lib/apt/lists/*

RUN --mount=type=cache,target=/root/.npm \
    npm install -g opencode-ai

RUN --mount=type=cache,target=/root/.cache/pip \
    python3 -m pip install --break-system-packages fortls fprettify

RUN mkdir -p /home/${USERNAME}/.config \
    /home/${USERNAME}/.ssh \
    /home/${USERNAME}/.local/share/opencode \
    && chmod 700 /home/${USERNAME}/.ssh \
    && chown -R ${USER_UID}:${USER_GID} /home/${USERNAME}/.config \
    /home/${USERNAME}/.ssh \
    /home/${USERNAME}/.local \
    && usermod -aG sudo ${USERNAME} \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" >/etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME}

ENV HOME=/home/${USERNAME}

USER ${USERNAME}
WORKDIR /home/${USERNAME}/workspaces/fdtd
