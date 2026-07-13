#!/usr/bin/env bash
set -u

expected_home="/home/developer"

warn() {
    printf 'devcontainer warning: %s\n' "$*" >&2
}

if [ "$(id -un)" != "developer" ]; then
    warn "running as $(id -un), expected developer"
fi

if [ "${HOME:-}" != "$expected_home" ]; then
    warn "HOME is ${HOME:-unset}, expected $expected_home"
fi

if [ "${XDG_CONFIG_HOME:-}" != "$expected_home/.config" ]; then
    warn "XDG_CONFIG_HOME is ${XDG_CONFIG_HOME:-unset}, expected $expected_home/.config"
fi

if [ "${XDG_DATA_HOME:-}" != "$expected_home/.local/share" ]; then
    warn "XDG_DATA_HOME is ${XDG_DATA_HOME:-unset}, expected $expected_home/.local/share"
fi

if ! command -v opencode >/dev/null 2>&1; then
    warn "opencode is not available in PATH"
fi

if ! command -v paraview >/dev/null 2>&1 && ! command -v pvpython >/dev/null 2>&1; then
    warn "ParaView is not available in PATH; rebuild the dev container image"
fi

if command -v git >/dev/null 2>&1; then
    git config --global --add safe.directory /home/developer/workspaces/fdtd >/dev/null 2>&1 || \
        warn "could not add workspace to git safe.directory"
fi

if [ ! -d "$expected_home/.config/opencode" ]; then
    warn "$expected_home/.config/opencode is not mounted"
fi

if [ ! -d "$expected_home/.agents/skills" ]; then
    warn "$expected_home/.agents/skills is not mounted"
fi

if [ ! -d "$expected_home/.local/share/opencode" ]; then
    warn "$expected_home/.local/share/opencode is not mounted"
fi

if [ ! -f "$expected_home/.local/share/opencode/auth.json" ]; then
    warn "OpenCode provider auth is missing at $expected_home/.local/share/opencode/auth.json"
    warn "run opencode /connect once on the host or in this container if providers are unavailable"
fi

exit 0
