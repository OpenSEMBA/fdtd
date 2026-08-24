#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s <program> [arguments ...]\n' "$0" >&2
}

if [[ "$#" -eq 0 ]]; then
    usage
    exit 2
fi

if ! command -v pidstat >/dev/null 2>&1; then
    printf 'pidstat is required. Install it with: sudo apt install sysstat\n' >&2
    exit 127
fi

if [[ "$1" == */* ]] && [[ ! -x "$1" ]]; then
    printf 'Program is missing or not executable: %s\n' "$1" >&2
    exit 2
fi

log_dir="${PWD}/logs"
timestamp="$(date +%Y%m%d-%H%M%S)"
log_file="${log_dir}/resource-usage-${timestamp}.log"
program_log="${log_dir}/program-output-${timestamp}.log"

mkdir -p "$log_dir"

printf 'Writing resource usage to %s\n' "$log_file"
printf 'Writing program output to %s\n' "$program_log"
PROGRAM_LOG="$program_log" pidstat -H -h -u -r -d 1 \
    -e bash -c 'exec "$@" > >(tee "$PROGRAM_LOG" >&2) 2>&1' monitor-resources "$@" \
    | tee "$log_file"
