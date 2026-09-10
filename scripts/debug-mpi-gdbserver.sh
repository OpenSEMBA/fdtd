#!/usr/bin/env bash
set -euo pipefail

state_dir="${TMPDIR:-/tmp}/semba-fdtd-mpi-gdbserver-${UID}"
pid_file="${state_dir}/pid"
log_file="${state_dir}/output.log"

print_log() {
    if [[ -f "$log_file" ]]; then
        while IFS= read -r line; do
            printf '%s\n' "$line" >&2
        done < "$log_file"
    fi
}

stop_job() {
    local pid cmdline

    if [[ ! -f "$pid_file" ]]; then
        return
    fi

    IFS= read -r pid < "$pid_file"
    if [[ "$pid" =~ ^[1-9][0-9]*$ ]] && kill -0 "$pid" 2>/dev/null; then
        cmdline="$(tr '\0' ' ' < "/proc/${pid}/cmdline" 2>/dev/null || true)"
        if [[ "$cmdline" != *mpirun* && "$cmdline" != *mpiexec* && "$cmdline" != *prterun* ]]; then
            printf 'Refusing to stop unexpected process %s: %s\n' "$pid" "$cmdline" >&2
            return 1
        fi

        kill "$pid" 2>/dev/null || true
        for _ in {1..30}; do
            if ! kill -0 "$pid" 2>/dev/null; then
                break
            fi
            sleep 0.1
        done
        if kill -0 "$pid" 2>/dev/null; then
            kill -KILL "$pid" 2>/dev/null || true
        fi
    fi

    rm -f "$pid_file"
}

port_is_listening() {
    local port_hex table slot local_address remote_address connection_state remainder

    printf -v port_hex '%04X' "$1"
    for table in /proc/net/tcp /proc/net/tcp6; do
        [[ -r "$table" ]] || continue
        while read -r slot local_address remote_address connection_state remainder; do
            if [[ "${local_address##*:}" == "$port_hex" && "$connection_state" == "0A" ]]; then
                return 0
            fi
        done < "$table"
    done
    return 1
}

mkdir -p "$state_dir"

if [[ "${1:-}" == "--stop" ]]; then
    stop_job
    exit 0
fi

if [[ "${1:-}" == "--workdir" ]]; then
    workdir="${2:-}"
    if [[ ! -d "$workdir" ]]; then
        printf 'MPI debug working directory does not exist: %s\n' "$workdir" >&2
        exit 2
    fi
    if [[ ! -w "$workdir" ]]; then
        printf 'MPI debug working directory is not writable: %s\n' "$workdir" >&2
        exit 2
    fi
    cd "$workdir"
    shift 2
fi

if [[ "${1:-}" == "--wait-for-port" ]]; then
    port="${2:-}"
    if ! [[ "$port" =~ ^[0-9]+$ ]] || ((port < 1 || port > 65535)); then
        printf 'Invalid port: %s\n' "$port" >&2
        exit 2
    fi

    for _ in {1..300}; do
        if port_is_listening "$port"; then
            printf 'MPI gdbserver ready on port %s\n' "$port"
            trap 'exit 0' INT TERM
            while true; do
                sleep 3600
            done
        fi
        sleep 0.1
    done

    printf 'Timed out waiting for MPI gdbserver on port %s.\n' "$port" >&2
    exit 1
fi

debug_rank=""
foreground=false
foreground_all=false
if [[ "${1:-}" == "--foreground-all" ]]; then
    foreground_all=true
    shift
elif [[ "${1:-}" == "--foreground-debug-rank" ]]; then
    foreground=true
    if [[ "$#" -lt 3 ]]; then
        printf 'Missing rank after --foreground-debug-rank.\n' >&2
        exit 2
    fi
    debug_rank="$2"
    shift 2
elif [[ "${1:-}" == "--debug-rank" ]]; then
    if [[ "$#" -lt 3 ]]; then
        printf 'Missing rank after --debug-rank.\n' >&2
        exit 2
    fi
    debug_rank="$2"
    shift 2
fi

if [ "$#" -lt 3 ]; then
    printf 'Usage: %s <ranks> <program> [arguments ...]\n' "$0" >&2
    printf '       %s --workdir <directory> <mode and arguments ...>\n' "$0" >&2
    printf '       %s --debug-rank <rank> <ranks> <program> [arguments ...]\n' "$0" >&2
    printf '       %s --foreground-all <ranks> <program> [arguments ...]\n' "$0" >&2
    printf '       %s --foreground-debug-rank <rank> <ranks> <program> [arguments ...]\n' "$0" >&2
    printf '       %s --wait-for-port <port>\n' "$0" >&2
    printf '       %s --stop\n' "$0" >&2
    exit 2
fi

ranks="$1"
program="$2"
shift 2

if ! [[ "$ranks" =~ ^[1-9][0-9]*$ ]]; then
    printf 'The number of MPI ranks must be a positive integer: %s\n' "$ranks" >&2
    exit 2
fi

if [[ -n "$debug_rank" ]] && {
    ! [[ "$debug_rank" =~ ^[0-9]+$ ]] || ((debug_rank >= ranks));
}; then
    printf 'The debug rank must be between 0 and %s: %s\n' "$((ranks - 1))" "$debug_rank" >&2
    exit 2
fi

if [[ ! -x "$program" ]]; then
    printf 'MPI debug executable is missing or not executable: %s\n' "$program" >&2
    exit 2
fi

stop_job

if [[ "$foreground_all" == true ]]; then
    printf '%s\n' "$$" > "$pid_file"
    printf 'MPI working directory: %s\n' "$PWD"
    printf 'Starting %s MPI gdbserver ranks on ports 20000-%s\n' \
        "$ranks" "$((20000 + ranks - 1))"
    exec mpirun --tag-output -np "$ranks" bash -c '
        port=$((20000 + OMPI_COMM_WORLD_RANK))
        exec gdbserver --no-disable-randomization --once ":${port}" "$@"
    ' gdbserver-rank "$program" "$@"
fi

if [[ "$foreground" == true ]]; then
    printf '%s\n' "$$" > "$pid_file"
    printf 'MPI working directory: %s\n' "$PWD"
    printf 'Starting %s MPI ranks with rank %s under gdbserver on port 20000\n' \
        "$ranks" "$debug_rank"
    exec mpirun --tag-output -np "$ranks" bash -c '
        debug_rank=$1
        shift
        if [[ "$OMPI_COMM_WORLD_RANK" == "$debug_rank" ]]; then
            exec gdbserver --no-disable-randomization --once :20000 "$@"
        fi
        exec "$@"
    ' mpi-debug-rank "$debug_rank" "$program" "$@"
fi

: > "$log_file"

if [[ -n "$debug_rank" ]]; then
    printf 'Starting %s MPI ranks with rank %s under gdbserver on port 20000\n' \
        "$ranks" "$debug_rank"
    expected_servers=1
    nohup mpirun --tag-output -np "$ranks" bash -c '
        debug_rank=$1
        shift
        if [[ "$OMPI_COMM_WORLD_RANK" == "$debug_rank" ]]; then
            exec gdbserver --no-disable-randomization --once :20000 "$@"
        fi
        exec "$@"
    ' mpi-debug-rank "$debug_rank" "$program" "$@" > "$log_file" 2>&1 < /dev/null &
else
    printf 'Starting %s MPI gdbserver ranks on ports 20000-%s\n' \
        "$ranks" "$((20000 + ranks - 1))"
    expected_servers="$ranks"
    nohup mpirun --tag-output -np "$ranks" bash -c '
        port=$((20000 + OMPI_COMM_WORLD_RANK))
        exec gdbserver --no-disable-randomization --once ":${port}" "$@"
    ' gdbserver-rank "$program" "$@" > "$log_file" 2>&1 < /dev/null &
fi
mpi_pid=$!
printf '%s\n' "$mpi_pid" > "$pid_file"

for _ in {1..150}; do
    ready=0
    while IFS= read -r line; do
        if [[ "$line" == *"Listening on port "* ]]; then
            ((ready += 1))
        fi
    done < "$log_file"

    if ((ready >= expected_servers)); then
        printf 'MPI gdbserver is ready.\n'
        exit 0
    fi

    if ! kill -0 "$mpi_pid" 2>/dev/null; then
        printf 'MPI gdbserver exited before all ranks were ready.\n' >&2
        print_log
        rm -f "$pid_file"
        exit 1
    fi

    sleep 0.1
done

printf 'Timed out waiting for all MPI gdbserver ranks.\n' >&2
print_log
stop_job
exit 1
