#!/usr/bin/env bash
# run_benchmarks.sh — Run semba-fdtd benchmarks and report results
# Usage: ./run_benchmarks.sh [cuda-fortran|cpu]
#   cuda-forfan: uses build-cuf-prof/bin/semba-fdtd with SEMBA_FDTD_ENABLE_CUF_RUNTIME=1
#   cpu: uses build-rls/bin/semba-fdtd (no GPU)
#   default: cuda-forfan

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

BINARY="${1:-cuda-forfan}"

case "$BINARY" in
  cuda-forfan)
    SEMBA_BIN="$ROOT_DIR/build-cuf-prof/bin/semba-fdtd"
    export SEMBA_FDTD_ENABLE_CUF_RUNTIME=1
    LABEL="CUDA-Fortran (RTX 5080)"
    ;;
  cpu)
    SEMBA_BIN="$ROOT_DIR/build-cpu-rls/bin/semba-fdtd"
    unset SEMBA_FDTD_ENABLE_CUF_RUNTIME
    LABEL="CPU (48-core Zen 4)"
    ;;
  *)
    echo "Usage: $0 [cuda-forfan|cpu]"
    exit 1
    ;;
esac

if [[ ! -x "$SEMBA_BIN" ]]; then
    echo "ERROR: Binary not found: $SEMBA_BIN"
    exit 1
fi

# Test cases: name, directory, input file, cells, steps
declare -a CASES=(
    "nodalSource:testData/cases/nodalSource:nodalSource.fdtd.json:445K:9000"
    "towelHanger:testData/cases/towelHanger:towelHanger.fdtd.json:216K:2000"
    "multipleAssigments:testData/cases/multipleAssigments:multipleDielectricMaterial.fdtd.json:500K:500"
  "sphere:testData/cases/sphere:sphere.fdtd.json:512K:100"
)

echo "=============================================="
echo " SEMBA-FDTD Benchmark Suite"
echo "=============================================="
echo " Binary : $SEMBA_BIN"
echo " Label  : $LABEL"
echo " GPU    : ${SEMBA_FDTD_ENABLE_CUF_RUNTIME:-off}"
echo " Date   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="
echo ""

# Header
printf "%-20s %8s %8s %10s %8s\n" "Case" "Cells" "Steps" "Time (s)" "Speedup"
printf "%-20s %8s %8s %10s %8s\n" "--------------------" "--------" "--------" "----------" "--------"

# CPU baseline times (from previous measurements)
declare -A CPU_BASELINE=(
    ["nodalSource"]="35.20"
    ["towelHanger"]="7.96"
    ["multipleAssigments"]="2.32"
    ["sphere"]="3.26"
)

TOTAL_TIME=0
TOTAL_STEPS=0

for entry in "${CASES[@]}"; do
    IFS=':' read -r NAME DIR FILE CELLS STEPS <<< "$entry"
    CASE_DIR="$ROOT_DIR/$DIR"

    if [[ ! -f "$CASE_DIR/$FILE" ]]; then
        echo "SKIP: $NAME (input not found: $CASE_DIR/$FILE)"
        continue
    fi

    # Clean output files
    rm -f "$CASE_DIR"/*.dat "$CASE_DIR"/*.bin "$CASE_DIR"/*.h5 "$CASE_DIR"/*.xdmf "$CASE_DIR"/*.vtk \
          "$CASE_DIR"/*.Report.txt "$CASE_DIR"/*.Warnings.txt "$CASE_DIR"/*.fdtd_Report.txt \
          "$CASE_DIR"/*.fdtd_tmpWarnings.txt "$CASE_DIR"/*.old 2>/dev/null || true

    # Run benchmark (3 iterations, take best)
    BEST_TIME=999999
    for i in 1 2 3; do
        start=$(date +%s%N)
        (cd "$CASE_DIR" && "$SEMBA_BIN" -i "$FILE" > /dev/null 2>&1)
        end=$(date +%s%N)
        elapsed=$(echo "scale=3; ($end - $start) / 1000000000" | bc)
        if (( $(echo "$elapsed < $BEST_TIME" | bc -l) )); then
            BEST_TIME=$elapsed
        fi
    done

    TOTAL_TIME=$(echo "$TOTAL_TIME + $BEST_TIME" | bc)
    TOTAL_STEPS=$((TOTAL_STEPS + STEPS))

    # Calculate speedup vs CPU baseline
    BASELINE="${CPU_BASELINE[$NAME]:-0}"
    if (( $(echo "$BASELINE > 0" | bc -l) )); then
        SPEEDUP=$(echo "scale=2; $BASELINE / $BEST_TIME" | bc)
    else
        SPEEDUP="N/A"
    fi

    printf "%-20s %8s %8s %10s %8s\n" "$NAME" "$CELLS" "$STEPS" "$BEST_TIME" "$SPEEDUP"x

done

echo ""
echo "----------------------------------------------"
echo " Total time : ${TOTAL_TIME}s"
echo " Total steps: $TOTAL_STEPS"
echo "----------------------------------------------"
echo ""
echo "Baseline times (CPU, 48-core Zen 4):"
for name in "${!CPU_BASELINE[@]}"; do
    echo "  $name: ${CPU_BASELINE[$name]}s"
done
echo ""
echo "=============================================="
