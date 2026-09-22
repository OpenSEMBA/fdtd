#!/usr/bin/env bash
# CPU (OMP=1) vs CUDA probe comparison for the P0/P1 CUDA slice.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EXE="${SEMBA_EXE:-$ROOT/build-rls-cuda/bin/semba-fdtd}"
CASE="${1:-/tmp/fdtd-cuda-profile/box_pml_smallpw_50.fdtd.json}"
WORKDIR="${TMPDIR:-/tmp}/fdtd-cuda-validate-$$"
export PATH="/usr/bin:/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu:/usr/local/cuda-13.3/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

if [[ ! -x "$EXE" ]]; then
  echo "Missing CUDA binary: $EXE (build with cmake --preset rls-cuda)" >&2
  exit 1
fi
if ! nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -q .; then
  : # free GPU
else
  echo "WARNING: other CUDA apps may be using the GPU; see 07-cuda-validation.md" >&2
fi

mkdir -p "$WORKDIR/cpu" "$WORKDIR/gpu"
cp "$CASE" "$WORKDIR/cpu/" "$WORKDIR/gpu/"
CASEDIR="$(dirname "$CASE")"
[[ -f "$CASEDIR/gauss.exc" ]] && cp "$CASEDIR/gauss.exc" "$WORKDIR/cpu/" "$WORKDIR/gpu/" || true
BN="$(basename "$CASE")"

echo "CPU golden -> $WORKDIR/cpu"
( cd "$WORKDIR/cpu" && OMP_NUM_THREADS=1 env -u SEMBA_FDTD_DEVICE "$EXE" -i "$BN" >cpu.log )
echo "CUDA run  -> $WORKDIR/gpu"
( cd "$WORKDIR/gpu" && OMP_NUM_THREADS=1 SEMBA_FDTD_DEVICE=cuda "$EXE" -i "$BN" >gpu.log )

python3 - "$WORKDIR" <<'PY'
import sys, numpy as np
from pathlib import Path
wd = Path(sys.argv[1])
cpu, gpu = wd/'cpu', wd/'gpu'
ok = True
rtol = 1e-3
for f in sorted(cpu.glob('*_tm.dat')):
    g = gpu/f.name
    def load(p):
        rows=[]
        for line in p.read_text().splitlines():
            line=line.strip()
            if not line or line[0] in '#%' or line[0].isalpha(): continue
            try: rows.append([float(x) for x in line.split()])
            except ValueError: continue
        return np.asarray(rows, float)
    a, b = load(f), load(g)
    af, bf = a[:, 1:], b[:, 1:]
    scale = float(np.max(np.abs(af))) if af.size else 0.0
    maxdiff = float(np.max(np.abs(af - bf))) if af.size else 0.0
    reldiff = maxdiff / scale if scale > 0 else 0.0
    close = scale > 1e-12 and reldiff <= rtol
    print(f'{f.name}: scale={scale:.3e} maxdiff={maxdiff:.3e} reldiff={reldiff:.3e} close={close}')
    ok &= close
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
PY
