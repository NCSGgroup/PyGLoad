#!/usr/bin/env bash
# PyGLoad environment setup (same scheme as PyGLDA v2: the code is not installed, only its environment).
#
#   bash installation/setup_env.sh check                       # what would happen; installs nothing
#   bash installation/setup_env.sh install                     # create the env, verify it, write the lock file
#   bash installation/setup_env.sh install --data-root DIR     # ... and register the data root
#
# Options:  --name NAME      environment name (default: pygload, from environment.yml)
#           --data-root DIR  write installation/pygload_paths.cfg with data_root = DIR (must lie outside the code)
#           --no-pygmt       leave pygmt out (it is only needed for figures)
#
# The code itself needs no installation: scripts put the repository folder on the path and import `src`
# (see doc/installation.md). After `install`:  conda activate <name>  and  python -m pytest -q tests
set -o pipefail

PHASE=${1:-check}; shift || true
ENV_NAME=pygload; DATA_ROOT=""; WITH_PYGMT=1
while [ $# -gt 0 ]; do
  case "$1" in
    --name) ENV_NAME=$2; shift 2 ;;
    --data-root) DATA_ROOT=$2; shift 2 ;;
    --no-pygmt) WITH_PYGMT=0; shift ;;
    *) echo "unknown argument: $1"; exit 2 ;;
  esac
done
case "$PHASE" in check|install) ;; *) echo "usage: bash installation/setup_env.sh check|install [options]"; exit 2 ;; esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # installation/
ROOT="$(dirname "$HERE")"                                # repository folder
YML="$HERE/environment.yml"; LOCK="$HERE/environment.lock.yml"
section () { echo; echo "----- $1 -----"; }
die ()     { echo "ERROR: $*"; exit 1; }

echo "=== PyGLoad env setup | phase=$PHASE | env=$ENV_NAME | pygmt=$WITH_PYGMT | $(date) ==="
echo "repository: $ROOT"

section "solver"
if command -v mamba >/dev/null 2>&1; then SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then SOLVER=conda
else die "no conda/mamba on PATH: install Miniforge (https://github.com/conda-forge/miniforge) and rerun"; fi
echo "using: $SOLVER ($($SOLVER --version 2>&1 | head -1))"

YML_USED="$YML"
if [ "$WITH_PYGMT" = 0 ] || [ "$ENV_NAME" != pygload ]; then
  YML_USED="$(mktemp --suffix=.yml)"
  sed -e "s/^name: .*/name: $ENV_NAME/" "$YML" | { [ "$WITH_PYGMT" = 0 ] && grep -v '^\s*- pygmt' || cat; } > "$YML_USED"
  echo "(adapted copy of environment.yml: $YML_USED)"
fi

ENV_DIR="$($SOLVER env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n{print $NF}')"

if [ -n "$DATA_ROOT" ]; then
  DATA_ROOT="$(cd "$DATA_ROOT" 2>/dev/null && pwd)" || die "data root $DATA_ROOT does not exist"
  case "$DATA_ROOT/" in "$ROOT"/*) die "the data root must lie outside the code folder ($ROOT)" ;; esac
fi

if [ "$PHASE" = check ]; then
  section "check"
  [ -n "$ENV_DIR" ] && echo "env '$ENV_NAME' exists at $ENV_DIR (install would only verify it)" \
                    || { echo "env '$ENV_NAME' does not exist; dry-run of its creation:"; $SOLVER env create -f "$YML_USED" --dry-run 2>&1 | tail -15; }
  [ -n "$DATA_ROOT" ] && echo "would write $HERE/pygload_paths.cfg with data_root = $DATA_ROOT"
  [ -f "$HERE/pygload_paths.cfg" ] && echo "current data root: $(grep -E '^\s*data_root' "$HERE/pygload_paths.cfg")" \
                                    || echo "no installation/pygload_paths.cfg yet (data root not set)"
  echo; echo "=== check done. Nothing was installed. Rerun with 'install' to proceed. ==="; exit 0
fi

section "create"
if [ -n "$ENV_DIR" ]; then
  echo "env '$ENV_NAME' already exists at $ENV_DIR: not recreated (remove it with '$SOLVER env remove -n $ENV_NAME' to rebuild)"
else
  $SOLVER env create -f "$YML_USED" -y 2>&1 | tail -5 || die "env creation failed"
  ENV_DIR="$($SOLVER env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n{print $NF}')"
fi
ENV_PY="$ENV_DIR/bin/python"; [ -x "$ENV_PY" ] || die "cannot find the env's python under $ENV_DIR"
echo "env dir: $ENV_DIR"

if [ -n "$DATA_ROOT" ]; then
  section "data root"
  sed "s|^data_root = .*|data_root = $DATA_ROOT|" "$HERE/pygload_paths.cfg.example" > "$HERE/pygload_paths.cfg"
  echo "wrote $HERE/pygload_paths.cfg (data_root = $DATA_ROOT)"
fi

section "verify: packages and the PyGLoad code"
(cd /tmp && "$ENV_PY" - "$ROOT" <<'PY') || die "verification failed"
import importlib, sys
print("python", sys.version.split()[0], sys.executable)
bad = []
for p in ["numpy", "scipy", "h5py", "netCDF4", "threadpoolctl", "geopandas", "shapely", "pytest", "pygmt"]:
    try:
        m = importlib.import_module(p); print(f"  OK   {p:14s} {getattr(m, '__version__', '?')}")
    except Exception as e:
        print(f"  {'SKIP' if p == 'pygmt' else 'FAIL'} {p:14s} {type(e).__name__}: {e}")
        if p != "pygmt": bad.append(p)
sys.path.insert(0, sys.argv[1])          # what every PyGLoad script does: the repository folder on the path
import src as pgw
print(f"  OK   PyGLoad code  {pgw.__version__} ({pgw.__file__})")
try:
    print("  data root:", pgw.FileTool.get_data_dir())
except Exception as e:
    print("  data root: not set yet:", e)
sys.exit(1 if bad else 0)
PY

section "tests"
(cd "$ROOT" && "$ENV_PY" -m pytest -q tests 2>&1 | tail -3) || die "tests failed"

section "lock"
$SOLVER env export -n "$ENV_NAME" --no-builds > "$LOCK" 2>/dev/null && echo "wrote $LOCK (exact versions of this install)"

echo; echo "=== install done | $(date) ==="
echo "next:  conda activate $ENV_NAME"
[ -f "$HERE/pygload_paths.cfg" ] || echo "       set the data root: bash installation/setup_env.sh install --data-root /path/to/PyGLoad_external_data"
