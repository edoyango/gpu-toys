#!/usr/bin/env bash
#
# Build the three NVHPC directive-form repro_tests binaries used by
# run_nvhpc_directive_comparison.sh.
#
# Usage:
#   ./build_nvhpc_directive_bins.sh [log-root]
#
# Optional environment overrides:
#   NVFORTRAN=/path/to/nvfortran
#   BASE_FLAGS="..."
#   EXTRA_FFLAGS="..."
#   MAKE=make

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

NVFORTRAN="${NVFORTRAN:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/bin/nvfortran}"
BASE_FLAGS="${BASE_FLAGS:--mp=gpu -acc=gpu -gpu=mem:separate -O4 -stdpar=gpu -Minline=name:flux_elem -Mnovect -Mnofma}"
EXTRA_FFLAGS="${EXTRA_FFLAGS:-}"
MAKE="${MAKE:-make}"
LOG_ROOT="${1:-nvhpc_directive_build_logs}"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${LOG_ROOT%/}/${STAMP}"

mkdir -p "$LOG_DIR"

if [[ ! -x "$NVFORTRAN" ]]; then
  echo "nvfortran not found or not executable: $NVFORTRAN" >&2
  exit 127
fi

{
  echo "timestamp=${STAMP}"
  echo "workdir=${SCRIPT_DIR}"
  echo "log_dir=${LOG_DIR}"
  echo "nvfortran=${NVFORTRAN}"
  echo "base_flags=${BASE_FLAGS}"
  echo "extra_fflags=${EXTRA_FFLAGS}"
  echo "make=${MAKE}"
  echo
  echo "nvfortran --version:"
  "$NVFORTRAN" --version 2>&1 || true
  echo
  echo "git:"
  git rev-parse HEAD 2>&1 || true
  git status --short 2>&1 || true
} > "$LOG_DIR/manifest.txt"

build_variant() {
  local variant="$1"
  local macro_flags="$2"
  local output_bin="$3"
  local log_file="$LOG_DIR/${variant}.build.log"
  local flags

  flags="${macro_flags} ${BASE_FLAGS} ${EXTRA_FFLAGS}"
  # Collapse repeated whitespace for a cleaner manifest/log command line.
  flags="$(printf '%s\n' "$flags" | tr '\n' ' ' | sed 's/[[:space:]][[:space:]]*/ /g; s/^ //; s/ $//')"

  {
    echo "variant=${variant}"
    echo "output=${output_bin}"
    echo "command=${MAKE} -B repro_tests FC=${NVFORTRAN} FFLAGS=\"${flags}\""
    echo
  } > "$log_file"

  echo "Building ${variant} -> ${output_bin}"
  rm -f repro_tests
  if "$MAKE" -B repro_tests FC="$NVFORTRAN" FFLAGS="$flags" >> "$log_file" 2>&1; then
    cp repro_tests "$output_bin"
    chmod +x "$output_bin"
    ls -lh "$output_bin" >> "$log_file"
    echo "Built ${output_bin}"
  else
    echo "Build failed for ${variant}; see ${log_file}" >&2
    tail -n 80 "$log_file" >&2 || true
    exit 1
  fi
}

build_variant "omp_long" ""       "repro_tests_nv_omp_long"
build_variant "omp_loop" "-DLOOP" "repro_tests_nv_omp_loop"
build_variant "acc"      "-DACC"  "repro_tests_nv_acc"

{
  echo "Built binaries:"
  ls -lh repro_tests_nv_omp_long repro_tests_nv_omp_loop repro_tests_nv_acc
  echo
  echo "Build logs:"
  ls -lh "$LOG_DIR"
} | tee "$LOG_DIR/summary.txt"

echo "Done. Build logs are in:"
echo "$LOG_DIR"
