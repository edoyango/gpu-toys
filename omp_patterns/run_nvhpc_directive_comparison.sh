#!/usr/bin/env bash
#
# Run the three NVHPC directive-form repro_tests binaries and collect enough
# raw output for later README analysis.
#
# Usage:
#   ./run_nvhpc_directive_comparison.sh [output-root]
#
# Optional environment overrides:
#   NVFORTRAN=/path/to/nvfortran
#   CUOBJDUMP=/path/to/cuobjdump
#   SM_ARCH=sm_80              # or leave unset to infer from nvidia-smi
#   TIMEOUT_SECONDS=3600       # optional per-binary timeout, if timeout exists

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

NVFORTRAN="${NVFORTRAN:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/bin/nvfortran}"
CUOBJDUMP="${CUOBJDUMP:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/cuda/13.1/bin/cuobjdump}"
OUT_ROOT="${1:-nvhpc_directive_results}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_ROOT%/}/${STAMP}"

mkdir -p "$OUT_DIR" || exit 1

variants=(
  "omp_long:repro_tests_nv_omp_long:OpenMP long form"
  "omp_loop:repro_tests_nv_omp_loop:OpenMP loop"
  "acc:repro_tests_nv_acc:OpenACC"
)

detect_sm_arch() {
  if command -v nvidia-smi >/dev/null 2>&1; then
    local cap
    if cap="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n 1)"; then
      cap="$(printf '%s\n' "$cap" | tr -d '[:space:].')"
      if [[ "$cap" =~ ^[0-9]+$ ]]; then
        printf 'sm_%s\n' "$cap"
        return 0
      fi
    fi
  fi
  printf 'sm_80\n'
}

SM_ARCH="${SM_ARCH:-$(detect_sm_arch)}"

{
  echo "timestamp=${STAMP}"
  echo "workdir=${SCRIPT_DIR}"
  echo "output_dir=${OUT_DIR}"
  echo "nvfortran=${NVFORTRAN}"
  echo "cuobjdump=${CUOBJDUMP}"
  echo "selected_sm_arch=${SM_ARCH}"
  echo
  echo "uname:"
  uname -a 2>&1 || true
  echo
  echo "hostname:"
  hostname 2>&1 || true
  echo
  echo "git:"
  git rev-parse HEAD 2>&1 || true
  git status --short 2>&1 || true
} > "$OUT_DIR/manifest.txt"

{
  echo "nvfortran --version"
  "$NVFORTRAN" --version 2>&1 || true
  echo
  echo "cuobjdump --version"
  "$CUOBJDUMP" --version 2>&1 || true
  echo
  echo "nvidia-smi"
  nvidia-smi 2>&1 || true
  echo
  echo "nvidia-smi --query-gpu"
  nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv 2>&1 || true
  echo
  echo "nvaccelinfo"
  if [[ -x "$(dirname "$NVFORTRAN")/nvaccelinfo" ]]; then
    "$(dirname "$NVFORTRAN")/nvaccelinfo" 2>&1 || true
  elif command -v nvaccelinfo >/dev/null 2>&1; then
    nvaccelinfo 2>&1 || true
  else
    echo "nvaccelinfo not found"
  fi
} > "$OUT_DIR/system_info.txt"

echo "Writing results under: $OUT_DIR"
echo "Selected resource-usage arch: $SM_ARCH"

for entry in "${variants[@]}"; do
  IFS=: read -r variant bin label <<< "$entry"

  {
    echo "variant=${variant}"
    echo "label=${label}"
    echo "binary=${bin}"
    if [[ -e "$bin" ]]; then
      ls -lh "$bin"
      echo
      echo "ldd:"
      ldd "$bin" 2>&1 || true
    else
      echo "missing binary"
    fi
  } > "$OUT_DIR/${variant}_binary_info.txt"

  if [[ -x "$CUOBJDUMP" && -x "$bin" ]]; then
    echo "Collecting cuobjdump resource usage for ${label}"
    "$CUOBJDUMP" -res-usage "$bin" > "$OUT_DIR/${variant}_cuobjdump_res_usage.txt" 2>&1
    echo "$?" > "$OUT_DIR/${variant}_cuobjdump.status"
  else
    echo "Skipping cuobjdump for ${label}; missing tool or binary" | tee "$OUT_DIR/${variant}_cuobjdump_res_usage.txt"
    echo "127" > "$OUT_DIR/${variant}_cuobjdump.status"
  fi
done

for entry in "${variants[@]}"; do
  IFS=: read -r variant bin label <<< "$entry"

  if [[ ! -x "$bin" ]]; then
    echo "Skipping run for ${label}; missing executable ${bin}"
    echo "127" > "$OUT_DIR/${variant}.status"
    continue
  fi

  echo "Running ${label}: ./${bin}"
  start_epoch="$(date +%s)"
  if [[ -n "${TIMEOUT_SECONDS:-}" ]] && command -v timeout >/dev/null 2>&1; then
    timeout "$TIMEOUT_SECONDS" "./$bin" > "$OUT_DIR/${variant}.out" 2> "$OUT_DIR/${variant}.err"
    status="$?"
  else
    "./$bin" > "$OUT_DIR/${variant}.out" 2> "$OUT_DIR/${variant}.err"
    status="$?"
  fi
  end_epoch="$(date +%s)"

  {
    echo "status=${status}"
    echo "elapsed_seconds=$((end_epoch - start_epoch))"
  } > "$OUT_DIR/${variant}.status"

  echo "Finished ${label} with status ${status}"
done

if command -v python3 >/dev/null 2>&1; then
  python3 - "$OUT_DIR" "$SM_ARCH" <<'PY'
from pathlib import Path
import re
import sys

out_dir = Path(sys.argv[1])
sm_arch = sys.argv[2]

variants = [
    ("omp_long", "OpenMP long form"),
    ("omp_loop", "OpenMP loop"),
    ("acc", "OpenACC"),
]

def status_for(variant):
    path = out_dir / f"{variant}.status"
    if not path.exists():
        return ""
    text = path.read_text(errors="replace")
    m = re.search(r"status=(\d+)", text)
    return m.group(1) if m else ""

timing_rows = []
for variant, label in variants:
    text_path = out_dir / f"{variant}.out"
    if not text_path.exists():
        continue
    text = text_path.read_text(errors="replace")
    status = status_for(variant)
    current = None
    current_test = None
    for line in text.splitlines():
        m = re.search(r"Size:\s*ni=\s*(\d+)\s+nj=\s*(\d+)\s+nz=\s*(\d+)\s+nteams=\s*(\d+)", line)
        if m:
            current = tuple(m.groups())
            current_test = None
            continue
        m = re.search(r"^\s*(Test\s+\d+[^:]*):\s*$", line)
        if m:
            current_test = m.group(1).strip()
            continue
        m = re.search(
            r"min=\s*([0-9.Ee+-]+)\s+max=\s*([0-9.Ee+-]+)\s+avg=\s*([0-9.Ee+-]+)\s+med=\s*([0-9.Ee+-]+)\s+s",
            line,
        )
        if m and current and current_test:
            min_s, max_s, avg_s, med_s = (float(x) for x in m.groups())
            timing_rows.append((
                variant, label, status, *current, current_test,
                min_s, max_s, avg_s, med_s,
                min_s * 1000.0, max_s * 1000.0, avg_s * 1000.0, med_s * 1000.0,
            ))

with (out_dir / "timings.tsv").open("w") as f:
    f.write(
        "variant\tlabel\tstatus\tni\tnj\tnz\tnteams\ttest\t"
        "min_s\tmax_s\tavg_s\tmed_s\tmin_ms\tmax_ms\tavg_ms\tmed_ms\n"
    )
    for row in timing_rows:
        f.write("\t".join(str(x) for x in row) + "\n")

def arch_blocks(text):
    matches = list(re.finditer(r"Fatbin elf code:\n=+\narch = (sm_\d+)\b", text))
    for i, match in enumerate(matches):
        start = match.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        yield match.group(1), text[start:end]

resource_rows = []
for variant, label in variants:
    path = out_dir / f"{variant}_cuobjdump_res_usage.txt"
    if not path.exists():
        continue
    text = path.read_text(errors="replace")
    for arch, block in arch_blocks(text):
        for m in re.finditer(
            r"Function ([^:]+):\n\s+REG:(\d+) STACK:(\d+) SHARED:(\d+) LOCAL:(\d+)",
            block,
        ):
            resource_rows.append((
                variant, label, arch, m.group(1),
                m.group(2), m.group(3), m.group(4), m.group(5),
            ))

with (out_dir / "resources_all_arches.tsv").open("w") as f:
    f.write("variant\tlabel\tarch\tfunction\treg\tstack\tshared\tlocal\n")
    for row in resource_rows:
        f.write("\t".join(row) + "\n")

with (out_dir / f"resources_{sm_arch}.tsv").open("w") as f:
    f.write("variant\tlabel\tarch\tfunction\treg\tstack\tshared\tlocal\n")
    for row in resource_rows:
        if row[2] == sm_arch:
            f.write("\t".join(row) + "\n")
PY
  echo "Wrote parsed summaries: timings.tsv, resources_all_arches.tsv, resources_${SM_ARCH}.tsv"
else
  echo "python3 not found; raw output files were still written"
fi

{
  for file in \
    manifest.txt \
    system_info.txt \
    omp_long.status \
    omp_loop.status \
    acc.status \
    timings.tsv \
    "resources_${SM_ARCH}.tsv" \
    omp_long.out \
    omp_long.err \
    omp_loop.out \
    omp_loop.err \
    acc.out \
    acc.err; do
    echo "===== ${file} ====="
    if [[ -f "$OUT_DIR/$file" ]]; then
      cat "$OUT_DIR/$file"
    else
      echo "missing"
    fi
    echo
  done
} > "$OUT_DIR/combined_results.txt"

echo "Wrote combined summary/raw file: combined_results.txt"
echo "Done. Archive or send this directory for analysis:"
echo "$OUT_DIR"
