#!/usr/bin/env bash
set -euo pipefail

# Verify required tools are available and print versions.
# Usage: ./scripts/verify_env.sh [conda_env_name]
# Defaults to snake.

ENV_NAME="${1:-snake}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH" >&2
  exit 2
fi

TOOLS=(python snakemake fastqc fastp spades.py mlst bwa bowtie2 samtools bcftools bedtools mash multiqc shellcheck)
missing=()

echo "Checking tools in env: ${ENV_NAME}" >&2
ENV_PREFIX="${CONDA_PREFIX_OVERRIDE:-$HOME/miniconda3/envs/${ENV_NAME}}"

check_tool() {
  local t="$1"
  # Try conda run
  if CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" bash -lc "command -v ${t}" >/dev/null 2>&1; then
    ver=$(CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" bash -lc "{ ${t} --version 2>/dev/null || ${t} -V 2>/dev/null || ${t} -v 2>/dev/null || echo version: unknown; } | head -n1")
    printf "%-14s %s\n" "${t}:" "${ver}"
    return 0
  fi
  # Fallback: direct bin path
  if [ -x "${ENV_PREFIX}/bin/${t}" ]; then
    ver=$("${ENV_PREFIX}/bin/${t}" --version 2>/dev/null || "${ENV_PREFIX}/bin/${t}" -V 2>/dev/null || "${ENV_PREFIX}/bin/${t}" -v 2>/dev/null || echo "version: unknown")
    ver=$(printf "%s" "$ver" | head -n1)
    printf "%-14s %s\n" "${t}:" "${ver}"
    return 0
  fi
  # Not found
  printf "%-14s %s\n" "${t}:" "MISSING"
  missing+=("${t}")
  return 1
}

for t in "${TOOLS[@]}"; do
  check_tool "$t"
done

if [ ${#missing[@]} -gt 0 ]; then
  printf "\nMissing tools: %s\n" "${missing[*]}" >&2
  printf "Install with: conda install -n %s -c conda-forge -c bioconda %s\n" "${ENV_NAME}" "${missing[*]}" >&2
  exit 1
fi

printf "\nAll required tools are available in %s.\n" "${ENV_NAME}"
