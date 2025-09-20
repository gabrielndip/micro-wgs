#!/usr/bin/env bash
set -euo pipefail

# End-to-end tiny run using synthetic FASTQs and toy reference.
# This exercises trimming, alignment, variant calling, and reports on a tiny dataset.

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

# Create tiny data and a temporary config
bash "${ROOT_DIR}/scripts/tests/make_tiny_data.sh"
CFG="${ROOT_DIR}/.tmp/tiny/config.yaml"

# Targets limited to the variant-calling path to keep runtime short
TARGETS=(
  results/variants/tiny.filtered.vcf.gz
  results/variants/tiny.filtered.vcf.gz.tbi
  results/reports/coverage.txt
  results/reports/variants_summary.tsv
  results/reports/vcf_stats.txt
)

echo "Running Snakemake tiny E2E on targets:"
printf ' - %s\n' "${TARGETS[@]}"

snakemake --use-conda --cores 2 --configfile "${CFG}" "${TARGETS[@]}"

echo "OK: tiny end-to-end run completed. Outputs under results/."

