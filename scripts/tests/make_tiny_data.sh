#!/usr/bin/env bash
set -euo pipefail

# Create tiny paired FASTQs and a minimal config file for an end-to-end test.
# - Uses the toy reference at data/species_refs/Toy_bacterium.fna
# - Produces a temporary config at .tmp/tiny/config.yaml

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="${ROOT_DIR}/.tmp/tiny"
mkdir -p "${TMP_DIR}" "${ROOT_DIR}/data" "${ROOT_DIR}/results"

# Ensure optional assets exist (toy reference)
if [ ! -f "${ROOT_DIR}/data/species_refs/Toy_bacterium.fna" ]; then
  bash "${ROOT_DIR}/scripts/setup_optional_assets.sh"
fi

R1="${TMP_DIR}/tiny_R1.fastq"
R2="${TMP_DIR}/tiny_R2.fastq"

# Minimal 4-read FASTQs (paired). Sequences loosely match the toy reference.
cat >"${R1}" << 'EOF'
@tiny1/1
ATGCGTACGTTAGCTAGCTAGGCT

+
IIIIIIIIIIIIIIIIIIIIIIII
@tiny2/1
GATCGATCGGCTAGCTAGCTAGCT

+
IIIIIIIIIIIIIIIIIIIIIIII
EOF

cat >"${R2}" << 'EOF'
@tiny1/2
AGC TAGCTA CGATCGATCGGATC

+
IIIIIIIIIIIIIIIIIIIIIIII
@tiny2/2
TAGCTAGCTAGCTAGCTAGATCG

+
IIIIIIIIIIIIIIIIIIIIIIII
EOF

# Remove spaces in the synthetic sequence lines (ensure contiguous bases)
sed -i.bak 's/ //g' "${R1}" && rm -f "${R1}.bak"
sed -i.bak 's/ //g' "${R2}" && rm -f "${R2}.bak"

CFG="${TMP_DIR}/config.yaml"
cat >"${CFG}" << EOF
sample: "tiny"
fastq_r1: "${R1}"
fastq_r2: "${R2}"
enable_variant_calling: true
reference_fasta: "${ROOT_DIR}/data/species_refs/Toy_bacterium.fna"
reference_sha256: ""
aligner: "bwa"
mark_duplicates: false
caller: "bcftools"
ploidy: 1

# Make filters permissive for tiny depth
filter_min_qual: 0
filter_min_dp: 1

fastp_extra_args: "--cut_front --cut_tail --length_required 20"

species_id_enabled: false
species_refs_dir: "${ROOT_DIR}/data/species_refs"
resistance_annotation_enabled: false
resistance_genes_bed: "${ROOT_DIR}/data/resistance_genes.bed"
EOF

echo "Created tiny dataset and config:"
echo "- R1: ${R1}"
echo "- R2: ${R2}"
echo "- CFG: ${CFG}"

