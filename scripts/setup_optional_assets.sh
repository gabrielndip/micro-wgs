#!/usr/bin/env bash
set -euo pipefail

# Scaffold optional assets under data/ without committing large files.
# - Creates data/species_refs/ with a README and .gitkeep
# - Creates data/resistance_genes.bed.example template

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

mkdir -p "${ROOT_DIR}/data/species_refs"

# .gitkeep to allow empty dir tracking (note: .gitignore already whitelists it)
if [ ! -f "${ROOT_DIR}/data/species_refs/.gitkeep" ]; then
  : > "${ROOT_DIR}/data/species_refs/.gitkeep"
fi

# README for species refs
SPECIES_README="${ROOT_DIR}/data/species_refs/README.md"
if [ ! -f "${SPECIES_README}" ]; then
  cat > "${SPECIES_README}" << 'EOF'
# species_refs

Place a small set of reference FASTA files here (e.g., representative genomes):
- Accepted extensions: .fa, .fna, .fasta
- Keep the set minimal to speed up Mash; these files are not tracked by git.

Example filenames (you provide the actual FASTA files):
- Escherichia_coli_K12_MG1655.fna
- Klebsiella_pneumoniae_NTUH_K2044.fna

These references are used to compute Mash distances against the assembled contigs
and identify the closest species/strain.
EOF
fi

# BED template for resistance genes
BED_EXAMPLE="${ROOT_DIR}/data/resistance_genes.bed.example"
if [ ! -f "${BED_EXAMPLE}" ]; then
  cat > "${BED_EXAMPLE}" << 'EOF'
# BED template: resistance genes on the SAME reference used for alignment
# Columns (BED3+): chrom  start  end  gene  [optional:score strand]
# Coordinates are 0-based, half-open per BED convention.
#
# Example entries (replace with your reference-specific coordinates):
#chrom	start	end	gene
chr1	10000	10500	blaTEM-1
chr1	250000	251200	gyrA
EOF
fi

# Create default BED if missing (points config to this path)
BED_DEFAULT="${ROOT_DIR}/data/resistance_genes.bed"
if [ ! -f "${BED_DEFAULT}" ]; then
  cp "${BED_EXAMPLE}" "${BED_DEFAULT}"
fi

# Tiny toy FASTA for Mash species ID (kept small, untracked by git)
TOY_FASTA="${ROOT_DIR}/data/species_refs/Toy_bacterium.fna"
if [ ! -f "${TOY_FASTA}" ]; then
  cat > "${TOY_FASTA}" << 'EOF'
>Toy_bacterium
ATGCGTACGTTAGCTAGCTAGGCTAGCTAGGATCGATCGATCGGATCCGATCGATCGTAGCTAGCTAGCTAGCTA
GATCGATCGGCTAGCTAGCTAGCTAGGATCCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGATCG
ATCGGATCCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGATCGATCGATCGATCGATCGGATCCA
TCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGATCGATCGATCGGATCCGATCGATCGATCGTAG
CTAGCTAGCTAGCTAGCTAGCTAGATCGATCGATCGATCGGATCCGATCGATCGTAGCTAGCTAGCTAGCTAGCT
EOF
fi

echo "Created/verified:"
echo "- data/species_refs/.gitkeep"
echo "- data/species_refs/README.md"
echo "- data/resistance_genes.bed.example"
echo "- data/species_refs/Toy_bacterium.fna (toy reference)"
echo "- data/resistance_genes.bed (default used by config.yaml)"
