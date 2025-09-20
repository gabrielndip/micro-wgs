#!/usr/bin/env bash
set -euo pipefail

# Minimal dry-run to verify rule graph and paths parse
snakemake -n --use-conda --cores 2 \
  results/qc/fastqc/${SAMPLE:-sample}_R1_fastqc.html \
  results/qc/fastqc/${SAMPLE:-sample}_R2_fastqc.html \
  results/qc/fastp/${SAMPLE:-sample}_fastp.html \
  results/assembly/${SAMPLE:-sample}/contigs.fasta \
  results/mlst/${SAMPLE:-sample}_mlst.tsv

# If variant calling is enabled in config.yaml, include variant targets in dry-run
if grep -Eq '^\s*enable_variant_calling:\s*true\s*$' config.yaml; then
  snakemake -n --use-conda --cores 2 \
    results/variants/${SAMPLE:-sample}.filtered.vcf.gz \
    results/variants/${SAMPLE:-sample}.filtered.vcf.gz.tbi
fi
