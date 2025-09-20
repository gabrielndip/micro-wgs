# Project 1: Pathogen Genome Variant Analysis Pipeline

## Objective
Develop a Snakemake pipeline to analyze a bacterial whole-genome sequencing (WGS) sample for species identification and variant calling, simulating a clinical diagnostics scenario.

This project focuses on microbial genomics by analyzing a single pathogen’s genome from raw sequencing reads to detected variants. Whole-genome sequencing is increasingly used in clinical microbiology to identify pathogens and genetic markers (e.g. antibiotic resistance mutations). The pipeline will map sequencing reads to a reference genome and call genomic variants, emulating how clinical labs might detect strain differences or resistance-associated mutations.

## Dataset & Tools
Dataset: Illumina short-read data from a public bacterial isolate. For example, a Klebsiella pneumoniae genome run from NCBI SRA or an E. coli Illumina sample. A reference genome (FASTA) for the species will be obtained from NCBI RefSeq (using a small subset of reads if needed to fit runtime).

Tools: Snakemake for workflow management; FastQC for raw read QC; Trimmomatic or fastp for adapter trimming; BWA for read alignment; Samtools for BAM processing; FreeBayes or bcftools for variant calling; (Optional: bedtools/BLAST to check if variants hit known resistance genes). All tools will be managed via conda environments for reproducibility.

## Methodology
- Data QC: Perform quality control on raw reads (FastQC report) and trimming of low-quality bases/adapters.
- Alignment: Align filtered reads to the reference genome using BWA, producing a sorted BAM file.
- Variant Calling: Call variants (SNPs/indels) from the BAM against the reference using FreeBayes or bcftools, generating a VCF file of differences.
- Analysis: Parse the VCF to identify notable variants. For example, highlight any variant in genes related to drug resistance or virulence (if known).
- Results Summary: Produce a brief report (e.g., a table of high-confidence variants with annotations). Possibly visualize coverage or variant effects if time permits.
- Reflection: Include a short discussion on the pipeline’s relevance to diagnostics (e.g., how WGS can rapidly identify species and resistance markers) and note any limitations (speed, coverage requirements, etc.).

## Project Repository Structure
Each project will use a structured Git repository to ensure clarity and reproducibility. For Project 1:

```
pathogen-genome-variant-analysis/
├── data/          # Input data (sample reads FASTQ, reference genome FASTA)
├── workflow/      # Snakemake workflow files and scripts
│   ├── Snakefile  # Main Snakemake pipeline definition
│   └── rules/     # (Optional) Additional Snakemake rule modules
├── envs/          # Conda environment YAMLs for tools (e.g., bwa.yaml, qc.yaml)
├── results/       # Output files (QC reports, BAM, VCF, summary tables)
└── README.md      # Documentation (overview, methods, usage, results, reflection)
```

Note: The Snakefile will define rules for each step (QC, alignment, variant calling, etc.), linking inputs to outputs. Using Snakemake’s integrated package management (`--use-conda`), each rule’s environment (software versions) is captured for reproducibility. The repository will include instructions to obtain or simulate the dataset (to avoid large files on GitHub).

## Expected Deliverables
- Fully documented GitHub repository with the above structure, containing all code (Snakefile, scripts) and a comprehensive README.md.
- Snakemake workflow that automates the WGS analysis end-to-end (from FASTQ to VCF).
- Example output files: trimmed reads, alignment BAM, variant VCF, and a short report or table of key variants.
- Project report in the README: including an overview of the approach, a description of the dataset, methods (with code snippets), summary of results (e.g. “X variants were found; one non-synonymous SNP in gene gyrA associated with fluoroquinolone resistance was detected”), and a brief reflection on the pipeline’s diagnostic utility.
- The README will also detail how to run the pipeline (e.g. `snakemake --use-conda --cores 4`) and how it ensures reproducibility (mentioning fixed versions in conda environments and references to standards).

