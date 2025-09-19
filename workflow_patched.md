# Snakemake Pipeline for Bacterial Genome QC, Trimming, Assembly, and MLST

This Snakemake workflow performs quality control (QC) on raw Illumina reads, read trimming, genome assembly, and MLST typing in a minimal, easy-to-run pipeline. It uses FastQC for read quality checks, fastp for trimming and filtering reads, SPAdes for de novo genome assembly, and Torsten Seemann’s mlst tool to determine the multilocus sequence type from the assembled contigs. Each major step runs in its own conda environment for reproducibility. The pipeline expects paired-end FASTQ files (e.g., SRR15334628_1.fastq and SRR15334628_2.fastq) and organizes outputs into logical folders as shown below.

## Project Folder Structure

Below is a recommended project structure. Raw data goes in the `data/` folder, and results are grouped by step (QC reports, trimmed reads, assembly output, MLST results). The conda environment YAMLs are stored in an `envs/` folder, and a `config.yaml` file holds sample-specific settings. Helper scripts and smoke tests live under `scripts/`.

```
proj01/
├── Snakefile                    # Snakemake pipeline
├── config.yaml                  # Configuration (sample name, input file paths)
├── envs/                        # Conda environments for each tool
│   ├── fastqc_env.yaml          # FastQC environment
│   ├── fastp_env.yaml           # fastp environment
│   ├── spades_env.yaml          # SPAdes environment
│   └── mlst_env.yaml            # mlst environment
├── data/                        # Input data folder
│   ├── samples.tsv              # Optional sample sheet (see example below)
│   ├── SRR15334628_1.fastq      # Raw reads (forward)
│   └── SRR15334628_2.fastq      # Raw reads (reverse)
├── scripts/
│   └── tests/                   # Smoke tests scaffolding
│       ├── test_help.sh
│       └── test_minirun.sh
└── results/                     # Output results folder
    ├── qc/
    │   ├── fastqc/              # FastQC reports
    │   │   ├── sample_R1_fastqc.html
    │   │   ├── sample_R1_fastqc.zip
    │   │   ├── sample_R2_fastqc.html
    │   │   └── sample_R2_fastqc.zip
    │   └── fastp/               # fastp reports
    │       ├── sample_fastp.html
    │       └── sample_fastp.json
    ├── trimmed/                 # Trimmed reads from fastp
    │   ├── sample_R1_trimmed.fastq.gz
    │   └── sample_R2_trimmed.fastq.gz
    ├── assembly/                # SPAdes assembly output
    │   └── sample/              # Folder per sample
    │       └── contigs.fasta    # Assembled contigs for the sample
    ├── mlst/                    # MLST typing results
    │   └── sample_mlst.tsv      # MLST output table for the sample
    ├── logs/                    # Per-rule or per-sample logs
    └── reports/                 # Summaries and provenance
```

## Configuration (config.yaml)

The `config.yaml` file defines user-editable parameters such as the sample name and input file paths. You can add multiple samples by expanding this configuration (e.g., using a list or additional entries); here we demonstrate a single sample for simplicity. Update these paths and the sample name as needed for your dataset.

```
sample: "sample"                         # Sample name (used for naming outputs)
fastq_r1: "data/SRR15334628_1.fastq"     # Path to raw forward reads (FASTQ)
fastq_r2: "data/SRR15334628_2.fastq"     # Path to raw reverse reads (FASTQ)
```

## Snakemake Workflow (Snakefile)

Below is the Snakefile defining the workflow. It includes four rules corresponding to the steps: FastQC, fastp trimming, SPAdes assembly, and MLST typing. The `rule all` aggregates final outputs so that running `snakemake --use-conda --cores 4` will execute all steps to produce the QC reports, assembly, and MLST result. Each rule uses conda environments and writes a log file under `results/logs/`.

### Load sample information from the config file

```
SAMPLE   = config["sample"]
FASTQ_R1 = config["fastq_r1"]
FASTQ_R2 = config["fastq_r2"]

# Default target: final QC reports (FastQC + fastp), assembly contigs, and MLST typing result
rule all:
    input:
        # FastQC HTML reports for forward and reverse reads
        f"results/qc/fastqc/{SAMPLE}_R1_fastqc.html",
        f"results/qc/fastqc/{SAMPLE}_R2_fastqc.html",
        # fastp HTML/JSON reports
        f"results/qc/fastp/{SAMPLE}_fastp.html",
        f"results/qc/fastp/{SAMPLE}_fastp.json",
        # SPAdes assembled contigs and MLST result table
        f"results/assembly/{SAMPLE}/contigs.fasta",
        f"results/mlst/{SAMPLE}_mlst.tsv"
```

### 1. Quality Control with FastQC

```
rule fastqc_raw_reads:
    """Run FastQC on raw reads to assess base quality and other QC metrics"""
    input:
        R1 = FASTQ_R1,
        R2 = FASTQ_R2
    output:
        R1_html = f"results/qc/fastqc/{SAMPLE}_R1_fastqc.html",
        R1_zip  = f"results/qc/fastqc/{SAMPLE}_R1_fastqc.zip",
        R2_html = f"results/qc/fastqc/{SAMPLE}_R2_fastqc.html",
        R2_zip  = f"results/qc/fastqc/{SAMPLE}_R2_fastqc.zip"
    log:
        f"results/logs/fastqc__{SAMPLE}.log"
    threads: 1
    conda: "envs/fastqc_env.yaml"
    shell:
        (
            "mkdir -p results/qc/fastqc results/logs && "
            "fastqc {input.R1} {input.R2} --outdir results/qc/fastqc > {log} 2>&1"
        )
```

### 2. Read Trimming with fastp

```
rule trim_reads_fastp:
    """Trim adapters and low-quality bases using fastp (paired-end)"""
    input:
        R1 = FASTQ_R1,
        R2 = FASTQ_R2
    output:
        trimmed_R1 = f"results/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        trimmed_R2 = f"results/trimmed/{SAMPLE}_R2_trimmed.fastq.gz",
        html_report = f"results/qc/fastp/{SAMPLE}_fastp.html",
        json_report = f"results/qc/fastp/{SAMPLE}_fastp.json"
    log:
        f"results/logs/fastp__{SAMPLE}.log"
    threads: 2
    conda: "envs/fastp_env.yaml"
    shell:
        (
            "mkdir -p results/trimmed results/qc/fastp results/logs && "
            "fastp -i {input.R1} -I {input.R2} "
            "-o {output.trimmed_R1} -O {output.trimmed_R2} "
            "-h {output.html_report} -j {output.json_report} > {log} 2>&1"
        )
```

### 3. Genome Assembly with SPAdes

```
rule assemble_genome_spades:
    """Assemble the trimmed reads into contigs using SPAdes"""
    input:
        R1 = f"results/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        R2 = f"results/trimmed/{SAMPLE}_R2_trimmed.fastq.gz"
    output:
        contigs = f"results/assembly/{SAMPLE}/contigs.fasta"
    params:
        outdir = f"results/assembly/{SAMPLE}"
    log:
        f"results/logs/spades__{SAMPLE}.log"
    threads: 4
    conda: "envs/spades_env.yaml"
    shell:
        (
            "mkdir -p {params.outdir} results/logs && "
            "spades.py -1 {input.R1} -2 {input.R2} -o {params.outdir} -t {threads} > {log} 2>&1"
        )
```

### 4. MLST Typing with mlst (PubMLST databases)

```
rule mlst_typing:
    """Identify MLST type by scanning assembled contigs against PubMLST schemes"""
    input:
        assembly = f"results/assembly/{SAMPLE}/contigs.fasta"
    output:
        f"results/mlst/{SAMPLE}_mlst.tsv"
    log:
        f"results/logs/mlst__{SAMPLE}.log"
    threads: 1
    conda: "envs/mlst_env.yaml"
    shell:
        (
            "mkdir -p results/mlst results/logs && "
            "mlst {input.assembly} 1> {output} 2> {log}"
        )
```

Note: The first time you run the mlst tool, it may download the PubMLST databases automatically. Subsequent runs will use the cached database. Ensure that the `data/` folder contains your input FASTQ files named as specified in `config.yaml`. After setting up the files, run the pipeline with `snakemake --use-conda --cores 4` (adjust `--cores` as needed).

## Conda Environment YAMLs

Below are the conda environment files for FastQC, fastp, SPAdes, and mlst. Using `--use-conda` will automatically create and use these isolated environments for each rule.

`envs/fastqc_env.yaml`

```
name: fastqc
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc=0.12.1
```

`envs/fastp_env.yaml`

```
name: fastp
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastp=0.23.4
```

`envs/spades_env.yaml`

```
name: spades
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - spades=3.15.5
```

`envs/mlst_env.yaml`

```
name: mlst
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - mlst=2.23.0
```

## Optional Sample Sheet (data/samples.tsv)

If you later expand to multiple samples, consider using a simple TSV file to manage inputs. Here is a minimal example you can place at `data/samples.tsv`:

```
sample_id	fastq1	fastq2
sample	data/SRR15334628_1.fastq	data/SRR15334628_2.fastq
```

## Smoke Tests Scaffolding (scripts/tests/)

Add the following files as quick smoke tests to validate the setup. They exercise Snakemake parsing and perform a tiny dry-run. You can adapt `test_minirun.sh` to run end-to-end on tiny fixtures.

`scripts/tests/test_help.sh`

```
#!/usr/bin/env bash
set -euo pipefail

echo "Snakemake version:" && snakemake --version
echo "Listing rules:" && snakemake -n --cores 1 --list
```

`scripts/tests/test_minirun.sh`

```
#!/usr/bin/env bash
set -euo pipefail

# Minimal dry-run to verify rule graph and paths parse
snakemake -n --use-conda --cores 2 \
  results/qc/fastqc/${SAMPLE:-sample}_R1_fastqc.html \
  results/qc/fastqc/${SAMPLE:-sample}_R2_fastqc.html \
  results/qc/fastp/${SAMPLE:-sample}_fastp.html \
  results/assembly/${SAMPLE:-sample}/contigs.fasta \
  results/mlst/${SAMPLE:-sample}_mlst.tsv

# For a real mini-run, ensure tiny FASTQs exist in data/ and run:
# snakemake --use-conda --cores 2
```

## Usage

- Create envs on demand via `--use-conda`, or pre-create with `mamba env create -f envs/<tool>_env.yaml`.
- Run: `snakemake --use-conda --cores 4`
- Outputs appear under `results/` in the subfolders outlined above.

## Notes and Corrections from Original Draft

- Title now reflects actual steps (no variant calling).
- Paired-end example fixed to `_1.fastq` and `_2.fastq`.
- QC outputs split into `results/qc/fastqc/` and `results/qc/fastp/`.
- Added `results/logs/` and `results/reports/` directories to align with repository guidelines.
- Added `envs/fastp_env.yaml` and `conda:` directive to the fastp rule.
- Added `log:` directives and ensured commands write logs via redirection.
- Removed stray artifacts and truncations present in the original text.

