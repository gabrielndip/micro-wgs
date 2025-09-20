# **Snakemake Pipeline for Bacterial Genome Variant Analysis**

Note: This document has been superseded by the root `README.md` and `workflow_patched.md` content merged there. It may be outdated; prefer `README.md` for the authoritative instructions.

This Snakemake workflow performs quality control (QC) on raw Illumina reads, read trimming, genome assembly, and MLST typing in a minimal, easy-to-run pipeline. It uses well-known tools: FastQC for read quality checks
bioconda.github.io, fastp for trimming and filtering reads github.com, SPAdes for de novo genome assembly bioconda.github.io, and Torsten Seemann’s mlst tool to determine the multilocus sequence type from the assembled contigs snakemake-wrappers.readthedocs.io. Each major step runs in its own conda environment for reproducibility (fastp is assumed to be pre-installed on the system). The pipeline expects paired-end FASTQ files (e.g. SRR15334628_1.fastq and SRR15334628_1.fastq) and organizes outputs into logical folders as shown below.

# **Project Folder Structure**

Below is a recommended project structure. Raw data goes in the `data/ folder`, and results are grouped by step (QC reports, trimmed reads, assembly output, MLST results). The conda environment YAMLs are stored in an `envs/` folder, and a `config.yaml file holds sample-specific settings.
proj01/  
├── Snakefile                   # Snakemake pipeline  
├── config.yaml                 # Configuration (sample name, input file paths, etc.)  
├── envs/                       # Conda environments for each tool  
│   ├── fastqc_env.yaml         # FastQC environment  
│   ├── spades_env.yaml         # SPAdes environment  
│   └── mlst_env.yaml           # mlst tool environment  
├── data/                       # Input data folder  
│   └── SRR15334628_1.fastq      # Raw reads (forward)  
│   └── SRR15334628_2.fastq      # Raw reads (reverse)  
└── results/                    # Output results folder  
    ├── qc/                     # FastQC reports (and fastp report)  
    │   ├── sample_R1_fastqc.html  
    │   ├── sample_R1_fastqc.zip  
    │   ├── sample_R2_fastqc.html  
    │   ├── sample_R2_fastqc.zip  
    │   ├── sample_fastp.html  
    │   └── sample_fastp.json  
    ├── trimmed/                # Trimmed reads from fastp  
    │   ├── sample_R1_trimmed.fastq.gz  
    │   └── sample_R2_trimmed.fastq.gz  
    ├── assembly/               # SPAdes assembly output  
    │   └── sample/             # (SPAdes creates a folder per sample)  
    │       └── contigs.fasta   # Assembled contigs for the sample  
    └── mlst/                   # MLST typing results  
        └── sample_mlst.tsv     # MLST output table for the sample  

# **Configuration (config.yaml)**

The config.yaml file defines user-editable parameters such as the sample name and input file paths. You can add multiple samples by expanding this configuration (e.g., using a list or additional entries), but here we demonstrate a single sample for simplicity. Update these paths and sample name as needed for my dataset.

## **config.yaml - pipeline configuration**
sample: "sample"                         # Sample name (used for naming outputs)
fastq_r1: "data/SRR15334628_1.fastq"      # Path to raw forward reads (FASTQ)
fastq_r2: "data/SRR15334628_2.fastq"      # Path to raw reverse reads (FASTQ)

# **Snakemake Workflow (Snakefile)**

Below is the Snakefile defining the workflow. It includes four rules corresponding to the steps: FastQC, fastp trimming, SPAdes assembly, and MLST typing. The rule all aggregates final outputs so that running snakemake `--use-conda --cores 4` will execute all steps to produce the QC reports, assembly, and MLST result. Each rule has a short description and uses conda environments (except fastp, which is assumed to be installed on the system PATH). Comments are provided throughout for clarity.

## **Snakefile: Bacterial variant analysis pipeline (QC, trimming, assembly, MLST)**

### **Load sample information from the config file**
SAMPLE   = config["sample"]
FASTQ_R1 = config["fastq_r1"]
FASTQ_R2 = config["fastq_r2"]

# Default target: final QC reports, assembly contigs, and MLST typing result
rule all:
    input:
        # FastQC HTML reports for forward and reverse reads
        f"results/qc/{SAMPLE}_R1_fastqc.html",
        f"results/qc/{SAMPLE}_R2_fastqc.html",
        # SPAdes assembled contigs and MLST result table
        f"results/assembly/{SAMPLE}/contigs.fasta",
        f"results/mlst/{SAMPLE}_mlst.tsv"

#################################################################
# 1. Quality Control with FastQC
#################################################################
rule fastqc_raw_reads:
    """Run FastQC on raw reads to assess base quality and other QC metrics"""
    input:
        R1=FASTQ_R1,
        R2=FASTQ_R2	
    output:
        # FastQC generates an HTML report and a ZIP archive per input FASTQ
        R1_html = f"results/qc/{SAMPLE}_R1_fastqc.html",
        R1_zip  = f"results/qc/{SAMPLE}_R1_fastqc.zip",
        R2_html = f"results/qc/{SAMPLE}_R2_fastqc.html",
        R2_zip  = f"results/qc/{SAMPLE}_R2_fastqc.zip"
    threads: 1
    conda: "envs/fastqc_env.yaml"
    shell:
        # The --outdir option directs FastQC output into the results/qc folder
        "fastqc {input.R1} {input.R2} --outdir results/qc"

#################################################################
# 2. Read Trimming with fastp
#################################################################
rule trim_reads_fastp:
    """Trim adapters and low-quality bases using fastp (paired-end)"""
    input:
        R1 = FASTQ_R1,
        R2 = FASTQ_R2
    output:
        # Trimmed FASTQ outputs (gzipped) and fastp's HTML/JSON reports
        trimmed_R1 = f"results/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        trimmed_R2 = f"results/trimmed/{SAMPLE}_R2_trimmed.fastq.gz",
        html_report = f"results/qc/{SAMPLE}_fastp.html",
        json_report = f"results/qc/{SAMPLE}_fastp.json"
    threads: 2
    # (No conda environment here; assumes system-installed fastp)
    shell:
        # fastp command with input (-i/-I) and output (-o/-O) for paired reads:contentRef
        # -h and -j options save the HTML/JSON reports to specified files
        "fastp -i {input.R1} -I {input.R2} "
        "-o {output.trimmed_R1} -O {output.trimmed_R2} "
        "-h {output.html_report} -j {output.json_report}"

#################################################################
# 3. Genome Assembly with SPAdes
#################################################################
rule assemble_genome_spades:
    """Assemble the trimmed reads into contigs using SPAdes"""
    input:
        R1 = f"results/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        R2 = f"results/trimmed/{SAMPLE}_R2_trimmed.fastq.gz"
    output:
        # The main assembly output is the contigs FASTA file
        contigs = f"results/assembly/{SAMPLE}/contigs.fasta"
    threads: 4
    conda: "envs/spades_env.yaml"
    params:
        outdir = f"results/assembly/{SAMPLE}"   # SPAdes output directory
    shell:
        # Run SPAdes with 4 threads (-t) on paired reads, outputting to the sample's assembly folder
        "spades.py -1 {input.R1} -2 {input.R2} -o {params.outdir} -t {threads}"

#################################################################
# 4. MLST Typing with mlst (PubMLST databases)
#################################################################
rule mlst_typing:
    """Identify MLST type by scanning assembled contigs against PubMLST schemes"""
    input:
        assembly = f"results/assembly/{SAMPLE}/contigs.fasta"
    output:
        f"results/mlst/{SAMPLE}_mlst.tsv"
    threads: 1
    conda: "envs/mlst_env.yaml"
    shell:
        # Run mlst on the assembly FASTA; output goes to TSV file:contentReference[oaicite:5]{index=5}
        "mlst {input.assembly} > {output}"
---
Note: The first time you run the mlst tool, it may download the PubMLST databases automatically. Subsequent runs will use the cached database. Also, ensure that the `data/` folder contains my input FASTQ files named as specified in config.yaml. After setting up the files, run the pipeline with `snakemake --use-conda --cores 4` (adjust --cores as needed).

# **Conda Environment YAMLs**

Below are the conda environment files for FastQC, SPAdes, and mlst. These files specify package dependencies and channels. Using `--use-conda` will automatically create and use these isolated environments for each rule. (No environment file is provided for fastp since we assume it’s already installed on my system.)

`envs/fastqc_env.yaml`
## **FastQC conda environment**
name: fastqc
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc=0.12.1   # FastQC tool for read quality control:contentReference[oaicite:6]{index=6}

`envs/spades_env.yaml`
## **SPAdes conda environment**
name: spades
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - spades=3.15.5   # SPAdes genome assembler:contentReference[oaicite:7]{index=7}

`envs/mlst_env.yaml`
## **mlst conda environment**
name: mlst
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - mlst=2.23.0     # Torsten Seemann's mlst tool for MLST typing:contentReference[oaicite:8]{index=8}

With these files in place, the pipeline is ready to run. It keeps each step simple and modular, and it produces outputs in organized directories for easy inspection. Running the command below will execute all steps using 4 CPU cores and create the results described above:

`snakemake --use-conda --cores 4
