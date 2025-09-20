import os


configfile: "config.yaml"

shell.executable("/bin/bash")

SAMPLE   = config["sample"]
FASTQ_R1 = config["fastq_r1"]
FASTQ_R2 = config["fastq_r2"]

# Get clean base names (e.g., SRR15334628_1)
R1_BASENAME = os.path.basename(FASTQ_R1).replace(".fastq.gz", "").replace(".fastq", "")
R2_BASENAME = os.path.basename(FASTQ_R2).replace(".fastq.gz", "").replace(".fastq", "")

rule all:
    input:
        f"results/qc/fastqc/{R1_BASENAME}_fastqc.html",
        f"results/qc/fastqc/{R2_BASENAME}_fastqc.html",
        f"results/qc/fastp/{SAMPLE}_fastp.html",
        f"results/qc/fastp/{SAMPLE}_fastp.json",
        f"results/assembly/{SAMPLE}/contigs.fasta",
        f"results/mlst/{SAMPLE}_mlst.tsv",
        "data/refgenome/ecoli_INF32-16-A_ref.fasta"

rule fastqc_raw_reads:
    input:
        R1=FASTQ_R1,
        R2=FASTQ_R2
    output:
        f"results/qc/fastqc/{R1_BASENAME}_fastqc.html",
        f"results/qc/fastqc/{R1_BASENAME}_fastqc.zip",
        f"results/qc/fastqc/{R2_BASENAME}_fastqc.html",
        f"results/qc/fastqc/{R2_BASENAME}_fastqc.zip"

    log:
        f"results/logs/fastqc__{SAMPLE}.log"
    threads: 1
    conda: "envs/fastqc_env.yaml"
    shell:
        (
            "mkdir -p results/qc/fastqc results/logs && "
            "fastqc {input.R1} {input.R2} --outdir results/qc/fastqc > {log} 2>&1"
        )

rule trim_reads_fastp:
    input:
        R1=FASTQ_R1,
        R2=FASTQ_R2
    output:
        trimmed_R1=f"results/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        trimmed_R2=f"results/trimmed/{SAMPLE}_R2_trimmed.fastq.gz",
        html_report=f"results/qc/fastp/{SAMPLE}_fastp.html",
        json_report=f"results/qc/fastp/{SAMPLE}_fastp.json"
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

rule assemble_genome_spades:
    input:
        R1=f"results/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        R2=f"results/trimmed/{SAMPLE}_R2_trimmed.fastq.gz"
    output:
        contigs=f"results/assembly/{SAMPLE}/contigs.fasta"
    params:
        outdir=f"results/assembly/{SAMPLE}"
    log:
        f"results/logs/spades__{SAMPLE}.log"
    threads: 4
    conda: "envs/spades_env.yaml"
    shell:
        (
            "mkdir -p {params.outdir} results/logs && "
            "spades.py -1 {input.R1} -2 {input.R2} -o {params.outdir} -t {threads} > {log} 2>&1"
        )

rule mlst_typing:
    input:
        assembly=f"results/assembly/{SAMPLE}/contigs.fasta"
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


rule download_reference:
    """Download and unpack the reference genome for E. coli INF32/16/A"""
    output:
        "data/refgenome/ecoli_INF32-16-A_ref.fasta"
    log:
        "results/logs/download_refgenome.log"
    conda:
        "envs/wget_env.yaml"
    shell:
        """
        mkdir -p data/refgenome results/logs
        wget -O data/refgenome/tmp_ref.fna.gz \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/634/245/GCF_017634245.1_ASM1763424v1/GCF_017634245.1_ASM1763424v1_genomic.fna.gz \
            >> {log} 2>&1
        gunzip -c data/refgenome/tmp_ref.fna.gz > {output} 2>> {log}
        rm data/refgenome/tmp_ref.fna.gz
        """
