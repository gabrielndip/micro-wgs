configfile: "config.yaml"

shell.executable("/bin/bash")

SAMPLE   = config["sample"]
FASTQ_R1 = config["fastq_r1"]
FASTQ_R2 = config["fastq_r2"]

rule all:
    input:
        f"results/qc/fastqc/{SAMPLE}_R1_fastqc.html",
        f"results/qc/fastqc/{SAMPLE}_R2_fastqc.html",
        f"results/qc/fastp/{SAMPLE}_fastp.html",
        f"results/qc/fastp/{SAMPLE}_fastp.json",
        f"results/assembly/{SAMPLE}/contigs.fasta",
        f"results/mlst/{SAMPLE}_mlst.tsv"

rule fastqc_raw_reads:
    input:
        R1=FASTQ_R1,
        R2=FASTQ_R2
    output:
        R1_html=f"results/qc/fastqc/{SAMPLE}_R1_fastqc.html",
        R1_zip=f"results/qc/fastqc/{SAMPLE}_R1_fastqc.zip",
        R2_html=f"results/qc/fastqc/{SAMPLE}_R2_fastqc.html",
        R2_zip=f"results/qc/fastqc/{SAMPLE}_R2_fastqc.zip"
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

