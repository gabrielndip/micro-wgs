import os


configfile: "config.yaml"

shell.executable("/bin/bash")

SAMPLE   = config["sample"]
FASTQ_R1 = config["fastq_r1"]
FASTQ_R2 = config["fastq_r2"]

# Variant-calling related config (with sensible defaults)
ENABLE_VARCALL   = bool(config.get("enable_variant_calling", False))
REFERENCE_FASTA  = config.get("reference_fasta", "data/refgenome/ecoli_INF32-16-A_ref.fasta")
ALIGNER          = config.get("aligner", "bwa")  # "bwa" or "bowtie2"
MARK_DUPLICATES  = bool(config.get("mark_duplicates", False))
CALLER           = config.get("caller", "bcftools")
PLOIDY           = int(config.get("ploidy", 1))

# Get clean base names (e.g., SRR15334628_1)
R1_BASENAME = os.path.basename(FASTQ_R1).replace(".fastq.gz", "").replace(".fastq", "")
R2_BASENAME = os.path.basename(FASTQ_R2).replace(".fastq.gz", "").replace(".fastq", "")

# Build default targets, optionally extend with variant-calling outputs
ALL_TARGETS = [
    f"results/qc/fastqc/{R1_BASENAME}_fastqc.html",
    f"results/qc/fastqc/{R2_BASENAME}_fastqc.html",
    f"results/qc/fastp/{SAMPLE}_fastp.html",
    f"results/qc/fastp/{SAMPLE}_fastp.json",
    f"results/assembly/{SAMPLE}/contigs.fasta",
    f"results/mlst/{SAMPLE}_mlst.tsv",
    REFERENCE_FASTA,
]

if ENABLE_VARCALL:
    ALL_TARGETS.extend([
        f"results/variants/{SAMPLE}.filtered.vcf.gz",
        f"results/variants/{SAMPLE}.filtered.vcf.gz.tbi",
    ])

rule all:
    input: ALL_TARGETS

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
    """Download and unpack the reference genome for E. coli INF32/16/A (GCA_018310205.1)"""
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
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/310/205/GCA_018310205.1_ASM1831020v1/GCA_018310205.1_ASM1831020v1_genomic.fna.gz \
            >> {log} 2>&1
        gunzip -c data/refgenome/tmp_ref.fna.gz > {output} 2>> {log}
        rm data/refgenome/tmp_ref.fna.gz
        """


# -----------------------
# Alignment & Variant Calling
# -----------------------

# Copy reference into results for indexing and reproducibility
rule copy_reference:
    input:
        REFERENCE_FASTA
    output:
        "results/ref/ref.fa"
    log:
        "results/logs/reference_copy.log"
    shell:
        (
            "mkdir -p results/ref results/logs && "
            "cp -f {input} {output} > {log} 2>&1"
        )

rule faidx_reference:
    input:
        fasta="results/ref/ref.fa"
    output:
        "results/ref/ref.fa.fai"
    log:
        "results/logs/faidx_reference.log"
    conda:
        "envs/samtools_env.yaml"
    shell:
        (
            "samtools faidx {input.fasta} > /dev/null 2>> {log}"
        )

rule bwa_index:
    input:
        fasta="results/ref/ref.fa"
    output:
        expand("results/ref/ref.fa.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"])
    log:
        "results/logs/bwa_index.log"
    conda:
        "envs/bwa_env.yaml"
    shell:
        (
            "bwa index {input.fasta} > {log} 2>&1"
        )

rule bowtie2_index:
    input:
        fasta="results/ref/ref.fa"
    output:
        expand("results/ref/bowtie2/ref.{i}.bt2", i=[1,2,3,4,5,6])
    log:
        "results/logs/bowtie2_index.log"
    conda:
        "envs/bowtie2_env.yaml"
    shell:
        (
            "mkdir -p results/ref/bowtie2 && "
            "bowtie2-build {input.fasta} results/ref/bowtie2/ref > {log} 2>&1"
        )

rule align_bwa:
    input:
        ref="results/ref/ref.fa",
        idx=expand("results/ref/ref.fa.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"]),
        R1=f"results/trimmed/{{sample}}_R1_trimmed.fastq.gz",
        R2=f"results/trimmed/{{sample}}_R2_trimmed.fastq.gz"
    output:
        "results/aln/{sample}.bwa.sam"
    params:
        sample=SAMPLE
    threads: 4
    log:
        "results/logs/bwa_align_{sample}.log"
    conda:
        "envs/bwa_env.yaml"
    shell:
        (
            "mkdir -p results/aln results/logs && "
            "bwa mem -t {threads} {input.ref} {input.R1} {input.R2} > {output} 2> {log}"
        )

rule align_bowtie2:
    input:
        idx=expand("results/ref/bowtie2/ref.{i}.bt2", i=[1,2,3,4,5,6]),
        R1=f"results/trimmed/{{sample}}_R1_trimmed.fastq.gz",
        R2=f"results/trimmed/{{sample}}_R2_trimmed.fastq.gz"
    output:
        "results/aln/{sample}.bowtie2.sam"
    params:
        index_prefix="results/ref/bowtie2/ref"
    threads: 4
    log:
        "results/logs/bowtie2_align_{sample}.log"
    conda:
        "envs/bowtie2_env.yaml"
    shell:
        (
            "mkdir -p results/aln results/logs && "
            "bowtie2 -x {params.index_prefix} -1 {input.R1} -2 {input.R2} -p {threads} -S {output} 2> {log}"
        )

rule select_alignment:
    input:
        lambda wildcards: f"results/aln/{SAMPLE}.{'bwa' if ALIGNER == 'bwa' else 'bowtie2'}.sam"
    output:
        f"results/aln/{SAMPLE}.sam"
    log:
        "results/logs/select_alignment.log"
    shell:
        (
            "ln -sf {input} {output} > {log} 2>&1 || cp -f {input} {output} >> {log} 2>&1"
        )

rule sort_bam:
    input:
        sam=f"results/aln/{SAMPLE}.sam"
    output:
        bam=f"results/aln/{SAMPLE}.sorted.bam"
    threads: 4
    log:
        f"results/logs/samtools_sort__{SAMPLE}.log"
    conda:
        "envs/samtools_env.yaml"
    shell:
        (
            "samtools view -@ {threads} -bS {input.sam} | samtools sort -@ {threads} -o {output.bam} - 2> {log}"
        )

rule index_bam:
    input:
        bam=f"results/aln/{SAMPLE}.sorted.bam"
    output:
        bai=f"results/aln/{SAMPLE}.sorted.bam.bai"
    log:
        f"results/logs/samtools_index__{SAMPLE}.log"
    conda:
        "envs/samtools_env.yaml"
    shell:
        (
            "samtools index {input.bam} > {log} 2>&1"
        )

rule mark_duplicates:
    input:
        bam=f"results/aln/{SAMPLE}.sorted.bam",
        bai=f"results/aln/{SAMPLE}.sorted.bam.bai"
    output:
        bam=f"results/aln/{SAMPLE}.mkdup.bam",
        bai=f"results/aln/{SAMPLE}.mkdup.bam.bai"
    log:
        f"results/logs/samtools_markdup__{SAMPLE}.log"
    conda:
        "envs/samtools_env.yaml"
    shell:
        (
            "samtools markdup -r {input.bam} {output.bam} > {log} 2>&1 && "
            "samtools index {output.bam} >> {log} 2>&1"
        )

def get_bam_for_calling():
    if MARK_DUPLICATES:
        return f"results/aln/{SAMPLE}.mkdup.bam"
    return f"results/aln/{SAMPLE}.sorted.bam"

rule bcftools_call_raw:
    input:
        bam=lambda wildcards: get_bam_for_calling(),
        bai=lambda wildcards: get_bam_for_calling() + ".bai",
        ref="results/ref/ref.fa",
        fai="results/ref/ref.fa.fai"
    output:
        vcf="results/variants/{sample}.raw.vcf.gz"
    params:
        ploidy=PLOIDY
    threads: 2
    log:
        "results/logs/bcftools_call_raw__{sample}.log"
    conda:
        "envs/bcftools_env.yaml"
    shell:
        (
            "mkdir -p results/variants results/logs && "
            "bcftools mpileup -Ou -f {input.ref} {input.bam} 2> {log} | "
            "bcftools call -mv --ploidy {params.ploidy} -Oz -o {output.vcf} >> {log} 2>&1"
        )

rule index_raw_vcf:
    input:
        vcf="results/variants/{sample}.raw.vcf.gz"
    output:
        tbi="results/variants/{sample}.raw.vcf.gz.tbi"
    log:
        "results/logs/index_raw_vcf__{sample}.log"
    conda:
        "envs/bcftools_env.yaml"
    shell:
        "bcftools index -t {input.vcf} > {log} 2>&1"

rule filter_vcf:
    input:
        vcf="results/variants/{sample}.raw.vcf.gz",
        tbi="results/variants/{sample}.raw.vcf.gz.tbi"
    output:
        vcf="results/variants/{sample}.filtered.vcf.gz",
        tbi="results/variants/{sample}.filtered.vcf.gz.tbi"
    log:
        "results/logs/filter_vcf__{sample}.log"
    conda:
        "envs/bcftools_env.yaml"
    shell:
        (
            "bcftools filter -e 'QUAL<20 || DP<5' -Oz -o {output.vcf} {input.vcf} > {log} 2>&1 && "
            "bcftools index -t {output.vcf} >> {log} 2>&1 && "
            "test -f {output.tbi} || ln -sf {output.vcf}.tbi {output.tbi} >> {log} 2>&1 || true"
        )
