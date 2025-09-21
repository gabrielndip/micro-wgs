## 1. Project Title
Reproducible Bacterial WGS Workflow for QC, Trimming, Assembly, MLST, and Optional Variant Calling

## 2. Abstract
Whole-genome sequencing (WGS) of bacterial isolates is increasingly used in clinical microbiology for rapid species/sequence-type confirmation, outbreak investigation, and detection of variants with potential impact on antimicrobial therapy or infection control. This project delivers a Snakemake-based, per-rule conda–pinned workflow that takes a single isolate from raw paired-end FASTQ through quality control (QC), adapter/quality trimming, de novo assembly, MLST typing, and optional read alignment with variant calling and clinical-style reporting. Integrated tools include FastQC/fastp (QC/trim), SPAdes (assembly), mlst (PubMLST typing), BWA/Bowtie2+Samtools (alignment/BAM), Bcftools (calling/filtering/statistics), and MultiQC (aggregated QC). Optional modules support Mash-based species identification, resistance-locus intersection (bedtools) against a user-curated panel, and SnpEff effect annotation when appropriate annotations are available.

Applied to the included sample (“sample”), alignment to an E. coli reference achieved high mapping rates (~98.5% mapped; ~74.7% properly paired). The workflow produced a robust assembly, an Achtman E. coli MLST assignment, and a filtered variant set summarized by site-level metrics. Aggregate statistics show 283 total variants (90 SNPs, 193 indels; Ts/Tv ~2.6). Outputs include coverage summaries, variant tables, QC dashboards, and tool/environment provenance suitable for technical review and audit.

Clinically relevant considerations are highlighted: explicit QC artifacts to support acceptance criteria (e.g., mapping rate, properly paired percentage, depth), optional resistance-locus reporting against a curated BED, and fully pinned software environments for traceability. This work demonstrates a concise, reproducible WGS pipeline aligned with laboratory medicine needs; it is intended for research and method-development use and would require formal validation, verification, and governance (e.g., CLIA/CAP) before diagnostic deployment.

## 3. Introduction
Bacterial WGS supports patient care and public health through rapid species confirmation, determination of sequence type for surveillance, and identification of variants relevant to antimicrobial resistance and transmission. In laboratory medicine, methods must be reproducible, auditable, and supported by transparent QC and provenance.

Research aims:
- Implement a transparent, end-to-end WGS pipeline for a single bacterial isolate, from raw reads to assembly, MLST, and optional variant calling.
- Ensure reproducibility via per-step, pinned conda environments and deterministic outputs.
- Provide configurable options for species identification, resistance-locus intersection, and variant effect annotation suitable for clinical method development.

Significance:
- Aligns with clinical needs for QCed, traceable, and reviewable WGS analyses.
- Produces outputs that map to common clinical questions (species/ST confirmation, coverage adequacy, curated resistance-locus checks).
- Serves as a foundation for validation and scale-up to multi-sample or HPC deployments.

## 4. Methods
Assumption: Short-read Illumina paired-end WGS for a single bacterial isolate with ploidy 1.

- Data Source and Pre-processing
  - Inputs: Paired FASTQs under `data/` (configured via `config.yaml`: `data/SRR15334628_1.fastq`, `data/SRR15334628_2.fastq`).
  - Initial QC: FastQC on raw reads to assess base quality, adapter content, and per-sequence metrics (rule `fastqc_raw_reads`).
  - Trimming: fastp for adapter and quality trimming; outputs gzipped paired reads plus HTML/JSON QC (rule `trim_reads_fastp`). A second FastQC run evaluates trimmed reads (rule `fastqc_trimmed_reads`).
  - Clinical acceptance checks (example guidance; adjust per SOP): mapping rate target >95%, properly paired >70–80%, and median depth sufficient for intended reporting (e.g., ≥30× for routine SNP detection). MultiQC aggregates these indicators.

- Core Analysis Pipeline
  - Assembly: SPAdes assembles trimmed reads to contigs (`results/assembly/<sample>/contigs.fasta`) (rule `assemble_genome_spades`).
  - MLST Typing: `mlst` scans contigs against PubMLST schemes; outputs TSV with scheme and allele calls (rule `mlst_typing`).
  - Reference Acquisition: If not provided, an E. coli reference FASTA is auto-downloaded (rule `download_reference`) into `data/refgenome/`; integrity verification via optional SHA256.
  - Alignment (optional, enabled here): BWA or Bowtie2 index and alignment, selected via `aligner` in `config.yaml` (rules `bwa_index`/`bowtie2_index` and `align_bwa`/`align_bowtie2`). Samtools converts/sorts/indexes BAM; optional duplicate marking (rules `sort_bam`, `index_bam`, `mark_duplicates`).
  - Variant Calling (optional, enabled here): Bcftools mpileup/call produces raw VCF; filtering excludes low-quality/low-depth sites (rules `bcftools_call_raw`, `index_raw_vcf`, `filter_vcf`). Default thresholds: `QUAL >= 20`, `DP >= 5` (configurable).
  - Reporting:
    - Coverage: `samtools flagstat` and `samtools coverage` summarized to `results/reports/coverage.txt` (rule `coverage_summary`).
    - Variant Summary: Bcftools query exports CHROM, POS, REF, ALT, QUAL, DP to `results/reports/variants_summary.tsv` (rule `variants_summary`).
    - VCF Stats: `bcftools stats` to `results/reports/vcf_stats.txt` (rule `vcf_stats`).
    - MultiQC aggregates QC artifacts into `results/reports/multiqc_report.html` (rule `multiqc`).
    - Versions: Snakemake version and pinned environment snippets to `results/reports/versions.txt` (rule `pipeline_versions`).

- Optional Analyses
  - Species Identification: Mash sketch of reference set in `data/species_refs/` and of assembled contigs; distances reported in `results/reports/species_id.tsv` and top hit in `results/reports/species_id.txt` (rules `mash_sketch_refs`, `mash_sketch_sample`, `mash_species_id`). Note: accuracy depends on the provided reference set.
  - Resistance Annotation: Convert filtered VCF to BED and intersect with a curated resistance BED (`data/resistance_genes.bed`), outputting putative hits in `results/reports/resistance_hits.tsv` (rules `variants_to_bed`, `validate_resistance_bed`, `resistance_intersect`). For clinical use, populate this BED from trusted sources (e.g., CARD/ResFinder/AMRFinderPlus-aligned panels) matched to the exact reference; optional SHA256 supports change control.
  - Variant Effect Annotation: If enabled, builds a local SnpEff database for custom genomes or uses a named genome; annotated VCF and TSV are written to `results/variants/` and `results/reports/` (rules `snpeff_build_db`, `snpeff_annotate`, `snpeff_summary`, `annotated_vcf_to_reports`). Disabled by default in `config.yaml`.

- Statistical Analysis & Key Parameters
  - Variant calling: `bcftools mpileup | bcftools call -mv --ploidy 1`.
  - Filtering: Exclude records failing `QUAL < 20 || DP < 5` (configurable: `filter_min_qual`, `filter_min_dp`).
  - Ploidy: Set to 1 (bacterial).
  - Alignment: BWA-MEM or Bowtie2 with multi-threading; duplicate marking optional.
  - QC Aggregation: MultiQC summarizes FastQC/fastp and other metrics.
  - Clinical validation notes (method-development guidance): define acceptance criteria prospectively (e.g., mapping rate, min depth, callable-genome fraction), establish LoD and reproducibility across runs, and include contamination controls and negative/positive controls as per SOPs.

- Code and Environment
  - Orchestration: Snakemake workflow in `Snakefile` with `--use-conda`.
  - Languages: Bash (rules), Python (utility script `scripts/gff_to_bed.py`).
  - Environments: Pinned conda YAMLs in `envs/` (e.g., `fastqc_env.yaml`, `fastp_env.yaml`, `spades_env.yaml`, `mlst_env.yaml`, `bwa_env.yaml`, `bowtie2_env.yaml`, `samtools_env.yaml`, `bcftools_env.yaml`, `multiqc_env.yaml`, `mash_env.yaml`, `snpeff_env.yaml`, `wget_env.yaml`). A minimal driver env is defined in `envs/snakemake.yml`.
  - Governance: per-rule environments and `results/reports/versions.txt` enable auditability and traceability; configure data retention and PHI/PII safeguards per institutional policy.

## 5. Results
- QC and Trimming
  - Raw and trimmed FastQC reports generated under `results/qc/fastqc/` and `results/qc/fastqc_trimmed/`.
  - fastp reports available at `results/qc/fastp/sample_fastp.html` and `.json`.
  - MultiQC summary: `results/reports/multiqc_report.html`.

- Assembly and Typing
  - Assembled contigs written to `results/assembly/sample/contigs.fasta`.
  - MLST result indicates an E. coli Achtman scheme assignment with a specific sequence type, reported in `results/mlst/sample_mlst.tsv`.

- Alignment and Coverage (enabled)
  - Sorted and indexed BAM: `results/aln/sample.sorted.bam` (+ `.bai`).
  - Coverage summary in `results/reports/coverage.txt` indicates high mapping (~98.5%) and a substantial fraction of properly paired reads (~74.7%). Depth across major contigs is high (dozens to >100×), consistent with robust coverage.

- Variant Calling and Summary (enabled)
  - Filtered VCF and index: `results/variants/sample.filtered.vcf.gz` (+ `.tbi`).
  - Variant counts from `results/reports/vcf_stats.txt`: 283 total variants (90 SNPs, 193 indels); transition/transversion ratio ~2.6.
  - Tabulated variant summary (position, alleles, QUAL, DP): `results/reports/variants_summary.tsv`.

- Optional Modules
  - Species ID: Outputs at `results/reports/species_id.tsv` and `results/reports/species_id.txt`. Note that, with the toy reference set provided in `data/species_refs/`, the top hit reflects the toy entry and is for demonstration only.
  - Resistance Intersect: `results/reports/resistance_hits.tsv` contained no hits with the example BED (`data/resistance_genes.bed`), which uses placeholder coordinates for demonstration.

- Provenance
  - Reference FASTA downloaded and tracked in `data/refgenome/ecoli_INF32-16-A_ref.fasta`.
  - Tool and environment versions recorded in `results/reports/versions.txt`.
  - Per-rule logs available in `results/logs/`.

## 6. Discussion and Conclusion
This workflow implements an end-to-end bacterial WGS analysis with modular, reproducible steps: rigorous QC, trimming, de novo assembly, MLST typing, optional alignment and variant calling, and comprehensive reporting. For the provided sample, high mapping rates and deep coverage enabled reliable variant detection, yielding 283 high-confidence variants after filtering. The MLST result identifies the isolate within a standard E. coli typing scheme, facilitating downstream epidemiological or comparative analyses. While resistance-locus intersection produced no hits using the example BED, this result reflects the demonstrative nature of the provided annotation; use of a curated, reference-matched resistance BED is recommended for biological interpretation. Similarly, Mash species identification depends on the supplied reference set; the toy references in `data/species_refs/` are intended only to exercise the method.

Clinical interpretation and practice: define laboratory acceptance criteria (e.g., mapping rate, depth, coverage uniformity), confirm species/ST using validated databases, and interpret variants in the context of curated knowledge bases and phenotypic AST when available. Reporting should clearly state methods, software versions, and limitations, and distinguish research-use-only outputs from clinically validated results.

Limitations include single-sample scope, dependence on reference choice (affecting callability/interpretation), and default filters that may require retuning for specific organisms or specimen types. Variant effect annotation requires appropriate genome annotations and curated annotations for resistance/virulence.

Future directions involve method validation (accuracy, precision, LoD), external quality assessment, scaling to multi-sample batches with controls, integrating assembly QC metrics, consensus genome generation, curated AMR panels, and deployment via institutional HPC/cloud with audited profiles.

Conclusion: The project delivers a clear, reproducible WGS pipeline that demonstrates strong bioinformatics engineering practices and produces interpretable bacterial genomics outputs (QC, assembly, typing, and variants) with robust reporting and provenance.

## 7. Bioinformatics Competencies Showcase
- Programming & Scripting
  - Bash (Snakemake shell directives, test scripts)
  - Python (argument-parsed utility: `scripts/gff_to_bed.py`)
- Bioinformatics Tools
  - QC/Trimming: FastQC, fastp, MultiQC
  - Assembly/Typing: SPAdes, mlst (PubMLST)
  - Alignment/BAM: BWA, Bowtie2, Samtools
  - Variant Calling/Processing: Bcftools (mpileup/call, filter, stats), Tabix
  - Optional: Mash (species ID), bedtools (intersect), SnpEff/SnpSift (annotation)
- Core Libraries
  - Python standard library (argparse, file I/O)
  - Conda-managed toolchain via Bioconda/conda-forge
- Key Concepts & Skills
  - NGS read QC, adapter/quality trimming
  - De novo assembly and sequence typing (MLST)
  - Reference indexing, short-read alignment, BAM processing
  - Variant calling, filtering strategies, and summary reporting
  - Workflow management (Snakemake), per-rule conda environments
  - Reproducibility, provenance tracking, and environment pinning
  - Input validation, integrity checks (optional SHA256)
  - Modular pipeline design with optional analyses and clear file organization

Paths referenced include: `Snakefile`, `config.yaml`, `envs/*.yaml`, `scripts/*.sh|*.py`, and outputs under `results/`. All counts and summaries are derived from generated reports such as `results/reports/coverage.txt`, `results/reports/vcf_stats.txt`, and `results/mlst/sample_mlst.tsv`. No raw sequences or full tables are displayed.
