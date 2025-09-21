Snakemake Pipeline: Bacterial WGS QC → Trimming → Assembly → MLST → Variant Calling (optional)

Overview
- Purpose: Analyze a single bacterial isolate from raw paired-end FASTQ through QC, trimming, assembly, MLST typing, and optional alignment + variant calling. Produces concise reports for coverage and variant calls.
- Tools (via conda envs): FastQC, fastp, SPAdes, mlst, BWA/Bowtie2, Samtools, Bcftools. Reference FASTA can be auto-downloaded.
- Reproducibility: Each rule runs in its own pinned conda environment under `envs/` using Snakemake’s `--use-conda`.

 Portfolio
- Clinical-focused project portfolio: see `portfolio_bacterial_wgs_pipeline.md` for a scientific, clinician-oriented summary of goals, methods, results, and interpretation.

Project Layout
- Workflow: `Snakefile` (root), per-tool env YAMLs under `envs/`.
- Inputs: place FASTQs under `data/`; reference fetched to `data/refgenome/` (or set in `config.yaml`).
- Outputs: under `results/` by step: `qc/`, `trimmed/`, `assembly/`, `mlst/`, `aln/`, `variants/`, `reports/`, and `logs/`.

Setup
- Create a minimal driver env: `mamba env create -f envs/snakemake.yml && mamba activate snake`
- Optional quick checks: `bash scripts/tests/test_help.sh`

 Configuration
- Edit `config.yaml` to point to your inputs and options:
  - `sample`: output prefix for this sample.
  - `fastq_r1`, `fastq_r2`: paths to paired FASTQ files under `data/`.
  - `enable_variant_calling`: `true|false` to toggle alignment/VCF steps and reports.
  - `reference_fasta`: reference genome FASTA path; or keep the default and let the pipeline download an E. coli reference.
  - `aligner`: `bwa` or `bowtie2`; `ploidy` (typically 1 for bacteria).
  - `filter_min_qual`, `filter_min_dp`: thresholds for bcftools filtering (defaults 20 and 5).
  - `fastp_extra_args`: optional extra flags appended to fastp (e.g., `--trim_front1 5 --length_required 50`).

Run
- Dry-run (graph/paths): `snakemake -n --use-conda --cores 2`
- Full run (example 4 cores): `snakemake --use-conda --cores 4`

Key Outputs (paths)
- QC: `results/qc/fastqc/<sample>_R1_fastqc.html`, `results/qc/fastqc/<sample>_R2_fastqc.html`
- Trimming: `results/trimmed/<sample>_R1_trimmed.fastq.gz`, `results/trimmed/<sample>_R2_trimmed.fastq.gz`
- QC after trimming: `results/qc/fastqc_trimmed/<sample>_R1_trimmed_fastqc.html`, `results/qc/fastqc_trimmed/<sample>_R2_trimmed_fastqc.html`
- Assembly: `results/assembly/<sample>/contigs.fasta`
- MLST: `results/mlst/<sample>_mlst.tsv`
- Alignment (if enabled): `results/aln/<sample>.sorted.bam` (+ `.bai`)
- Variants (if enabled): `results/variants/<sample>.filtered.vcf.gz` (+ `.tbi`)
- Reports (if enabled):
  - Coverage summary: `results/reports/coverage.txt` (flagstat + coverage)
  - Variant summary (TSV): `results/reports/variants_summary.tsv` (CHROM, POS, REF, ALT, QUAL, DP)
  - VCF stats: `results/reports/vcf_stats.txt` (bcftools stats)
  - MultiQC: `results/reports/multiqc_report.html` (aggregated QC for FastQC/fastp and more)
  - Versions: `results/reports/versions.txt` (Snakemake version and pinned env dependencies)
  - Annotated VCF (if SnpEff enabled): `results/reports/variants_annotated.vcf.gz` (+ `.tbi`)

Optional Features
- Species identification (Mash):
  - Enable via `species_id_enabled: true` in `config.yaml` and provide small reference FASTAs in `species_refs_dir` (default `data/species_refs/`).
  - Outputs: `results/reports/species_id.tsv` (Mash distances, sorted), `results/reports/species_id.txt` (top hit).
- Resistance annotation (BED):
  - Enabled by default. Uses `resistance_genes_bed` from `config.yaml` (default `data/resistance_genes.bed`).
  - Output: `results/reports/resistance_hits.tsv` (VCF-derived variant loci intersected with BED features).
  - Building a BED from public sources (CARD/ResFinder):
    - Download an annotation (GFF3) for your exact reference genome and a list of resistance genes of interest (e.g., from CARD/ResFinder).
    - Convert GFF to BED and optionally filter by gene names using the included helper:
      - Example: `python scripts/gff_to_bed.py data/refgenome/your_ref.gff3 data/resistance_genes.bed --feature gene --name-attr Name --names data/gene_list.txt`
      - Or filter by regex: `--regex "bla|tet|aac\(6'\)-Ib"`
    - Set `resistance_genes_bed` in `config.yaml` to the resulting BED. Optionally set `resistance_genes_bed_sha256` to lock content and enable checksum validation.
- Reference integrity check:
  - Optionally set `reference_sha256` in `config.yaml`; the download step verifies checksum after fetch.
  - Optionally set `resistance_genes_bed_sha256` to verify the BED used for resistance annotation.

Optional Assets Setup
- Scaffold helper files (no data included): `bash scripts/setup_optional_assets.sh`
  - Creates `data/species_refs/` with `.gitkeep` and a README describing expected FASTAs.
  - Creates `data/resistance_genes.bed.example` you can copy/edit to your reference.
  - Creates `data/resistance_genes.bed` (default path used by config) with the example contents.

Data Notes
- Do not commit large data. `.gitignore` already excludes common bioinformatics files and `data/` contents by default.
- Provide instructions or small fixtures instead of committing raw datasets. For public reads, reference NCBI SRA accessions.

Smoke Tests
- Help and rule listing: `bash scripts/tests/test_help.sh`
- Minimal dry-run targets: `bash scripts/tests/test_minirun.sh`
 - Tiny end-to-end test: `bash scripts/tests/test_e2e_tiny.sh` (creates small FASTQs and runs a real mini pipeline against the toy reference)

Reproducibility
- Each rule specifies `conda:` to pin tool versions from `envs/*.yaml`.
- Run with `--use-conda` so Snakemake creates and uses per-rule environments.
 - Optional checksums: set `reference_sha256` to verify downloaded reference integrity.

CI
- A minimal GitHub Actions workflow `.github/workflows/ci.yml` runs shellcheck, Snakemake lint, and a dry-run on pushes/PRs.
 - Optional manual E2E test via Actions → CI → Run workflow (set `run_e2e=true`).

Notes on Data
- Large raw data should not be committed. Use instructions or scripts to fetch small fixtures locally. The repository’s `.gitignore` excludes typical large bioinformatics files.
- If running air-gapped, set `reference_fasta` in `config.yaml` to a locally available FASTA and optionally set `reference_sha256` for integrity verification.

Performance & Reruns
- Threads are configurable per-rule under `threads:` in `config.yaml`. Defaults are tuned for a laptop.
- For faster code-only iterations, prefer: `snakemake --rerun-triggers mtime`.

Profiles
- Example local profile: `profiles/local/config.yaml`
  - Use via: `snakemake --profile profiles/local`
- Add an HPC profile (e.g., `profiles/slurm/`) tailored to your environment.

Troubleshooting
- Conda solve issues: try `mamba` and clear caches; ensure bioconda + conda-forge channels are configured.
- Permissions/Exec: ensure scripts are executable (`chmod +x scripts/**/*.sh`); on macOS, allow terminal full disk access if needed.
- Java required for SnpEff (provided by env). If memory errors occur, set `_JAVA_OPTIONS` (e.g., `-Xmx2g`).

Notes
- Variant calling and reports are optional; toggle with `enable_variant_calling` in `config.yaml`.
- Reference download rule (`download_reference`) retrieves an E. coli reference if you keep the default.
 - Variant effect annotation (SnpEff):
  - Enable via `snpeff_enabled: true` in `config.yaml`.
  - Choose `snpeff_genome` to match your reference. For a custom reference, set `snpeff_genome: custom` and provide `snpeff_gff` (GFF3 matching your reference). The pipeline builds a local DB under `snpeff_data_dir` (default `results/snpeff`).
  - Outputs: `results/variants/<sample>.annotated.vcf.gz` (+ `.tbi`), a copy in `results/reports/variants_annotated.vcf.gz` (+ `.tbi`), and `results/reports/variants_annotated.tsv` (effect, impact, gene, HGVS.p).
