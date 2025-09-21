## Phase 1 — Workaround: Run Apptainer in a Linux VM on macOS (ARM64)

Why: Apptainer does not run natively on macOS. I run a lightweight ARM64 Ubuntu VM and do all Apptainer builds/runs inside it. This avoids cross‑arch headaches and keeps performance good on an Apple Silicon Mac.

- Option A (what I use): Lima (lightweight, simple host‑file mounts)
- Option B: Multipass (Canonical’s Ubuntu VM manager)

I show Lima below; I include Multipass commands afterward.

1. Install Lima on macOS
- `brew install lima qemu`
- `limactl start --name=apptainer-ubuntu template://ubuntu-lts`

2. Enter the VM shell
- `limactl shell apptainer-ubuntu`

3. Install Apptainer and build deps (inside the VM)
- `sudo apt-get update`
- `sudo apt-get install -y apptainer squashfs-tools git build-essential`

Note: On Ubuntu 22.04+ the package is named `apptainer`. If it’s not in the mirror, I install from Sylabs packages or build from source. Using `apt` on Ubuntu LTS is the simplest.

4. Make my project directory visible inside the VM
- Lima auto‑mounts `$HOME`. In the VM, my macOS `$HOME` appears under `/Users/<me>`. I navigate to my repo:
  - `cd /Users/<me>/Projects/bioinfo-pipelines/micro-wgs`
- For faster IO, I often copy the repo into a VM‑local path:
  - `rsync -av --delete /Users/<me>/Projects/bioinfo-pipelines/micro-wgs/ ~/micro-wgs/`
  - `cd ~/micro-wgs`

Multipass alternative
- `brew install --cask multipass`
- `multipass launch --name apptainer --mem 6G --disk 30G 22.04`
- `multipass shell apptainer`
- Inside VM: `sudo apt-get update && sudo apt-get install -y apptainer squashfs-tools git`
- Mount my project: `multipass transfer /path/to/micro-wgs apptainer:/home/ubuntu/micro-wgs`


## Phase 2 — Containerization Strategy (ARM64‑conscious)

I have a working Snakemake pipeline and an `environment.yml`. The easiest path for me is a single “runner” container that:
- Includes Snakemake and my tools resolved from `environment.yml` (using micromamba)
- Executes Snakemake inside the container
- Binds my project’s working directory for I/O

This avoids per‑rule container declarations and complex refactors. I can add per‑rule containers later if needed.

Caveats:
- On ARM64, some bio‑tools may not have linux‑aarch64 conda packages. If any tool fails to resolve for ARM, I build and run on x86_64 via CI for the HPC target (see Phase 6).


## Phase 3 — Create an Apptainer Definition (from environment.yml)

I create `containers/pipeline.def` at the repo root with:

- Bootstrap from a multi-arch base (micromamba is multi-arch)
- Install Snakemake + pipeline deps from `environment.yml` at build time
- Set PATH and locales
- Create a non-root workdir

Example definition (assumes `environment.yml` at the project root):

```
Bootstrap: docker
From: mambaorg/micromamba:1.5.8

%labels
    Maintainer my.name@example.com
    Description Snakemake pipeline (ARM64-ready) built via micromamba

%environment
    export MAMBA_ROOT_PREFIX=/opt/micromamba
    export PATH=/opt/micromamba/bin:$PATH
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%post
    set -euo pipefail
    # Create env from environment.yml (name is taken from the file or set here)
    micromamba install -y -n pipeline -f /work/environment.yml -c conda-forge -c bioconda
    micromamba clean -a -y
    # Snakemake (pin a version compatible with your workflow)
    micromamba install -y -n pipeline snakemake=7.32.4
    micromamba clean -a -y

    # Convenience symlink to activate env by default
    ln -s /opt/micromamba/envs/pipeline /opt/env

    # Runtime user/work setup (owner remains the caller in Apptainer)
    mkdir -p /opt/work
    chmod 1777 /opt/work

%runscript
    # Default entrypoint: run snakemake in the bound workdir
    exec /opt/env/bin/snakemake "$@"

%apprun snakemake
    exec /opt/env/bin/snakemake "$@"

%apprun bash
    exec /bin/bash "$@"

%files
    # The project directory will be bound at runtime; we only need environment.yml at build time if desired.
```

Notes:
- Using `/work` as the working directory (see run commands) keeps paths simple.
- If my pipeline uses many envs (per‑rule `.yaml` files), I first try merging the packages into one `environment.yml`. I keep it simple, and only refactor to per‑rule images if needed.


## Phase 4 — Build and Smoke‑Test the ARM64 SIF (inside the VM)

1. Build SIF
- `cd ~/micro-wgs`
- `mkdir -p containers`
- Move the definition above into `containers/pipeline.def`
- `apptainer build containers/pipeline_arm64.sif containers/pipeline.def`

2. Smoke test (Snakemake version)
- `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --version`

3. Dry‑run my workflow inside the container
- `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake -n --cores 4`

4. Full run (adjust cores and any required CLI flags)
- `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --cores 4`

Tips:
- I bind more paths as needed (e.g., scratch, reference volumes):
  - `apptainer exec -B /scratch:/scratch -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --cores 4`


## Phase 5 — (Optional) Use Snakemake’s Apptainer Integration Later

If/when I migrate to per‑rule containers:
- I add `container: "docker://..."` or `"oras://..."` lines per rule
- I run with Snakemake’s Apptainer integration:
  - `snakemake --use-apptainer --cores 4`
- I pass Apptainer args via:
  - `snakemake --use-apptainer --apptainer-args "-B $PWD:/work -W /work --containall"`

For now, the simplest approach I use is “run Snakemake inside the monolithic container,” as shown above.


## Phase 6 — Produce an x86_64 SIF for HPC via CI (if my cluster is x86)

If my deployment target is an x86_64 cluster, I produce the x86_64 SIF in CI on Linux/AMD64 (fast and native), even if I develop on ARM64 locally.

1. Add a GitHub Actions workflow (e.g., `.github/workflows/sif.yml`):

```
name: Build SIF

on:
  push:
    branches: [ main, chore/**, feat/**, fix/** ]
  workflow_dispatch: {}

jobs:
  build-sif-amd64:
    runs-on: ubuntu-22.04
    steps:
      - uses: actions/checkout@v4
      - name: Install Apptainer
        run: |
          sudo apt-get update
          sudo apt-get install -y apptainer squashfs-tools
      - name: Build SIF
        run: |
          mkdir -p containers
          apptainer build containers/pipeline_amd64.sif containers/pipeline.def
      - name: Upload artifact
        uses: actions/upload-artifact@v4
        with:
          name: pipeline_amd64.sif
          path: containers/pipeline_amd64.sif
```

2. I download the artifact and copy it to my cluster.

3. Verify and run on the cluster:
- `apptainer exec pipeline_amd64.sif snakemake --version`
- `apptainer exec -B "$PWD":/work -W /work pipeline_amd64.sif snakemake --cores 8`

Signing (recommended for provenance)
- `apptainer key newpair`
- `apptainer sign containers/pipeline_amd64.sif`
- `apptainer verify containers/pipeline_amd64.sif`


## Phase 7 — Data Binding, Permissions, and Reproducibility

- I bind working directories explicitly:
  - `-B "$PWD":/work -W /work`
- I keep host data read‑only where appropriate:
  - `-B "$PWD":/work:ro` (if I only need reads)
- I use containment for cleanliness:
  - `--containall --no-home`
- I record exact Snakemake and environment versions (my pipeline already writes `results/reports/versions.txt`).
- I avoid writing to `/` in the container; I use `-W /work` and write under the repo’s `results/` folder.


## Phase 8 — ARM64 Caveats and Fallbacks

- Package availability: Some Bioconda packages lack linux‑aarch64 builds. If micromamba fails to solve for ARM64, I either:
  1) Swap/pin packages to ARM‑friendly alternatives (e.g., use bwa‑mem2 vs bwa if supported), or
  2) Build the SIF on x86_64 in CI (Phase 6) and run on x86_64 targets (HPC). For local ARM testing, I still use the monolithic container but limit tool usage to what resolves.
- Performance: A native ARM64 VM (Lima) is efficient. I avoid cross‑arch emulation locally for routine dev and use CI for x86 builds.
- File paths: On Lima, my macOS home is mounted under `/Users/<me>`. I prefer syncing the repo into `~/micro-wgs` in the VM for consistent Linux paths.


## Phase 9 — Security & Best Practices

- I keep `environment.yml` pinned (channels: conda‑forge, bioconda).
- I use Apptainer in unprivileged mode (default) for safer execution.
- I sign SIF images and verify on target systems.
- I avoid embedding secrets in the image; I mount credentials at runtime if needed (e.g., read‑only tokens).
- I use read‑only binds for inputs and write outputs to a dedicated `results/` directory.


## Phase 10 — Quickstart Summary

- Start VM and enter it:
  - `limactl start --name=apptainer-ubuntu template://ubuntu-lts`
  - `limactl shell apptainer-ubuntu`
- Install Apptainer:
  - `sudo apt-get update && sudo apt-get install -y apptainer squashfs-tools`
- Build SIF (in `~/micro-wgs`):
  - `apptainer build containers/pipeline_arm64.sif containers/pipeline.def`
- Dry‑run pipeline:
  - `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake -n --cores 4`
- Run pipeline:
  - `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --cores 4`
- For x86_64 HPC: I use GitHub Actions to build `pipeline_amd64.sif` and run it on the cluster.

This plan keeps things practical for me: I develop and test locally in an ARM64 VM, package everything in a single Apptainer image from my `environment.yml`, and, if needed, produce an x86_64 SIF in CI for deployment to typical HPC environments.
