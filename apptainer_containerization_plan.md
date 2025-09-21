## Phase 1 — Workaround: Run Apptainer in a Linux VM on macOS (ARM64)

Why: Apptainer does not run natively on macOS. We’ll run a lightweight ARM64 Ubuntu VM and do all Apptainer builds/runs inside it. This avoids cross-arch headaches and keeps performance good on an Apple Silicon Mac.

- Option A (recommended): Lima (lightweight, simple host-file mounts)
- Option B: Multipass (Canonical’s Ubuntu VM manager)

Pick one. Below shows Lima; Multipass commands are included after.

1. Install Lima on macOS
- `brew install lima qemu`
- `limactl start --name=apptainer-ubuntu template://ubuntu-lts`

2. Enter the VM shell
- `limactl shell apptainer-ubuntu`

3. Install Apptainer and build deps (inside the VM)
- `sudo apt-get update`
- `sudo apt-get install -y apptainer squashfs-tools git build-essential`

Note: On Ubuntu 22.04+ the package is named `apptainer`. If unavailable in your mirror, install from Sylabs packages or build from source. Keep it simple: `apt` is fine on current Ubuntu LTS.

4. Ensure your project directory is visible inside the VM
- Lima auto-mounts `$HOME`. In the VM, your macOS `$HOME` appears under `/Users/<you>`. Navigate to your repo:
  - `cd /Users/<you>/Projects/bioinfo-pipelines/micro-wgs`
- If needed, copy the repo into a VM-local path for faster IO:
  - `rsync -av --delete /Users/<you>/Projects/bioinfo-pipelines/micro-wgs/ ~/micro-wgs/`
  - `cd ~/micro-wgs`

Multipass alternative
- `brew install --cask multipass`
- `multipass launch --name apptainer --mem 6G --disk 30G 22.04`
- `multipass shell apptainer`
- Inside VM: `sudo apt-get update && sudo apt-get install -y apptainer squashfs-tools git`
- Mount your project: `multipass transfer /path/to/micro-wgs apptainer:/home/ubuntu/micro-wgs`


## Phase 2 — Containerization Strategy (ARM64-conscious)

You have a working Snakemake pipeline and an `environment.yml`. The easiest path for a container novice is a single “runner” container that:
- Includes Snakemake and your tools resolved from `environment.yml` (using micromamba)
- Executes Snakemake inside the container
- Binds your project’s working directory for I/O

This avoids per-rule container declarations and complex refactors. You can add per-rule containers later.

Caveats:
- On ARM64, some bio-tools may not have linux-aarch64 conda packages. If any tool fails to resolve for ARM, build and run on x86_64 via CI for the HPC target (see Phase 6).


## Phase 3 — Create an Apptainer Definition (from environment.yml)

Create `containers/pipeline.def` at repo root with:

- Bootstrap from a multi-arch base (micromamba is multi-arch)
- Install Snakemake + pipeline deps from `environment.yml` at build time
- Set PATH and locales
- Create a non-root workdir

Example definition (assumes `environment.yml` at project root):

```
Bootstrap: docker
From: mambaorg/micromamba:1.5.8

%labels
    Maintainer your.name@example.com
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
- If your pipeline uses many envs (per-rule `.yaml` files), consider merging the packages into one `environment.yml` first. Start simple; refactor per-rule images later if needed.


## Phase 4 — Build and Smoke-Test the ARM64 SIF (inside the VM)

1. Build SIF
- `cd ~/micro-wgs`
- `mkdir -p containers`
- Move the definition above into `containers/pipeline.def`
- `apptainer build containers/pipeline_arm64.sif containers/pipeline.def`

2. Smoke test (Snakemake version)
- `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --version`

3. Dry-run your workflow inside the container
- `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake -n --cores 4`

4. Full run (adjust cores and any required CLI flags)
- `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --cores 4`

Tips:
- Bind more paths as needed (e.g., scratch, reference volumes):
  - `apptainer exec -B /scratch:/scratch -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --cores 4`


## Phase 5 — (Optional) Use Snakemake’s Apptainer Integration Later

If/when you migrate to per-rule containers:
- Add `container: "docker://..."` or `"oras://..."` lines per rule
- Run on host (or inside your runner container) with Snakemake’s Apptainer integration:
  - `snakemake --use-apptainer --cores 4`
- Pass Apptainer args via:
  - `snakemake --use-apptainer --apptainer-args "-B $PWD:/work -W /work --containall"`

For now, the simplest approach is “run Snakemake inside the monolithic container,” as shown above.


## Phase 6 — Produce an x86_64 SIF for HPC via CI (if your cluster is x86)

If your deployment target is an x86_64 cluster, produce the x86_64 SIF in CI on Linux/AMD64 (fast and native), even if you develop on ARM64 locally.

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

2. Download the artifact and copy to your cluster.

3. Verify and run on the cluster:
- `apptainer exec pipeline_amd64.sif snakemake --version`
- `apptainer exec -B "$PWD":/work -W /work pipeline_amd64.sif snakemake --cores 8`

Signing (recommended for provenance)
- `apptainer key newpair`
- `apptainer sign containers/pipeline_amd64.sif`
- `apptainer verify containers/pipeline_amd64.sif`


## Phase 7 — Data Binding, Permissions, and Reproducibility

- Bind working directories explicitly:
  - `-B "$PWD":/work -W /work`
- Keep host data read-only where appropriate:
  - `-B "$PWD":/work:ro` (if you only need reads)
- Use containment for cleanliness:
  - `--containall --no-home`
- Record exact Snakemake and environment versions (your pipeline already writes `versions.txt`; keep doing that).
- Avoid writing to `/` in the container; use `-W /work` and write under the repo’s `results/` folder.


## Phase 8 — ARM64 Caveats and Fallbacks

- Package availability: Some Bioconda packages lack linux-aarch64 builds. If micromamba fails to solve for ARM64, you have two options:
  1) Swap/Pin packages to ARM-friendly alternatives (e.g., use bwa-mem2 vs bwa if supported, etc.).
  2) Build the SIF on x86_64 in CI (Phase 6) and run on x86_64 targets (HPC). For local ARM testing, you can still use the monolithic container but limit tool usage to what resolves.
- Performance: Native ARM64 VM (Lima) is efficient. Avoid cross-arch emulation locally for routine dev; use CI for x86 builds.
- File paths: On Lima, your macOS home is mounted under `/Users/<you>`. Prefer syncing the repo into `~/micro-wgs` in the VM for consistent Linux paths.


## Phase 9 — Security & Best Practices

- Keep `environment.yml` pinned (channels: conda-forge, bioconda).
- Use Apptainer in unprivileged mode (default) for safer execution.
- Sign SIF images and verify on target systems.
- Avoid embedding secrets in the image; mount credentials at runtime if needed (e.g., read-only tokens).
- Use read-only binds for inputs; write outputs to a dedicated `results/` directory.


## Phase 10 — Quickstart Summary

- Start VM and enter it:
  - `limactl start --name=apptainer-ubuntu template://ubuntu-lts`
  - `limactl shell apptainer-ubuntu`
- Install Apptainer:
  - `sudo apt-get update && sudo apt-get install -y apptainer squashfs-tools`
- Build SIF (in `~/micro-wgs`):
  - `apptainer build containers/pipeline_arm64.sif containers/pipeline.def`
- Dry-run pipeline:
  - `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake -n --cores 4`
- Run pipeline:
  - `apptainer exec -B "$PWD":/work -W /work containers/pipeline_arm64.sif snakemake --cores 4`
- For x86_64 HPC: use GitHub Actions to build `pipeline_amd64.sif` and run it on the cluster.

This plan keeps things practical: you develop and test locally in an ARM64 VM, package everything in a single Apptainer image from your `environment.yml`, and, if needed, produce an x86_64 SIF in CI for deployment to typical HPC environments.
