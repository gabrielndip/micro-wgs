#!/usr/bin/env bash
set -euo pipefail

# Minimal dry-run to verify rule graph and that default targets resolve
snakemake -n --use-conda --cores 2
