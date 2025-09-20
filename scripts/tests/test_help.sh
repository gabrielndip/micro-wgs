#!/usr/bin/env bash
set -euo pipefail

echo "Snakemake version:" && snakemake --version
echo "Listing rules:" && snakemake -n --cores 1 --list

