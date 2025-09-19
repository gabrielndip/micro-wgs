# Repository Guidelines

## Project Structure & Module Organization
- `data/` — Input datasets and references. Treat as read-only; do not commit large binaries or secrets. Example: `data/samples.tsv`.
- `envs/` — Environment specs (e.g., Conda YAML). Pin versions and document usage. Example: `envs/rnaseq.yml`.
- `scripts/` — Pipeline steps and utilities. Prefer small, composable scripts. Example names: `qc_fastqc.sh`, `align_star.py`.
- `results/` — Generated outputs. Never hand-edit; scripts should write here in reproducible subfolders.

## Build, Test, and Development Commands
- Create env: `mamba env create -f envs/<name>.yml` (or `conda env create ...`).
- Activate: `mamba activate <name>`.
- Run a script: `bash scripts/<step>.sh` or `python scripts/<tool>.py`.
- Lint shell: `shellcheck scripts/*.sh`.

## Coding Style & Naming Conventions
- Bash: start scripts with `set -euo pipefail`; 2-space indent; functions and files use `snake_case`.
- Python: prefer `black` (line length 88), `isort`, and `flake8`; use `argparse` for CLIs; `snake_case` for functions/vars.
- Filenames: lowercase, `snake_case` (or hyphenated for CLIs); no spaces. Inputs read from `data/`, outputs written to `results/`.
- Script I/O: accept input/output paths via flags/env vars; avoid hard-coded absolute paths.

## Testing Guidelines
- Add smoke tests for each script (e.g., `scripts/tests/`) that exercise `--help` and a minimal run.
- Shell: run `shellcheck` and test with a tiny fixture in `data/`.
- Python: if present, add `pytest` tests and keep unit tests close to code. Target critical steps and parsers.

## Commit & Pull Request Guidelines
- Use Conventional Commits: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`.
- Keep PRs small and focused. Include: purpose, commands to reproduce (with paths), expected outputs (paths under `results/`), and any env changes (`envs/` diffs).
- Link issues and note data requirements (sample sizes, required references) in the PR description.

## Security & Configuration Tips
- Do not commit PHI/PII or large raw data. Add large patterns to `.gitignore` under `data/` and `results/`.
- Pin tool versions in `envs/`; prefer reproducible environments (Conda/Mamba or containers).
- Validate inputs before processing (e.g., sample sheet schema, file existence) and fail fast with clear messages.

## Agent-Specific Instructions
- Scope: entire repository. Preserve directory layout.
- Do not modify contents of `data/` or `results/` manually; produce changes via scripts in `scripts/`.
- Prefer additive changes (new scripts/configs) over destructive edits; avoid removing user data.
