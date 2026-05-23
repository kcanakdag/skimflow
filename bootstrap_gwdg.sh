#!/usr/bin/env bash
# Self-contained bootstrap for running the genome-skim pipeline on GWDG.
# Clones (or updates) the pipeline repo into ~/projects, then hands off to
# scripts/run_gwdg.sh which builds the container and launches the run.
#
# One-liner:
#   curl -fsSL https://raw.githubusercontent.com/kcanakdag/skimflow/main/bootstrap_gwdg.sh \
#       | bash -s -- --input my.csv
#
# Any flags after `bash -s --` are forwarded to scripts/run_gwdg.sh
# (see `--help` there for the full list).
#
# Env var overrides (all optional):
#   PIPELINE_REPO    Git URL of the pipeline. Default: kcanakdag/skimflow on GitHub.
#   PIPELINE_REF     Branch/tag/commit to check out. Default: main.
#   PIPELINE_DIR     Where to clone. Default: $HOME/projects/skimflow.

set -euo pipefail

DEFAULT_REPO="https://github.com/kcanakdag/skimflow.git"

PIPELINE_REPO="${PIPELINE_REPO:-$DEFAULT_REPO}"
PIPELINE_REF="${PIPELINE_REF:-main}"
PIPELINE_DIR="${PIPELINE_DIR:-$HOME/projects/skimflow}"

log() { printf '[bootstrap] %s\n' "$*"; }
die() { printf '[bootstrap] ERROR: %s\n' "$*" >&2; exit 1; }

command -v git >/dev/null 2>&1 || die "git is not on PATH"

mkdir -p "$(dirname "$PIPELINE_DIR")"

if [[ -d "$PIPELINE_DIR/.git" ]]; then
    log "updating existing checkout at $PIPELINE_DIR"
    git -C "$PIPELINE_DIR" fetch --quiet --tags origin
    git -C "$PIPELINE_DIR" checkout --quiet "$PIPELINE_REF"
    # Only fast-forward if we're tracking a branch (skip for detached tags/commits).
    if git -C "$PIPELINE_DIR" symbolic-ref -q HEAD >/dev/null; then
        git -C "$PIPELINE_DIR" pull --ff-only --quiet
    fi
else
    log "cloning $PIPELINE_REPO into $PIPELINE_DIR"
    git clone --quiet --branch "$PIPELINE_REF" "$PIPELINE_REPO" "$PIPELINE_DIR"
fi

log "handing off to $PIPELINE_DIR/scripts/run_gwdg.sh"
exec bash "$PIPELINE_DIR/scripts/run_gwdg.sh" "$@"
