#!/usr/bin/env bash
# One-shot launcher for the genome-skim pipeline on GWDG HPC (Göttingen).
# Uses Apptainer (the only container runtime on the cluster). Idempotent —
# safe to re-run; existing workspace and .sif are reused.
#
# What it does:
#   1. Loads modules (apptainer, java 17, nextflow if available).
#   2. Falls back to installing Nextflow into ~/bin if there's no module.
#   3. Finds (or allocates) a ceph-ssd workspace.
#   4. Builds the RESPECT Apptainer image. SCC policy forbids container
#      builds on login nodes, so on scc-* partitions the build runs inside
#      an interactive srun job.
#   5. Checks the Gurobi licence file is present.
#   6. Runs `nextflow run . -profile gwdg` with the chosen samplesheet.
#
# Usage:
#   ./scripts/run_gwdg.sh [options] [-- extra nextflow args]
#
# Options:
#   -i, --input PATH           Samplesheet CSV. Default: assets/samplesheet_test.csv
#   -l, --gurobi-lic PATH      Path to Gurobi WLS licence. Default: ~/gurobi.lic
#   -p, --partition NAME       SLURM partition. Default: scc-cpu
#                              (NHR users: standard96 / standard96s)
#   -f, --filesystem NAME      Workspace filesystem. Default: ceph-ssd
#                              (use ceph-hdd for >90 day runs)
#   -w, --workspace-name NAME  Workspace name. Default: genome-skim
#   -d, --workspace-days N     Initial allocation (days). Default: 30
#   -a, --account ID           Slurm --account=<id>. Usually unnecessary on
#                              project-specific GWDG usernames; only set if
#                              you have multiple project accounts.
#   -h, --help                 Show this help and exit
#
# Anything after a literal `--` is forwarded verbatim to `nextflow run`,
# e.g. `-resume`, `--outdir foo`, `-with-trace`.

set -euo pipefail

# --- Defaults ------------------------------------------------------------
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="${REPO_DIR}/assets/samplesheet_test.csv"
ACCOUNT="${GWDG_ACCOUNT:-}"
GUROBI_LIC="${GUROBI_LIC:-$HOME/gurobi.lic}"
PARTITION="${GWDG_PARTITION:-scc-cpu}"
WS_FS="${GWDG_WS_FS:-ceph-ssd}"
WS_NAME="${GWDG_WS_NAME:-genome-skim}"
WS_DAYS="${GWDG_WS_DAYS:-30}"
NF_EXTRA=()

log() { printf '[run_gwdg] %s\n' "$*"; }
die() { printf '[run_gwdg] ERROR: %s\n' "$*" >&2; exit 1; }
usage() { sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'; }

# --- Parse flags ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)            INPUT="$2";       shift 2 ;;
        -a|--account)          ACCOUNT="$2";     shift 2 ;;
        -l|--gurobi-lic)       GUROBI_LIC="$2";  shift 2 ;;
        -p|--partition)        PARTITION="$2";   shift 2 ;;
        -f|--filesystem)       WS_FS="$2";       shift 2 ;;
        -w|--workspace-name)   WS_NAME="$2";     shift 2 ;;
        -d|--workspace-days)   WS_DAYS="$2";     shift 2 ;;
        -h|--help)             usage; exit 0 ;;
        --)                    shift; NF_EXTRA=("$@"); break ;;
        -*)                    die "unknown option: $1 (try --help)" ;;
        *)                     die "unexpected positional arg: $1 (use --input PATH)" ;;
    esac
done

# --- Pre-flight ----------------------------------------------------------
[[ -f "$INPUT" ]]      || die "samplesheet not found: $INPUT"
[[ -f "$GUROBI_LIC" ]] || die "Gurobi licence not found at: $GUROBI_LIC
       Get a free academic WLS licence at https://www.gurobi.com/academia/
       and save it to ~/gurobi.lic (or pass --gurobi-lic /path/to/file)."

# --- Module loads --------------------------------------------------------
if ! type module >/dev/null 2>&1; then
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
type module >/dev/null 2>&1 || die "'module' is not available — are you on a GWDG login/compute node?"

log "loading modules"
# openjdk 17 — Nextflow ≥24.x needs Java 17. Bare `module load openjdk` may
# silently pick Java 11; pin to a 17.x version and fall back if absent.
module load openjdk/17.0.11_9 >/dev/null 2>&1 \
    || module load openjdk/17.0.8.1_1 >/dev/null 2>&1 \
    || module load openjdk >/dev/null 2>&1 \
    || die "no openjdk module available"
module load apptainer >/dev/null 2>&1 || die "failed to load apptainer module"

# Nextflow: prefer the cluster module, fall back to ~/bin install.
if ! module load nextflow/24.10.0 >/dev/null 2>&1 \
    && ! module load nextflow >/dev/null 2>&1; then
    mkdir -p "$HOME/bin"
    export PATH="$HOME/bin:$PATH"
    if ! command -v nextflow >/dev/null 2>&1; then
        log "no nextflow module available, installing into ~/bin"
        ( cd "$HOME/bin" && curl -fsSL https://get.nextflow.io | bash )
    fi
fi
command -v nextflow >/dev/null 2>&1 || die "nextflow is not on PATH after module/install attempts"

# --- Workspace -----------------------------------------------------------
if ws_list 2>/dev/null | awk '{print $1}' | grep -qx "$WS_NAME"; then
    GWDG_WORKSPACE="$(ws_find "$WS_NAME")"
    log "reusing workspace: $GWDG_WORKSPACE"
else
    log "allocating $WS_FS workspace '$WS_NAME' for $WS_DAYS days"
    ws_allocate -F "$WS_FS" "$WS_NAME" "$WS_DAYS" >/dev/null
    GWDG_WORKSPACE="$(ws_find "$WS_NAME")"
    log "allocated: $GWDG_WORKSPACE"
fi

APPTAINER_CACHEDIR="$GWDG_WORKSPACE/apptainer_cache"
mkdir -p "$APPTAINER_CACHEDIR"

# Env vars consumed by conf/gwdg.config.
export GWDG_WORKSPACE GWDG_PARTITION="$PARTITION"
[[ -n "$ACCOUNT" ]] && export GWDG_ACCOUNT="$ACCOUNT"
export APPTAINER_CACHEDIR GUROBI_LIC
export NXF_APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR"

# --- RESPECT container ---------------------------------------------------
RESPECT_SIF="$APPTAINER_CACHEDIR/respect_0.2.sif"
build_respect() {
    apptainer build "$RESPECT_SIF" "$REPO_DIR/containers/respect/respect.def"
}
if [[ ! -f "$RESPECT_SIF" ]]; then
    if [[ "$PARTITION" == scc-* ]]; then
        # SCC policy: don't build containers on login nodes.
        log "building RESPECT image inside an srun job (SCC policy)"
        srun --partition="$PARTITION" --time=00:30:00 \
             --cpus-per-task=4 --mem=8G \
             bash -c "$(declare -f build_respect); build_respect"
    else
        log "building RESPECT image (NHR — direct build on login node OK)"
        build_respect
    fi
else
    log "RESPECT image cached: $RESPECT_SIF"
fi

# --- Launch --------------------------------------------------------------
cd "$REPO_DIR"
log "launching pipeline"
log "  samplesheet : $INPUT"
log "  partition   : $PARTITION"
log "  workspace   : $GWDG_WORKSPACE ($WS_FS)"
log "  gurobi lic  : $GUROBI_LIC"
[[ -n "$ACCOUNT" ]]           && log "  account     : $ACCOUNT (explicit override)"
[[ ${#NF_EXTRA[@]} -gt 0 ]]   && log "  extra nf args: ${NF_EXTRA[*]}"
echo

nextflow run . \
    -profile gwdg \
    --input "$INPUT" \
    --respect_container "$RESPECT_SIF" \
    "${NF_EXTRA[@]}"

echo
log "done."
log "archive results:  rsync -av results/ \"\$PROJECT/genome-skim-results/\""
log "free workspace:   ws_release $WS_NAME"
