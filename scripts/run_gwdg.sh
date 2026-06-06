#!/usr/bin/env bash
# One-shot launcher for the genome-skim pipeline on GWDG HPC (Göttingen).
# Uses Apptainer (the only container runtime on the cluster). Idempotent -
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
# Input options (pick ONE):
#   --r1 PATH                  R1 FASTQ for a single-sample run.
#   --r2 PATH                  R2 FASTQ (omit for single-end).
#   --sample-id NAME           Override sample_id (default: derived from R1).
#   --species NAME             Override species name (default: sample_id).
#   --expected-size BP         Optional expected genome size (bp), shown in the report.
#   --long-reads PATH          Long-read FASTQ for a single long-read sample.
#                              Long-only direct mode; do not combine with --r1/--input.
#   --lr-type NAME             Long-read tech: nanopore (default) or pacbio.
#   -i, --input PATH           Samplesheet CSV (multi-sample alternative).
#                              Default if neither --r1 nor --input is set:
#                              assets/samplesheet_test.csv
#
# Run options:
#   -l, --gurobi-lic PATH      Path to Gurobi WLS licence. Default: ~/gurobi.lic
#   --kraken2-db PATH          kraken2 DB directory for short-read decontam (optional).
#                              Known-good GWDG default (Mateo): /mnt/ceph-ssd/workspaces/
#                              ws/scc_ubet_bleidorn/u16307-genome_size/dodecaceria/
#                              k2_pluspfp_16_GB_20260226
#   -p, --partition NAME       SLURM partition. Default: scc-cpu
#                              (NHR users: standard96 / standard96s)
#   -f, --filesystem NAME      Workspace filesystem. Default: ceph-ssd
#                              (use ceph-hdd if you need >90 days total -
#                              ceph-ssd hard caps at 30d × 2 extensions = 90d)
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
INPUT=""
R1=""
R2=""
SAMPLE_ID=""
SPECIES=""
EXPECTED_SIZE=""
LONG_READS=""
LR_TYPE=""
KRAKEN2_DB=""
ACCOUNT="${GWDG_ACCOUNT:-}"
GUROBI_LIC="${GUROBI_LIC:-$HOME/gurobi.lic}"
PARTITION="${GWDG_PARTITION:-scc-cpu}"
WS_FS="${GWDG_WS_FS:-ceph-ssd}"
WS_NAME="${GWDG_WS_NAME:-genome-skim}"
WS_DAYS="${GWDG_WS_DAYS:-30}"
NF_EXTRA=()

log() { printf '[run_gwdg] %s\n' "$*"; }
die() { printf '[run_gwdg] ERROR: %s\n' "$*" >&2; exit 1; }
usage() { sed -n '2,51p' "$0" | sed 's/^# \{0,1\}//'; }

# --- Parse flags ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1)                  R1="$2";          shift 2 ;;
        --r2)                  R2="$2";          shift 2 ;;
        --sample-id)           SAMPLE_ID="$2";   shift 2 ;;
        --species)             SPECIES="$2";     shift 2 ;;
        --expected-size)       EXPECTED_SIZE="$2"; shift 2 ;;
        --long-reads)          LONG_READS="$2";  shift 2 ;;
        --lr-type)             LR_TYPE="$2";     shift 2 ;;
        --kraken2-db)          KRAKEN2_DB="$2";  shift 2 ;;
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
        *)                     die "unexpected positional arg: $1 (use --r1 or --input)" ;;
    esac
done

# --- Pre-flight ----------------------------------------------------------
# Exactly one input mode: --r1 (short), --long-reads (long), or a samplesheet.
# The committed test CSV is the implicit default when none is given.
[[ -n "$LONG_READS" && ( -n "$R1" || -n "$INPUT" ) ]] && \
    die "--long-reads is long-only; do not combine with --r1 or --input"

if [[ -n "$R1" ]]; then
    [[ -f "$R1" ]] || die "R1 not found: $R1"
    [[ -z "$R2" || -f "$R2" ]] || die "R2 not found: $R2"
elif [[ -n "$LONG_READS" ]]; then
    [[ -f "$LONG_READS" ]] || die "long reads not found: $LONG_READS"
elif [[ -z "$INPUT" ]]; then
    INPUT="${REPO_DIR}/assets/samplesheet_test.csv"
fi
[[ -z "$INPUT" || -f "$INPUT" ]] || die "samplesheet not found: $INPUT"
[[ -z "$KRAKEN2_DB" || -d "$KRAKEN2_DB" ]] || die "kraken2 DB dir not found: $KRAKEN2_DB"
if [[ -z "$LONG_READS" ]]; then
    [[ -f "$GUROBI_LIC" ]] || die "Gurobi licence not found at: $GUROBI_LIC
       Get a free academic WLS licence at https://www.gurobi.com/academia/
       and save it to ~/gurobi.lic (or pass --gurobi-lic /path/to/file)."
fi

# --- Module loads --------------------------------------------------------
if ! type module >/dev/null 2>&1; then
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
type module >/dev/null 2>&1 || die "'module' is not available - are you on a GWDG login/compute node?"

log "loading modules"
# openjdk 17 - Nextflow ≥24.x needs Java 17. We don't pin a patch version
# because GWDG bumps them and an exact pin breaks silently; Lmod resolves
# `openjdk/17` to the latest 17.x. Fall through to bare `openjdk` as a
# last resort and check the version after.
module load openjdk/17 >/dev/null 2>&1 \
    || module load openjdk >/dev/null 2>&1 \
    || die "no openjdk module available (run 'module spider openjdk' on a login node)"
module load apptainer >/dev/null 2>&1 || die "failed to load apptainer module"

# Nextflow: GWDG doesn't publish a documented `nextflow` module (no entry
# in the software stacks docs as of 2026-05). We probe for one anyway in
# case it shows up later, then fall through to the ~/bin curl install.
if ! module load nextflow >/dev/null 2>&1; then
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
if [[ -n "$LONG_READS" ]]; then
    log "long-read-only run: skipping RESPECT licence check and image build (RESPECT not used)"
elif [[ ! -f "$RESPECT_SIF" ]]; then
    if [[ "$PARTITION" == scc-* ]]; then
        # SCC policy: don't build containers on login nodes.
        log "building RESPECT image inside an srun job (SCC policy)"
        srun --partition="$PARTITION" --time=00:30:00 \
             --cpus-per-task=4 --mem=8G \
             bash -c "$(declare -f build_respect); build_respect"
    else
        log "building RESPECT image (NHR - direct build on login node OK)"
        build_respect
    fi
else
    log "RESPECT image cached: $RESPECT_SIF"
fi

# --- Launch --------------------------------------------------------------
cd "$REPO_DIR"
log "launching pipeline"
if [[ -n "$R1" ]]; then
    log "  reads       : $R1${R2:+ + $R2}"
    [[ -n "$SAMPLE_ID" ]]     && log "  sample_id   : $SAMPLE_ID"
    [[ -n "$SPECIES" ]]       && log "  species     : $SPECIES"
    [[ -n "$EXPECTED_SIZE" ]] && log "  expected bp : $EXPECTED_SIZE"
elif [[ -n "$LONG_READS" ]]; then
    log "  long reads  : $LONG_READS"
    log "  lr_type     : ${LR_TYPE:-nanopore (default)}"
    [[ -n "$SAMPLE_ID" ]]     && log "  sample_id   : $SAMPLE_ID"
    [[ -n "$SPECIES" ]]       && log "  species     : $SPECIES"
    [[ -n "$EXPECTED_SIZE" ]] && log "  expected bp : $EXPECTED_SIZE"
else
    log "  samplesheet : $INPUT"
fi
[[ -n "$KRAKEN2_DB" ]]        && log "  kraken2 db  : $KRAKEN2_DB"
log "  partition   : $PARTITION"
log "  workspace   : $GWDG_WORKSPACE ($WS_FS)"
log "  gurobi lic  : $GUROBI_LIC"
[[ -n "$ACCOUNT" ]]           && log "  account     : $ACCOUNT (explicit override)"
[[ ${#NF_EXTRA[@]} -gt 0 ]]   && log "  extra nf args: ${NF_EXTRA[*]}"
echo

NF_ARGS=( -profile gwdg --respect_container "$RESPECT_SIF" )
[[ -n "$KRAKEN2_DB" ]] && NF_ARGS+=( --kraken2_db "$KRAKEN2_DB" )
if [[ -n "$R1" ]]; then
    NF_ARGS+=( --r1 "$R1" )
    [[ -n "$R2" ]]            && NF_ARGS+=( --r2 "$R2" )
    [[ -n "$SAMPLE_ID" ]]     && NF_ARGS+=( --sample_id "$SAMPLE_ID" )
    [[ -n "$SPECIES" ]]       && NF_ARGS+=( --species "$SPECIES" )
    [[ -n "$EXPECTED_SIZE" ]] && NF_ARGS+=( --expected_size "$EXPECTED_SIZE" )
elif [[ -n "$LONG_READS" ]]; then
    NF_ARGS+=( --long_reads "$LONG_READS" )
    [[ -n "$LR_TYPE" ]]       && NF_ARGS+=( --lr_type "$LR_TYPE" )
    [[ -n "$SAMPLE_ID" ]]     && NF_ARGS+=( --sample_id "$SAMPLE_ID" )
    [[ -n "$SPECIES" ]]       && NF_ARGS+=( --species "$SPECIES" )
    [[ -n "$EXPECTED_SIZE" ]] && NF_ARGS+=( --expected_size "$EXPECTED_SIZE" )
else
    NF_ARGS+=( --input "$INPUT" )
fi

nextflow run . "${NF_ARGS[@]}" "${NF_EXTRA[@]}"

echo
log "done."
# Persistent project storage on GWDG resolves through ~/.project/<projectid>/
# rather than a guaranteed $PROJECT env var (see GWDG storage map docs).
log "archive results:  rsync -av results/ ~/.project/<projectid>/genome-skim-results/"
log "free workspace:   ws_release $WS_NAME"
