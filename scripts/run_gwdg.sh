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
#   5. Checks the Gurobi licence file when short-read RESPECT is needed.
#   6. Runs `nextflow run . -profile gwdg` with the chosen input.
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
#   --skip-respect             Do not run RESPECT; no Gurobi licence is needed.
#   --kraken2-db PATH          kraken2 DB directory for short-read decontam (optional).
#   --mitos2-refdir PATH       Pre-downloaded MITOS2 reference data directory (optional).
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
KRAKEN2_DB="${KRAKEN2_DB:-}"
MITOS2_REFDIR="${MITOS2_REFDIR:-}"
ACCOUNT="${GWDG_ACCOUNT:-}"
GUROBI_LIC="${GUROBI_LIC:-$HOME/gurobi.lic}"
SKIP_RESPECT="${SKIP_RESPECT:-false}"
PARTITION="${GWDG_PARTITION:-scc-cpu}"
if [[ -v GWDG_CLUSTER_OPTS ]]; then
    GWDG_CLUSTER_OPTS="${GWDG_CLUSTER_OPTS}"
else
    GWDG_CLUSTER_OPTS="--constraint=inet"
fi
GWDG_PROXY="${GWDG_PROXY:-http://www-cache.gwdg.de:3128}"
WS_FS="${GWDG_WS_FS:-ceph-ssd}"
WS_NAME="${GWDG_WS_NAME:-genome-skim}"
WS_DAYS="${GWDG_WS_DAYS:-30}"
NF_EXTRA=()

log() { printf '[run_gwdg] %s\n' "$*"; }
die() { printf '[run_gwdg] ERROR: %s\n' "$*" >&2; exit 1; }
abs_path() {
    local p="$1"
    [[ -z "$p" ]] && return 0
    realpath "$p"
}
samplesheet_has_short_reads() {
    local csv="$1"
    awk -F',' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
                gsub(/^"|"$/, "", $i)
                if ($i == "fastq_1") col = i
            }
            next
        }
        col {
            v = $col
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
            gsub(/^"|"$/, "", v)
            if (v != "") { found = 1; exit }
        }
        END { exit found ? 0 : 1 }
    ' "$csv"
}
nf_extra_skips_respect() {
    local prev=""
    local arg
    for arg in "${NF_EXTRA[@]}"; do
        if [[ "$prev" == "--skip_respect" && "$arg" == "true" ]]; then
            return 0
        fi
        if [[ "$arg" == "--skip_respect=true" ]]; then
            return 0
        fi
        prev="$arg"
    done
    return 1
}
usage() {
    awk '
        NR == 1 { next }
        /^#/ { sub(/^# ?/, ""); print; next }
        { exit }
    ' "$0"
}

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
        --skip-respect)        SKIP_RESPECT=true; shift ;;
        --kraken2-db)          KRAKEN2_DB="$2";  shift 2 ;;
        --mitos2-refdir)       MITOS2_REFDIR="$2"; shift 2 ;;
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
[[ -z "$MITOS2_REFDIR" || -d "$MITOS2_REFDIR" ]] || die "MITOS2 refdir not found: $MITOS2_REFDIR"

HAS_SHORT_SAMPLES=false
if [[ -n "$R1" ]]; then
    HAS_SHORT_SAMPLES=true
elif [[ -n "$INPUT" ]] && samplesheet_has_short_reads "$INPUT"; then
    HAS_SHORT_SAMPLES=true
fi
if nf_extra_skips_respect; then
    SKIP_RESPECT=true
fi

if [[ "$HAS_SHORT_SAMPLES" == true && "$SKIP_RESPECT" != true ]]; then
    [[ -f "$GUROBI_LIC" ]] || die "Gurobi licence not found at: $GUROBI_LIC
       Get a free academic WLS licence at https://www.gurobi.com/academia/
       and save it to ~/gurobi.lic (or pass --gurobi-lic /path/to/file)."
fi

[[ -n "$R1" ]]            && R1="$(abs_path "$R1")"
[[ -n "$R2" ]]            && R2="$(abs_path "$R2")"
[[ -n "$LONG_READS" ]]    && LONG_READS="$(abs_path "$LONG_READS")"
[[ -n "$INPUT" ]]         && INPUT="$(abs_path "$INPUT")"
[[ -n "$KRAKEN2_DB" ]]    && KRAKEN2_DB="$(abs_path "$KRAKEN2_DB")"
[[ -n "$MITOS2_REFDIR" ]] && MITOS2_REFDIR="$(abs_path "$MITOS2_REFDIR")"
[[ -f "$GUROBI_LIC" ]]    && GUROBI_LIC="$(abs_path "$GUROBI_LIC")"

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
export GWDG_WORKSPACE GWDG_PARTITION="$PARTITION" GWDG_CLUSTER_OPTS
[[ -n "$ACCOUNT" ]] && export GWDG_ACCOUNT="$ACCOUNT"
export APPTAINER_CACHEDIR GUROBI_LIC
export NXF_APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR"
if [[ -n "$GWDG_PROXY" ]]; then
    export http_proxy="${http_proxy:-$GWDG_PROXY}"
    export https_proxy="${https_proxy:-$GWDG_PROXY}"
    export ftp_proxy="${ftp_proxy:-$GWDG_PROXY}"
    export APPTAINERENV_http_proxy="${APPTAINERENV_http_proxy:-$GWDG_PROXY}"
    export APPTAINERENV_https_proxy="${APPTAINERENV_https_proxy:-$GWDG_PROXY}"
    export APPTAINERENV_ftp_proxy="${APPTAINERENV_ftp_proxy:-$GWDG_PROXY}"
fi

# --- RESPECT container ---------------------------------------------------
RESPECT_SIF="$APPTAINER_CACHEDIR/respect_0.2.sif"
export REPO_DIR RESPECT_SIF
build_respect() {
    apptainer build "$RESPECT_SIF" "$REPO_DIR/containers/respect/respect.def"
}
if [[ "$HAS_SHORT_SAMPLES" != true ]]; then
    log "no short-read samples: skipping RESPECT image build (RESPECT not used)"
elif [[ "$SKIP_RESPECT" == true ]]; then
    log "skipping RESPECT image build (--skip-respect)"
elif [[ ! -f "$RESPECT_SIF" ]]; then
    if [[ "$PARTITION" == scc-* ]]; then
        # SCC policy: don't build containers on login nodes.
        log "building RESPECT image inside an srun job (SCC policy)"
        SRUN_EXTRA=()
        [[ -n "$GWDG_CLUSTER_OPTS" ]] && read -r -a SRUN_EXTRA <<< "$GWDG_CLUSTER_OPTS"
        srun "${SRUN_EXTRA[@]}" --partition="$PARTITION" --time=00:30:00 \
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
[[ -n "$MITOS2_REFDIR" ]]     && log "  MITOS2 ref  : $MITOS2_REFDIR"
[[ "$SKIP_RESPECT" == true ]] && log "  RESPECT     : skipped"
log "  partition   : $PARTITION"
log "  workspace   : $GWDG_WORKSPACE ($WS_FS)"
[[ "$SKIP_RESPECT" != true ]]  && log "  gurobi lic  : $GUROBI_LIC"
[[ -n "$ACCOUNT" ]]           && log "  account     : $ACCOUNT (explicit override)"
[[ ${#NF_EXTRA[@]} -gt 0 ]]   && log "  extra nf args: ${NF_EXTRA[*]}"
echo

NF_ARGS=( -profile gwdg --respect_container "$RESPECT_SIF" )
[[ "$SKIP_RESPECT" == true ]] && NF_ARGS+=( --skip_respect true )
[[ -n "$KRAKEN2_DB" ]] && NF_ARGS+=( --kraken2_db "$KRAKEN2_DB" )
[[ -n "$MITOS2_REFDIR" ]] && NF_ARGS+=( --mitos2_refdir "$MITOS2_REFDIR" )
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
