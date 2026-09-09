#!/usr/bin/env bash
# Coverage-titration benchmark for skimflow (standalone; not a Nextflow entry).
#
# Down-samples one high-coverage sample to a ladder of target coverages (several
# random seeds each), runs the full pipeline at every level, and collects each
# level's per-sample metrics row into one tidy titration table. Use it to find
# the minimum skim depth for reliable mitogenome / gene recovery (the classic
# genome-skimming experiment; see docs and the thesis benchmarking plan).
#
# Down-sampling uses rasusa (coverage-aware, seedable, paired-aware) run inside
# the same container engine the pipeline uses, so no host tools are required
# beyond the engine and nextflow.
#
# Example:
#   benchmarks/run_titration.sh \
#       --r1 test_data/pmisa_small_R1.fastq.gz \
#       --r2 test_data/pmisa_small_R2.fastq.gz \
#       --genome-size 1200000000 \
#       --coverages 0.05,0.1,0.25,0.5,1,2,5 \
#       --seeds 3 --engine podman \
#       --organelle-db results/mitogenome/getorganelle_db \
#       --busco-db results/markers/busco_downloads \
#       -- -resume
#
# Everything after `--` is forwarded verbatim to `nextflow run`.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

R1=""; R2=""; GENOME_SIZE=""
COVERAGES="0.05,0.1,0.25,0.5,1,2,5"
SEEDS=1
ENGINE="podman"
SAMPLE_ID="bench"
OUTDIR="$REPO_DIR/benchmarks/titration_out"
ORGANELLE_DB=""; BUSCO_DB=""
RASUSA_IMG="quay.io/biocontainers/rasusa:5.1.0--hfa8f182_0"
NF_PASSTHROUGH=()

usage() {
    sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1)            R1="$2"; shift 2 ;;
        --r2)            R2="$2"; shift 2 ;;
        --genome-size)   GENOME_SIZE="$2"; shift 2 ;;
        --coverages)     COVERAGES="$2"; shift 2 ;;
        --seeds)         SEEDS="$2"; shift 2 ;;
        --engine)        ENGINE="$2"; shift 2 ;;
        --sample-id)     SAMPLE_ID="$2"; shift 2 ;;
        --outdir)        OUTDIR="$2"; shift 2 ;;
        --organelle-db)  ORGANELLE_DB="$2"; shift 2 ;;
        --busco-db)      BUSCO_DB="$2"; shift 2 ;;
        --rasusa-image)  RASUSA_IMG="$2"; shift 2 ;;
        -h|--help)       usage 0 ;;
        --)              shift; NF_PASSTHROUGH=("$@"); break ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

[[ -n "$R1" ]]           || { echo "ERROR: --r1 is required" >&2; usage 1; }
[[ -n "$GENOME_SIZE" ]]  || { echo "ERROR: --genome-size (bp) is required for coverage titration" >&2; usage 1; }
R1="$(cd "$(dirname "$R1")" && pwd)/$(basename "$R1")"
[[ -f "$R1" ]] || { echo "ERROR: R1 not found: $R1" >&2; exit 1; }
if [[ -n "$R2" ]]; then
    R2="$(cd "$(dirname "$R2")" && pwd)/$(basename "$R2")"
    [[ -f "$R2" ]] || { echo "ERROR: R2 not found: $R2" >&2; exit 1; }
fi

mkdir -p "$OUTDIR"
TITRATION="$OUTDIR/titration.csv"
: > "$TITRATION"
header_written=0

# Run a command inside the chosen engine, binding the dirs it needs.
# Usage: run_in_engine <dir>... -- <command>...
run_in_engine() {
    local -a binds=()
    local d
    while [[ $# -gt 0 && "$1" != "--" ]]; do
        binds+=("$1"); shift
    done
    shift  # drop the "--" separator; "$@" is now the command
    case "$ENGINE" in
        podman|docker)
            local -a flags=()
            for d in "${binds[@]}"; do flags+=(-v "$d:$d"); done
            "$ENGINE" run --rm "${flags[@]}" "$RASUSA_IMG" "$@" ;;
        apptainer|singularity)
            local -a flags=()
            for d in "${binds[@]}"; do flags+=(-B "$d"); done
            "$ENGINE" exec "${flags[@]}" "docker://$RASUSA_IMG" "$@" ;;
        *) echo "ERROR: unsupported --engine '$ENGINE'" >&2; exit 1 ;;
    esac
}

IFS=',' read -r -a COV_LIST <<< "$COVERAGES"

for cov in "${COV_LIST[@]}"; do
    for seed in $(seq 1 "$SEEDS"); do
        run_id="${SAMPLE_ID}_c${cov}_s${seed}"
        run_dir="$OUTDIR/$run_id"
        mkdir -p "$run_dir"
        sub_r1="$run_dir/${run_id}_R1.fastq.gz"
        echo ">>> [$run_id] down-sampling to ${cov}x (seed $seed)"

        # Directories rasusa must see: the reads' dir(s) and the run dir.
        bind_dirs=("$(dirname "$R1")" "$run_dir")
        [[ -n "$R2" ]] && bind_dirs+=("$(dirname "$R2")")

        if [[ -n "$R2" ]]; then
            sub_r2="$run_dir/${run_id}_R2.fastq.gz"
            run_in_engine "${bind_dirs[@]}" -- \
                rasusa reads --coverage "$cov" --genome-size "$GENOME_SIZE" \
                    --seed "$seed" -o "$sub_r1" -o "$sub_r2" "$R1" "$R2"
        else
            run_in_engine "${bind_dirs[@]}" -- \
                rasusa reads --coverage "$cov" --genome-size "$GENOME_SIZE" \
                    --seed "$seed" -o "$sub_r1" "$R1"
        fi

        echo ">>> [$run_id] running pipeline"
        nf_args=(run "$REPO_DIR" -profile "$ENGINE"
                 --r1 "$sub_r1"
                 --sample_id "$run_id"
                 --expected_size "$GENOME_SIZE"
                 --outdir "$run_dir/results")
        [[ -n "$R2" ]]           && nf_args+=(--r2 "$sub_r2")
        [[ -n "$ORGANELLE_DB" ]] && nf_args+=(--organelle_db "$ORGANELLE_DB")
        [[ -n "$BUSCO_DB" ]]     && nf_args+=(--busco_db "$BUSCO_DB")
        [[ ${#NF_PASSTHROUGH[@]} -gt 0 ]] && nf_args+=("${NF_PASSTHROUGH[@]}")
        nextflow "${nf_args[@]}"

        metrics="$run_dir/results/summary/skimflow_metrics.csv"
        if [[ -f "$metrics" ]]; then
            if [[ "$header_written" -eq 0 ]]; then
                printf 'coverage,seed,%s\n' "$(head -1 "$metrics")" >> "$TITRATION"
                header_written=1
            fi
            tail -n +2 "$metrics" | while IFS= read -r row; do
                [[ -n "$row" ]] && printf '%s,%s,%s\n' "$cov" "$seed" "$row" >> "$TITRATION"
            done
        else
            echo "WARN: no metrics for $run_id (pipeline produced no summary)" >&2
        fi
    done
done

echo
echo "Titration table: $TITRATION"
[[ -s "$TITRATION" ]] && column -t -s, "$TITRATION" || echo "(no rows collected)"
