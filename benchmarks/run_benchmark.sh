#!/usr/bin/env bash
# One-command benchmark for skimflow: public data in, PASS/FAIL out.
#
#   1. Fetches the five public annelid skims of assets/test_manifest.tsv into
#      test_data/ (skipped for files that are already there).
#   2. Runs the pipeline on benchmarks/samplesheet_benchmark.csv. No sample has
#      a mito_seed, so every mitogenome comes from the default self-seeding.
#   3. Compares each mitogenome with its published GenBank record
#      (benchmarks/validate_mitogenomes.py, blastn in a container).
#   4. Checks the numbers against benchmarks/expected.tsv and the species trees
#      against benchmarks/expected_splits.tsv (benchmarks/check_benchmark.py).
#
# Exit code 0 means every check passed.
#
# Usage:
#   benchmarks/run_benchmark.sh [--engine podman|docker|apptainer] [--profile P]
#                               [--outdir DIR] [--skip-respect] [-- nextflow args]
#
# --engine picks the container runtime (default podman). --profile is the
# Nextflow profile list without the data profile (default: the engine name;
# on GWDG use --engine apptainer --profile gwdg inside an sbatch job).
# RESPECT is skipped automatically when no Gurobi licence is found
# (GUROBI_LIC or ~/gurobi.lic); the mitogenome checks do not need it.
# Needs: nextflow, the container engine, python3, curl and internet access
# (ENA for reads, NCBI for the reference records).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
ENGINE=podman
PROFILE=""
OUTDIR="$REPO/results_benchmark"
SKIP_RESPECT=0
BLAST_IMG=quay.io/biocontainers/blast:2.17.0--hb02a186_1

usage() { sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --engine) ENGINE="$2"; shift 2 ;;
        --profile) PROFILE="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --skip-respect) SKIP_RESPECT=1; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        *) echo "[benchmark] unknown option: $1 (see --help)" >&2; exit 2 ;;
    esac
done
case "$ENGINE" in podman|docker|apptainer) ;; *)
    echo "[benchmark] --engine must be podman, docker or apptainer" >&2; exit 2 ;;
esac
PROFILE="${PROFILE:-$ENGINE}"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
cd "$REPO"

echo "[benchmark] 1/4 test data"
# Fetch only the samples whose reads are missing, so files that are already
# there are never downloaded again.
missing=0
while IFS=$'\t' read -r acc pairs sample _rest <&3; do
    [[ -z "$acc" || "$acc" == accession || "$acc" == \#* ]] && continue
    if [[ -s "test_data/${sample}_R1.fastq.gz" && -s "test_data/${sample}_R2.fastq.gz" ]]; then
        continue
    fi
    missing=1
    ./scripts/fetch_test_data.sh "$acc" "$pairs" "$sample"
done 3< assets/test_manifest.tsv
(( missing )) || echo "[benchmark]     all reads present in test_data/"

EXTRA=()
lic="${GUROBI_LIC:-$HOME/gurobi.lic}"
if (( SKIP_RESPECT )) || [[ ! -f "$lic" ]]; then
    (( SKIP_RESPECT )) || echo "[benchmark]     no Gurobi licence at $lic, running with --skip_respect"
    EXTRA+=(--skip_respect)
fi

echo "[benchmark] 2/4 pipeline (profile $PROFILE)"
nextflow run . -profile "$PROFILE" \
    --input "$REPO/benchmarks/samplesheet_benchmark.csv" \
    --outdir "$OUTDIR" "${EXTRA[@]}" "$@"

echo "[benchmark] 3/4 validation against GenBank"
case "$ENGINE" in
    podman|docker)
        BLASTN="$ENGINE run --rm -v $REPO:$REPO -v $OUTDIR:$OUTDIR -w $REPO $BLAST_IMG blastn" ;;
    apptainer)
        BLASTN="apptainer exec -B $REPO -B $OUTDIR docker://$BLAST_IMG blastn" ;;
esac
python3 benchmarks/validate_mitogenomes.py \
    --results "$OUTDIR" --refs benchmarks/validation_refs.tsv \
    --ref-dir "$OUTDIR/benchmark_refs" --blastn "$BLASTN" \
    --out "$OUTDIR/validation.tsv"

echo "[benchmark] 4/4 checks"
set +e
python3 benchmarks/check_benchmark.py \
    --validation "$OUTDIR/validation.tsv" --expected benchmarks/expected.tsv \
    --phylogeny "$OUTDIR/phylogeny" --splits benchmarks/expected_splits.tsv \
    | tee "$OUTDIR/benchmark_checks.txt"
status=${PIPESTATUS[0]}
set -e
if (( status == 0 )); then
    echo "[benchmark] PASS (details in $OUTDIR/benchmark_checks.txt)"
else
    echo "[benchmark] FAIL (details in $OUTDIR/benchmark_checks.txt)"
fi
exit "$status"
