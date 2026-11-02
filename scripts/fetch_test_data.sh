#!/usr/bin/env bash
# Subsample real public paired-end FASTQ pairs from ENA into test_data/.
#
# Usage:
#   fetch_test_data.sh [ACC] [N] [NAME]     one run (default: DRR148567, 50000 pairs, demo)
#   fetch_test_data.sh --manifest <tsv>     one run per manifest row
#
# Default accession: DRR148567 - Bostrychus sinensis (goby) low-cov WGS, ~214k pairs total,
# used by the smoke profile. The manifest (e.g. assets/test_manifest.tsv) is a TSV with a
# header and the columns accession, pairs, sample_id, species_id, expected_genome_size_bp;
# each row writes test_data/<sample_id>_R{1,2}.fastq.gz from the first <pairs> read pairs.

set -euo pipefail

OUT_DIR="$(cd "$(dirname "$0")/.." && pwd)/test_data"
mkdir -p "$OUT_DIR"

# fetch_pair ACC N_READS OUT_NAME
fetch_pair() {
    local ACC="$1" N_READS="$2" OUT_NAME="$3"
    local URLS R1_URL R2_URL LINES n1 n2
    local R1_OUT="$OUT_DIR/${OUT_NAME}_R1.fastq.gz"
    local R2_OUT="$OUT_DIR/${OUT_NAME}_R2.fastq.gz"

    echo "[fetch] querying ENA for ${ACC}"
    URLS=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp&format=tsv" \
           | tail -n1 | awk '{print $NF}' | tr ';' '\n')

    # Pick mates by suffix: runs with unpaired reads also list <acc>.fastq.gz.
    R1_URL=$(echo "$URLS" | grep -E '_1\.f(ast)?q\.gz$' | head -n1 || true)
    R2_URL=$(echo "$URLS" | grep -E '_2\.f(ast)?q\.gz$' | head -n1 || true)

    if [[ -z "${R1_URL:-}" || -z "${R2_URL:-}" ]]; then
        echo "[fetch] ENA did not return _1/_2 paired FASTQ URLs for ${ACC}" >&2
        echo "$URLS" >&2
        return 1
    fi

    # ENA fastq_ftp paths usually have no scheme; force https:// (ENA serves both).
    R1_URL="https://${R1_URL#*://}"
    R2_URL="https://${R2_URL#*://}"

    echo "[fetch] R1: $R1_URL"
    echo "[fetch] R2: $R2_URL"
    echo "[fetch] subsampling first ${N_READS} pairs"

    LINES=$((N_READS * 4))

    # `head` closes its input after enough lines, which sends SIGPIPE upstream.
    # Disable pipefail just for these two pipelines so curl/zcat exiting 141/23 is fine.
    # curl then prints a harmless "(23) Failure writing output"; a real download
    # failure is caught by the line-count check below.
    set +o pipefail
    curl -fsSL "$R1_URL" | zcat | head -n "$LINES" | gzip > "$R1_OUT"
    curl -fsSL "$R2_URL" | zcat | head -n "$LINES" | gzip > "$R2_OUT"
    set -o pipefail

    n1=$(zcat "$R1_OUT" | wc -l)
    n2=$(zcat "$R2_OUT" | wc -l)
    if (( n1 == 0 || n1 % 4 != 0 || n1 != n2 )); then
        echo "[fetch] download failed or truncated for ${ACC} (R1=${n1} R2=${n2} lines)" >&2
        rm -f "$R1_OUT" "$R2_OUT"
        return 1
    fi
    if (( n1 < LINES )); then
        echo "[fetch] note: ${ACC} has only $((n1 / 4)) pairs (< ${N_READS})" >&2
    fi

    echo "[fetch] wrote:"
    ls -lh "$R1_OUT" "$R2_OUT"
}

if [[ "${1:-}" == "--manifest" ]]; then
    MANIFEST="${2:?--manifest needs a TSV path}"
    [[ -f "$MANIFEST" ]] || { echo "[fetch] manifest not found: $MANIFEST" >&2; exit 1; }
    # Read rows on fd 3 so nothing inside the loop can consume the manifest.
    while IFS=$'\t' read -r acc pairs sample_id _rest <&3; do
        [[ -z "${acc:-}" || "$acc" == "accession" || "$acc" == \#* ]] && continue
        fetch_pair "$acc" "$pairs" "$sample_id"
    done 3< "$MANIFEST"
else
    fetch_pair "${1:-DRR148567}" "${2:-50000}" "${3:-demo}"
fi
