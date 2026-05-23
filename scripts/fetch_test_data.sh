#!/usr/bin/env bash
# Subsample a real public paired-end FASTQ pair from ENA into test_data/.
# Default accession: DRR148567 - Bostrychus sinensis (goby) low-cov WGS, ~214k pairs total.
# We take the first N pairs to keep the smoke test fast (~few MB on disk).

set -euo pipefail

ACC="${1:-DRR148567}"
N_READS="${2:-50000}"          # number of read PAIRS to keep
OUT_NAME="${3:-demo}"           # produces ${OUT_NAME}_R1.fastq.gz and _R2.fastq.gz

OUT_DIR="$(cd "$(dirname "$0")/.." && pwd)/test_data"
mkdir -p "$OUT_DIR"

echo "[fetch] querying ENA for ${ACC}"
URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp&format=tsv" \
       | tail -n1 | awk '{print $NF}' | tr ';' '\n')

R1_URL=$(echo "$URLS" | sed -n '1p')
R2_URL=$(echo "$URLS" | sed -n '2p')

if [[ -z "${R1_URL:-}" || -z "${R2_URL:-}" ]]; then
    echo "[fetch] ENA did not return two FASTQ URLs for ${ACC}" >&2
    echo "$URLS" >&2
    exit 1
fi

# ENA returns ftp:// paths without a scheme; prepend https:// (their server speaks both).
R1_URL="https://${R1_URL#ftp://}"
R2_URL="https://${R2_URL#ftp://}"
R1_URL="${R1_URL#https://https://}"; R1_URL="https://${R1_URL#https://}"
R2_URL="${R2_URL#https://https://}"; R2_URL="https://${R2_URL#https://}"

echo "[fetch] R1: $R1_URL"
echo "[fetch] R2: $R2_URL"
echo "[fetch] subsampling first ${N_READS} pairs"

LINES=$((N_READS * 4))

# `head` closes its input after enough lines, which sends SIGPIPE upstream.
# Disable pipefail just for these two pipelines so curl/zcat exiting 141/23 is fine.
set +o pipefail
curl -fsSL "$R1_URL" 2>/dev/null | zcat | head -n "$LINES" | gzip > "$OUT_DIR/${OUT_NAME}_R1.fastq.gz"
curl -fsSL "$R2_URL" 2>/dev/null | zcat | head -n "$LINES" | gzip > "$OUT_DIR/${OUT_NAME}_R2.fastq.gz"
set -o pipefail

echo "[fetch] wrote:"
ls -lh "$OUT_DIR/${OUT_NAME}_R1.fastq.gz" "$OUT_DIR/${OUT_NAME}_R2.fastq.gz"
