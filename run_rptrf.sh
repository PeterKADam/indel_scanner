#!/bin/bash
set -euo pipefail

# Example:
# ./run_rptrf.sh \
#   --fasta "/path/to/reference.fa" \
#   --sample "ph" \
#   --data-dir "/home/peterkad/pkadmaster/data" \
#   --rptrf "/home/peterkad/RPTRF/RPTRF2" \
#   --max-motif 100 \
#   --min-length 10

FASTA_FILE=""
SAMPLE_NAME=""
DATA_DIR=""
RPTRF_BIN="/home/peterkad/RPTRF/RPTRF2"
MAX_MOTIF=100
MIN_LENGTH=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta)
      FASTA_FILE="$2"
      shift 2
      ;;
    --sample)
      SAMPLE_NAME="$2"
      shift 2
      ;;
    --data-dir)
      DATA_DIR="$2"
      shift 2
      ;;
    --rptrf)
      RPTRF_BIN="$2"
      shift 2
      ;;
    --max-motif)
      MAX_MOTIF="$2"
      shift 2
      ;;
    --min-length)
      MIN_LENGTH="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$FASTA_FILE" || -z "$SAMPLE_NAME" || -z "$DATA_DIR" ]]; then
  echo "Usage: $0 --fasta FASTA --sample NAME --data-dir DIR [--rptrf PATH] [--max-motif N] [--min-length N]" >&2
  exit 1
fi

OUT_DIR="${DATA_DIR}/${SAMPLE_NAME}/repeatregions"
mkdir -p "$OUT_DIR"

if [[ -d "$RPTRF_BIN" ]]; then
  RPTRF_BIN="${RPTRF_BIN%/}/RPTRF"
fi

if [[ ! -x "$RPTRF_BIN" ]]; then
  echo "ERROR: RPTRF binary not executable at $RPTRF_BIN" >&2
  exit 1
fi

(
  cd "$OUT_DIR"
  "$RPTRF_BIN" -s "$FASTA_FILE" -m "$MAX_MOTIF" -t "$MIN_LENGTH"
)
echo "RPTRF output written to $OUT_DIR"
