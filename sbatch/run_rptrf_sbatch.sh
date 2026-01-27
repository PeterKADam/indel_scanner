#!/bin/bash
#SBATCH --job-name=rptrf
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=02:00:00

set -euo pipefail

# Example:
# sbatch sbatch/run_rptrf_sbatch.sh \
#   --fasta "/path/to/reference.fa" \
#   --sample "ph" \
#   --data-dir "/home/peterkad/pkadmaster/data" \
#   --rptrf "/home/peterkad/RPTRF/RPTRF2" \
#   --max-motif 100 \
#   --min-length 10

FASTA_FILE="$2"
SAMPLE_NAME="$4"
DATA_DIR="$6"
RPTRF_BIN="$8"
MAX_MOTIF="${10}"
MIN_LENGTH="${12}"

OUT_DIR="${DATA_DIR}/${SAMPLE_NAME}/repeatregions"
mkdir -p "$OUT_DIR"

cd "$OUT_DIR"
"${RPTRF_BIN}" -s "${FASTA_FILE}" -m "${MAX_MOTIF}" -t "${MIN_LENGTH}"
