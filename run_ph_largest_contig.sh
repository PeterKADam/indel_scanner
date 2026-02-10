#!/bin/bash
set -euo pipefail

# Define input files and output location
# BAM file (must be indexed with .bai)
BAM_FILE="/home/peterkad/pkadmaster/data/mutationalscanning_bam/ph/diploid_assembly/ph_plus_unmapped_diploid_v2.bam"
# Reference FASTA (must be indexed with .fai)
FASTA_FILE="/home/peterkad/pkadmaster/data/ph/ph_diploid.fa"
# Base results directory: output becomes {OUTPUT_FILE}/{sample}/{timestamp}/
OUTPUT_FILE="/home/peterkad/pkadmaster/indel_scanner/results"
# STR results directory for the sample (RPTRF output per contig)
STR_DIR="/home/peterkad/pkadmaster/data/ph/repeatregions"
# Pipeline configuration YAML
CONFIG_FILE="config.yaml"

# Resolve the largest contig from BAM header
LARGEST_CONTIG="$(python - "$BAM_FILE" <<'PY'
import sys
import pysam

bam_path = sys.argv[1]
with pysam.AlignmentFile(bam_path, "rb") as bam:
    refs = list(zip(bam.header.references, bam.header.lengths))
if not refs:
    raise SystemExit("No contigs found in BAM header.")
print(max(refs, key=lambda x: x[1])[0])
PY
)"

echo "Largest contig in BAM: ${LARGEST_CONTIG}"

# Run the unified pipeline on only the largest contig
# This assumes you run it from the root of the 'indel_scanner' project directory
python -m indel_scanner.main \
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --strdir "$STR_DIR" \
    --output "$OUTPUT_FILE" \
    --config "$CONFIG_FILE" \
    --contigs "$LARGEST_CONTIG"
