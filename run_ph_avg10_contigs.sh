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

# Run pipeline on ~10 contigs closest to average contig length.
# This uses the existing test-run selection logic in configurator/parallel code.
python -m indel_scanner.main \
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --strdir "$STR_DIR" \
    --output "$OUTPUT_FILE" \
    --config "$CONFIG_FILE" \
    --test-run \
    --test-contigs-count 10
