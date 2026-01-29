#!/bin/bash

# Activate your python environment if you use one
# source /path/to/your/venv/bin/activate

# Define input files and output location
# BAM file (must be indexed with .bai)
BAM_FILE=""
# Reference FASTA (must be indexed with .fai)
FASTA_FILE=""
# Base results directory: output becomes {OUTPUT_FILE}/{sample}/{timestamp}/
OUTPUT_FILE="results/"
# STR results directory for the sample (RPTRF output per contig)
STR_DIR=""
# Pipeline configuration YAML
CONFIG_FILE="config.yaml"

# Run the unified pipeline
# This assumes you run it from the root of the 'indel_scanner' project directory
python -m indel_scanner.main \
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --strdir "$STR_DIR" \
    --output "$OUTPUT_FILE" \
    --config "$CONFIG_FILE"