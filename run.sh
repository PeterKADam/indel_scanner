#!/bin/bash

# Activate your python environment if you use one
# source /path/to/your/venv/bin/activate

# Define input files and output location
BAM_FILE="/home/peterkad/pkadmaster/data/mutationalscanning_bam/ph/diploid_assembly/ph_plus_unmapped_diploid_v2.bam"
FASTA_FILE="/home/peterkad/pkadmaster/data/ph/ph_diploid.fa"
OUTPUT_FILE="/home/peterkad/pkadmaster/indel_scanner/results/test/indels.tsv"
CONFIG_FILE="config.yaml"



# Run the scanner
# This assumes you run it from the root of the 'indel_scanner' project directory
python -m indel_scanner.main scan\
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --output "$OUTPUT_FILE" \
    --config "$CONFIG_FILE"