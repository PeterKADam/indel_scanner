#!/bin/bash

# Activate your python environment if you use one
# source /path/to/your/venv/bin/activate

# Define input files and output location
BAM_FILE="/home/peterkad/pkadmaster/data/mutationalscanning_bam/ph/diploid_assembly/ph_plus_unmapped_diploid_v2.bam"
FASTA_FILE="/home/peterkad/pkadmaster/data/ph/ph_diploid.fa"
OUTPUT_FILE="/home/peterkad/pkadmaster/indel_scanner/results/test/ph/"
CONFIG_FILE="config.yaml"
STR_DIR="repeatregions"

# Find the largest contig by length (2nd column from idxstats output).
LARGEST_CONTIG=$(samtools idxstats "$BAM_FILE" | sort -k2,2nr | head -1 | cut -f1)

if [ -z "$LARGEST_CONTIG" ]; then
    echo "ERROR: Could not determine largest contig from $BAM_FILE"
    exit 1
fi

echo "Running on largest contig: $LARGEST_CONTIG"

# Run the unified pipeline
# This assumes you run it from the root of the 'indel_scanner' project directory
python -m indel_scanner.main \
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --strdir "$STR_DIR" \
    --output "$OUTPUT_FILE" \
    --config "$CONFIG_FILE" \
    --contigs "$LARGEST_CONTIG"
