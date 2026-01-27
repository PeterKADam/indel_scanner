#!/bin/bash

# Activate your python environment if you use one
# source /path/to/your/venv/bin/activate

# Define input files and output location
BAM_FILE="/home/peterkad/pkadmaster/data/mutationalscanning_bam/ph/diploid_assembly/ph_plus_unmapped_diploid_v2.bam"
FASTA_FILE="/home/peterkad/pkadmaster/data/ph/ph_diploid.fa"
# Base results directory: output becomes {OUTPUT_FILE}/{sample}/{timestamp}/
OUTPUT_FILE="/home/peterkad/pkadmaster/indel_scanner/results"
# STR results directory for the sample
STR_DIR="/home/peterkad/pkadmaster/data/ph/repeatregions"
CONFIG_FILE="config.yaml"

#BAM_FILE_tr="/home/peterkad/pkadmaster/data/mutationalscanning_bam/tr/diploid_assembly/tr_plus_unmapped_diploid_v2.bam"
#FASTA_FILE_tr="/home/peterkad/pkadmaster/data/tr/tr_diploid.fa"
#OUTPUT_FILE_tr="/home/peterkad/pkadmaster/indel_scanner/results/test"
#CONFIG_FILE="config.yaml"


# Run the unified pipeline
# This assumes you run it from the root of the 'indel_scanner' project directory
python -m indel_scanner.main \
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --strdir "$STR_DIR" \
    --output "$OUTPUT_FILE" \
    --config "$CONFIG_FILE"