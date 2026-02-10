#!/bin/bash
set -euo pipefail

# Usage:
#   ./extract_igv_read_bundle.sh [OUTPUT_DIR]
#
# Optional environment overrides:
#   BAM=/path/to/file.bam
#   FASTA=/path/to/file.fa
#   READ="read_name"
#   REGION="contig:start-end"
#
# Produces:
#   - read_of_interest.sam
#   - read_by_name.sam
#   - locus_ref.fa
#   - read_summary.txt

OUTPUT_DIR="${1:-./igv_read_bundle}"
mkdir -p "$OUTPUT_DIR"

BAM="${BAM:-/home/peterkad/pkadmaster/data/mutationalscanning_bam/ph/diploid_assembly/ph_plus_unmapped_diploid_v2.bam}"
FASTA="${FASTA:-/home/peterkad/pkadmaster/data/ph/ph_diploid.fa}"
READ="${READ:-m84108_250704_184541_s4/179770115/ccs}"
REGION="${REGION:-h1tg000033l:835294-835499}"

echo "BAM:    $BAM"
echo "FASTA:  $FASTA"
echo "READ:   $READ"
echo "REGION: $REGION"
echo "OUT:    $OUTPUT_DIR"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi
if [[ ! -f "$FASTA" ]]; then
  echo "ERROR: FASTA not found: $FASTA" >&2
  exit 1
fi

samtools view -h "$BAM" "$REGION" \
  | awk -v q="$READ" '(/^@/) || ($1==q)' \
  > "$OUTPUT_DIR/read_of_interest.sam"

samtools view -h "$BAM" \
  | awk -v q="$READ" '(/^@/) || ($1==q)' \
  > "$OUTPUT_DIR/read_by_name.sam"

samtools faidx "$FASTA" "$REGION" > "$OUTPUT_DIR/locus_ref.fa"

samtools view "$BAM" "$REGION" \
  | awk -v q="$READ" '$1==q{print "QNAME="$1"\tFLAG="$2"\tRNAME="$3"\tPOS="$4"\tMAPQ="$5"\tCIGAR="$6}' \
  > "$OUTPUT_DIR/read_summary.txt"

echo "Done. Files written to: $OUTPUT_DIR"
ls -1 "$OUTPUT_DIR"
