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

"/home/peter/pkadmaster/indel_scanner/sbatch/run_rptrf_sbatch.sh" "$@"
