#!/bin/bash
# Submit RPTRF jobs and dependent pipeline jobs per sample.
# Usage: sbatch/submit_rptrf_and_pipeline.sh
set -euo pipefail

SAMPLES=(chk da1 da4 la la4 ph tr tased)

DATA_DIR="/home/peterkad/pkadmaster/data"
BAM_ROOT="/home/peterkad/pkadmaster/data/mutationalscanning_bam"
OUTPUT_DIR="/home/peterkad/pkadmaster/indel_scanner/results"
CONFIG_FILE="config.yaml"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

fasta_path() {
  local sample="$1"
  echo "${DATA_DIR}/${sample}/${sample}_diploid.fa"
}

bam_path() {
  local sample="$1"
  case "$sample" in
    da4)
      echo "${BAM_ROOT}/da4/diploid_assembly/da4_diploid_plus_revio.bam"
      ;;
    *)
      echo "${BAM_ROOT}/${sample}/diploid_assembly/${sample}_plus_unmapped_diploid_v2.bam"
      ;;
  esac
}

for sample in "${SAMPLES[@]}"; do
  fasta="$(fasta_path "$sample")"
  bam="$(bam_path "$sample")"
  strdir="${DATA_DIR}/${sample}/repeatregions"

  rptrf_jobid="$(sbatch <<EOF | awk '{print \$4}'
#!/bin/bash
#SBATCH --account pkadmaster
#SBATCH -c 10
#SBATCH --mem 300g
#SBATCH --job-name=rptrf_${sample}

set -euo pipefail
cd "${REPO_DIR}"
./run_rptrf.sh --fasta "${fasta}" --sample "${sample}" --data-dir "${DATA_DIR}"
EOF
)"

  pipeline_jobid="$(sbatch --dependency=afterok:${rptrf_jobid} <<EOF | awk '{print \$4}'
#!/bin/bash
#SBATCH --account pkadmaster
#SBATCH -c 30
#SBATCH --mem 100g
#SBATCH --job-name=pipe_${sample}

set -euo pipefail
cd "${REPO_DIR}"
python -m indel_scanner.main \
  --bam "${bam}" \
  --fasta "${fasta}" \
  --strdir "${strdir}" \
  --output "${OUTPUT_DIR}" \
  --config "${CONFIG_FILE}"
EOF
)"

  echo "sample=${sample} rptrf_jobid=${rptrf_jobid} pipeline_jobid=${pipeline_jobid}"
done
