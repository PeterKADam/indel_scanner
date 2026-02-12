#!/bin/bash
# Submit pipeline-only jobs per sample using existing RPTRF outputs.
# Usage: sbatch/submit_pipeline_only.sh
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

  if [[ ! -d "${strdir}" ]]; then
    echo "WARN: skipping ${sample}; missing STR directory: ${strdir}" >&2
    continue
  fi

  pipeline_jobid="$(sbatch <<EOF | awk '{print $4}'
#!/bin/bash
#SBATCH --account pkadmaster
#SBATCH -c 30
#SBATCH --mem 100g
#SBATCH --time 10:00:00
#SBATCH --job-name=pipe_${sample}

set -euo pipefail
if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: conda.sh not found under ${HOME}/miniconda3 or ${HOME}/anaconda3" >&2
  exit 1
fi
conda activate indel_scanner
cd "${REPO_DIR}"
python -m indel_scanner.main \
  --bam "${bam}" \
  --fasta "${fasta}" \
  --strdir "${strdir}" \
  --output "${OUTPUT_DIR}" \
  --config "${CONFIG_FILE}"
EOF
)"

  echo "sample=${sample} pipeline_jobid=${pipeline_jobid}"
done
