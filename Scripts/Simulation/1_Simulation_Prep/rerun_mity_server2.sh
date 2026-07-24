#!/usr/bin/env bash
#SBATCH --job-name=rerun
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=1:00:00
#SBATCH --output=/vast/scratch/users/chen.k/%x_%j.log

set -euo pipefail
IFS=$'\n\t'

module purge

# -----------------------------
# Args
# -----------------------------
sample="${1:-}"

if [[ -z "${sample}" ]]; then
  echo 'Usage: sbatch 0_rerun_mity_server2.sh <sample>'
  echo 'Example: sbatch 0_rerun_mity_server2.sh HG01881'
  exit 1
fi

# -----------------------------
# Paths & directories
# -----------------------------
DATA_DIR="/vast/projects/bahlo_mtDNA/Results_1000G/${sample}"
# DATA_DIR="/vast/scratch/users/$USER/tmp"

OUT_ROOT=/vast/projects/bahlo_mtDNA/Simulation
MITY_OUT="${OUT_ROOT}/1000G/${sample}/mity"
MTDNA2_OUT="${OUT_ROOT}/1000G/${sample}/mtdna-server-2"

bam_file="${sample}.chrM.bam"

# -----------------------------
# Tool locations
# -----------------------------
MITY_DIR=/vast/projects/bahlo_mtDNA/software

# -----------------------------
# Pre-flight checks
# -----------------------------
[[ -f "$DATA_DIR/$bam_file" ]] || { echo "Missing BAM: $DATA_DIR/$bam_file" >&2; exit 1; }
[[ -f "$MITY_DIR/mity_latest.sif" ]] || { echo "Missing mity SIF: $MITY_DIR/mity_latest.sif" >&2; exit 1; }
[[ -f /vast/projects/bahlo_mtDNA/Simulation/mtdna-server-2/mtdna-server-2.config ]] || { echo "Missing mtDNA-Server 2 config" >&2; exit 1; }

mkdir -p "$MITY_OUT" "$MTDNA2_OUT"

# -----------------------------
# mity (call + report in one container)
# -----------------------------
module load apptainer

echo "Running mity for sample: $sample"
MITY_LOG="$MITY_OUT/${sample}.mity.log"

# Bind only what we need: chrM BAM directory and output dir
apptainer exec \
  --bind "$DATA_DIR":/data \
  --bind "$MITY_OUT":/out \
  "$MITY_DIR/mity_latest.sif" bash -c "
    set -euo pipefail
    echo 'Running mity call for $sample…'
    mity call \
      --reference hg38 \
      --prefix \"$sample\" \
      --output-dir /out \
      --min-alternate-fraction 0.005 \
      --min-alternate-count 4 \
      --normalise \"/data/$bam_file\"

    echo 'Running mity report for $sample…'
    cd /out
    mity report \
      --prefix \"$sample\" \
      --min_vaf 0 \
      --contig chrM \
      \"${sample}.normalise.vcf.gz\"

    echo 'mity completed for $sample.'
  " >"$MITY_LOG" 2>&1 
  
# -----------------------------
# mtDNA-Server 2 (Nextflow)
# -----------------------------

module load nextflow

# Use Apptainer cache for Nextflow's Singularity backend
mkdir -p /vast/scratch/users/$USER/apptainer_cache
export APPTAINER_CACHEDIR=/vast/scratch/users/$USER/apptainer_cache

echo "Starting mtDNA-Server 2 in $MTDNA2_OUT"
cd "$MTDNA2_OUT"

cp /vast/projects/bahlo_mtDNA/Simulation/mtdna-server-2/mtdna-server-2.config .
sed -i "s|^ *files *=.*|    files = \"${DATA_DIR}/$bam_file\"|" mtdna-server-2.config

# Expect a local config (mtdna-server-2.config) to provide inputs/workdir
# Add -resume to reuse completed work if re-running
nextflow run genepi/mtdna-server-2 -r v2.1.16 \
  -c mtdna-server-2.config \
  -profile singularity \
  | tee "$MTDNA2_OUT/nextflow.log"

echo "mtDNA-Server 2 launch complete."
