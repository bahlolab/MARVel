#!/usr/bin/env bash
#SBATCH --job-name=RunSim
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output= /vast/projects/bahlo_mtDNA/Simulation/scripts/logs/%x_%j.log

set -euo pipefail
IFS=$'\n\t'

module purge

# -----------------------------
# Paths & directories
# -----------------------------
sample="${1:-}"
VAF="${2:-}"

if [[ -z "${sample}" || -z "${VAF}" ]]; then
  echo 'Usage: sbatch 2_simulation.sh <sample> <VAF>'
  echo 'Example: sbatch 2_simulation.sh HG01881 10'
  exit 1
fi

DATA_DIR="/vast/projects/bahlo_mtDNA/Simulation/1000G/${sample}/AF${VAF}_DS2k"
WGS_DIR="$DATA_DIR/WGS"
CHR_DIR="$DATA_DIR/chrM"              

OUT_ROOT=/vast/projects/bahlo_mtDNA/Simulation
MITY_OUT="$OUT_ROOT/mity/${sample}_AF${VAF}"
MTDNA2_OUT="$OUT_ROOT/mtdna-server-2/${sample}_AF${VAF}"
MITOHPC_OUT="$OUT_ROOT/mitoHPC/${sample}_AF${VAF}"
MTSWIRL_OUT="$OUT_ROOT/mtSwirl/${sample}_AF${VAF}"

# -----------------------------
# bam & cram
# -----------------------------
bam_file="${sample}.chrM.10indel.40snv.af${VAF}.bam"
cram_file="${sample}.af${VAF}.cram"

# -----------------------------
# Tool locations
# -----------------------------
MITY_DIR=/vast/projects/bahlo_mtDNA/software
HP_SDIR=/stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/scripts

# -----------------------------
# Pre-flight checks
# -----------------------------
[[ -f "$CHR_DIR/$bam_file" ]] || { echo "Missing BAM: $CHR_DIR/$bam_file" >&2; exit 1; }
[[ -f "$WGS_DIR/$cram_file" ]] || { echo "Missing CRAM: $WGS_DIR/$cram_file" >&2; exit 1; }
[[ -f "$WGS_DIR/${cram_file}.crai" ]] || { echo "Missing CRAI: $WGS_DIR/${cram_file}.crai" >&2; exit 1; }

[[ -f "$MITY_DIR/mity_latest.sif" ]] || { echo "Missing mity SIF: $MITY_DIR/mity_latest.sif" >&2; exit 1; }
[[ -f /vast/projects/bahlo_mtDNA/Simulation/mtdna-server-2/mtdna-server-2.config ]] || { echo "Missing mtDNA-Server 2 config" >&2; exit 1; }

mkdir -p "$MITY_OUT" "$MTDNA2_OUT" "$MITOHPC_OUT" "$MTSWIRL_OUT"

# -----------------------------
# mity (call + report in one container)
# -----------------------------
module load apptainer

echo "Running mity for sample: $sample"
# Log file per sample
MITY_LOG="$MITY_OUT/${sample}.mity.log"

# Bind only what we need: chrM BAM directory and output dir
apptainer exec \
  --bind "$CHR_DIR":/data \
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
sed -i "s|^ *files *=.*|    files = \"${CHR_DIR}/*.bam\"|" mtdna-server-2.config

# Expect a local config (mtdna-server-2.config) to provide inputs/workdir
# Add -resume to reuse completed work if re-running
nextflow run genepi/mtdna-server-2 -r v2.1.16 \
  -c mtdna-server-2.config \
  -profile singularity \
  | tee "$MTDNA2_OUT/nextflow.log"

echo "mtDNA-Server 2 launch complete."

# -----------------------------
# MitoHPC
# -----------------------------
cd "$MITOHPC_OUT"

rm -f in.txt

# Source the init file to load environment variables
cp $OUT_ROOT/mitoHPC/init.sh .
sed -i "s|^export HP_ADIR=.*|export HP_ADIR=$WGS_DIR|" init.sh
grep HP_ADIR init.sh

. ./init.sh

# Generate commands from canonical run.sh
[[ -f "$HP_SDIR/run.sh" ]] || { echo "Cannot find MitoHPC run.sh at $HP_SDIR" >&2; exit 1; }

"$HP_SDIR/run.sh" > run.all.sh

echo "Executing MitoHPC (logging to output.log)"
bash ./run.all.sh > output.log 2>&1

echo "MitoHPC complete."

# -----------------------------
# mtSwirl
# -----------------------------
export APPTAINER_TMPDIR=/vast/scratch/users/$USER/tmp 
export APPTAINER_CACHEDIR=/vast/scratch/users/$USER/scache 
[ ! -d $APPTAINER_TMPDIR ] && mkdir $APPTAINER_TMPDIR # creates the tmp folder on vast scratch if it doesn't exist 
[ ! -d $APPTAINER_CACHEDIR ] && mkdir $APPTAINER_CACHEDIR # same as above, except for cached containers

cd $MTSWIRL_OUT
cp /vast/projects/bahlo_mtDNA/Simulation/mtSwirl/input.json .

jq --arg sample "$sample" --arg cram_file "$WGS_DIR/${cram_file}" '
  .["MitochondriaPipeline.sample_name"] = $sample |
  .["MitochondriaPipeline.wgs_aligned_input_bam_or_cram"] = $cram_file |
  .["MitochondriaPipeline.wgs_aligned_input_bam_or_cram_index"] = ($cram_file + ".crai")
' input.json > tmp.json && mv tmp.json input.json

#screen
MTSWIRL_SDIR=/vast/scratch/users/$USER/mtSwirl_HPC
source $MTSWIRL_SDIR/miniwdl_env/bin/activate
module unload apptainer
module load apptainer/1.3.5
module list # make sure apptainer loaded properly

miniwdl run $MTSWIRL_SDIR/WDL/v2.5_MongoSwirl_Single/fullMitoPipeline_v2_5_Single.wdl --input input.json | tee "$MTSWIRL_OUT/miniwdl.log"

