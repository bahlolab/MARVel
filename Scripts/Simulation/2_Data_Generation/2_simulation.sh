#!/usr/bin/env bash
#SBATCH --job-name=sim
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/vast/projects/bahlo_mtDNA/Simulation/scripts/logs/%x_%j.log

set -euo pipefail

module purge #unload all loaded environment modules
module load samtools/1.22.1
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /vast/scratch/users/$USER/conda_envs/bamsurgeon

BG="${1:-}"
VAF="${2:-}"

if [[ -z "${BG}" || -z "${VAF}" ]]; then
  echo 'Usage: sbatch 2_simulation.sh <BG> <VAF>'
  echo 'Example: sbatch 2_simulation.sh HG01881 10'
  exit 1
fi

REF=/vast/projects/bahlo_mtDNA/1000G/reference/hg38_v0_Homo_sapiens_assembly38.fasta
BG_CRAM=/vast/projects/bahlo_mtDNA/1000G/WGS/${BG}.final.cram
N1=10; N2=40

SIM_DIR="/vast/projects/bahlo_mtDNA/Simulation/1000G/${BG}"
OUTDIR=${SIM_DIR}/AF${VAF}_DS2k
TMPDIR=/vast/scratch/users/$USER/tmp
WORKDIR=$TMPDIR/bamsurgeon_${BG}_${SLURM_JOB_ID:-$$}

mkdir -p "$OUTDIR"/{chrM,WGS,vcf,logs} "$WORKDIR"

THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Python: $(which python)"; python --version
echo "addsnv:  $(which addsnv.py)"
echo "addindel:  $(which addindel.py)"

# 1) Extract chrM-only BAM from ONE background sample 
BG_CHRM_BAM="$TMPDIR/${BG}.chrM.bam"
if [[ ! -f "$BG_CHRM_BAM" ]]; then
  samtools view -@ "$THREADS" -T "$REF" -b "$BG_CRAM" chrM \
    | samtools sort -@ "$THREADS" -o "$BG_CHRM_BAM"
  samtools index -@ "$THREADS" "$BG_CHRM_BAM"
fi

# 2) Downsample to x2000 and x100
BG_CHRM_DS2k="$TMPDIR/${BG}.chrM.ds2000.bam"
BG_CHRM_DS100="$TMPDIR/${BG}.chrM.ds100.bam"

if [[ ! -f "$BG_CHRM_DS2k" || ! -f "$BG_CHRM_DS100" ]]; then

  MeanDepth=$(samtools depth -a "$BG_CHRM_BAM" \
    | awk '{sum+=$3; n++} END {printf "%.2f", sum/n}')
  echo "MeanDepth=${MeanDepth}"
  
  DS2k=$(awk -v md="$MeanDepth" 'BEGIN {
    v = 2000/md; if (v > 1) v = 1;
    printf "%.2f", v
  }')
  echo "DS2k=${DS2k}"

  DS100=$(awk -v md="$MeanDepth" 'BEGIN {
    v = 100/md; if (v > 1) v = 1;
    printf "%.2f", v
  }')
  echo "DS100=${DS100}"

  seed=42
  DS2k_frac=${DS2k#0.}
  DS100_frac=${DS100#0.}

  echo "Downsampling to ~2000×"
   samtools view -@ "$THREADS" -s "${seed}.${DS2k_frac}" -b "$BG_CHRM_BAM" \
    | samtools sort -@ "$THREADS" -o "$BG_CHRM_DS2k"
  samtools index "$BG_CHRM_DS2k"
  
  echo "Downsampling to ~100×"
  samtools view -@ "$THREADS" -s "${seed}.${DS100_frac}" -b "$BG_CHRM_BAM" \
    | samtools sort -@ "$THREADS" -o "$BG_CHRM_DS100"
  samtools index "$BG_CHRM_DS100"
  
  echo -n "Mean depth DS2k: "
  samtools depth -a "$BG_CHRM_DS2k" \
    | awk '{sum+=$3; n++} END {printf "%.2f\n", sum/n}'

  echo -n "Mean depth DS100: "
  samtools depth -a "$BG_CHRM_DS100" \
    | awk '{sum+=$3; n++} END {printf "%.2f\n", sum/n}'

fi

# 3) Targeted variants
indel="${SIM_DIR}/target.10indel.txt"
indelFILE="$OUTDIR/target.indel.${VAF}.txt"

SNV="${SIM_DIR}/target.40snv.txt"
SNVFILE="$OUTDIR/target.snv.${VAF}.txt"

awk -v vaf="$VAF" 'BEGIN{OFS="\t"}
{
  $4 = sprintf("%.2f", (vaf+0)/100)
  if ($NF=="NA") NF--
  print
}' "$indel" > "$indelFILE"

awk -v vaf="$VAF" 'BEGIN{OFS="\t"} {$4=sprintf("%.2f",(vaf+0)/100); print}' "$SNV" > "$SNVFILE"

# Quick sanity checks
echo "indelFILE lines: $(wc -l < "$indelFILE")"
head -n 10 "$indelFILE"

echo "SNVFILE lines: $(wc -l < "$SNVFILE")"
head -n 6 "$SNVFILE"

# 4) Add indels + SNVs
### Indels
INDEL_RAW="$WORKDIR/${BG}.chrM.${N1}indel.af${VAF}.raw.bam"
INDEL_SORT="$WORKDIR/${BG}.chrM.${N1}indel.af${VAF}.sort.bam"
LOG="$OUTDIR/logs/${BG}_${N1}indel_af${VAF}.log"

cd "$WORKDIR"

addindel.py \
-v "$indelFILE" \
-f "$BG_CHRM_DS2k" \
-r "$REF" \
-o "$INDEL_RAW" \
--maxdepth 5000 \
--seed 1234 \
> "$LOG" 2>&1

grep -nE 'WARNING|dropped for outcover/incover' "$LOG" || true

# sanity check: ensure addindel produced output
[[ -s "$INDEL_RAW" ]] || { echo "ERROR: addindel failed, no BAM produced"; exit 1; }

samtools sort -@ "$THREADS" -o "$INDEL_SORT" "$INDEL_RAW"
samtools index -@ "$THREADS" "$INDEL_SORT"

### SNVs
SPIKED_RAW="$WORKDIR/${BG}.chrM.${N2}snv.af${VAF}.raw.bam"
LOG="$OUTDIR/logs/${BG}_${N2}snv_af${VAF}.log"

addsnv.py \
-v "$SNVFILE" \
-f "$INDEL_SORT" \
-r "$REF" \
-o "$SPIKED_RAW" \
--maxdepth 5000 \
--seed 1234 \
> "$LOG" 2>&1

grep -nE 'WARNING|dropped for outcover/incover' "$LOG" || true

cp "$WORKDIR"/*.addindel*.vcf "$OUTDIR/vcf/" 2>/dev/null || true
cp "$WORKDIR"/*.addsnv*.vcf  "$OUTDIR/vcf/" 2>/dev/null || true

# sort/index
SPIKED_BAM="$OUTDIR/chrM/${BG}.chrM.${N1}indel.${N2}snv.af${VAF}.bam"
samtools sort -@ "$THREADS" -o "$SPIKED_BAM" "$SPIKED_RAW"
samtools index -@ "$THREADS" "$SPIKED_BAM"

echo "Spiked BAM: $SPIKED_BAM"
echo "Log:        $LOG"

samtools depth -a -r chrM "$BG_CHRM_DS2k" > "$OUTDIR/chrM.incover.txt"
samtools depth -a -r chrM "$SPIKED_BAM" > "$OUTDIR/chrM.outcover.txt"

# 5) Merge nDNA + mtDNA
# --- Convert BAM to CRAM ---
echo "[$(date)] Converting BAM to CRAM..."
samtools view -C -T "$REF" \
  --output-fmt-option version=3.0 \
  -o "$WORKDIR/${BG}.chrM.af${VAF}.cram" \
  "$SPIKED_BAM"

samtools index "$WORKDIR/${BG}.chrM.af${VAF}.cram"

# --- Make A nuclear-only CRAM (drop chrM, 20 mins) ---
if [[ ! -f "$SIM_DIR/${BG}_nochrM.cram" ]]; then
    echo "[$(date)] Creating nuclear-genome-only CRAM for A..."
  samtools view -@ "$THREADS" -T "$REF" -C \
  --output-fmt-option version=3.0 \
  -e 'rname!="chrM"' \
  -o "$SIM_DIR/${BG}_nochrM.cram" \
  "$BG_CRAM"
fi

# Merge nuclear genome + new mixed chrM (20 mins)
echo "[$(date)] Combining nuclear genome with simulated chrM..."

samtools merge -f -@ "$THREADS" -O CRAM \
  "$WORKDIR/${BG}.af${VAF}.cram" \
  "$SIM_DIR/${BG}_nochrM.cram" \
  "$WORKDIR/${BG}.chrM.af${VAF}.cram"

# Force final WGS CRAM to version 3.0 for old GATK/Picard compatibility
samtools view -@ "$THREADS" -T "$REF" -C \
  --output-fmt-option version=3.0 \
  -o "$OUTDIR/WGS/${BG}.af${VAF}.cram" \
  "$WORKDIR/${BG}.af${VAF}.cram"

samtools index -@ "$THREADS" "$OUTDIR/WGS/${BG}.af${VAF}.cram"
hexdump -C "$OUTDIR/WGS/${BG}.af${VAF}.cram" | head -n 1
