#!/usr/bin/bash -l
#SBATCH --job-name=ssu_align_16S
#SBATCH --output=ssu_align_%j.out
#SBATCH --error=ssu_align_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -euo pipefail
set -x

# --- Activate environment ---
source /usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh
conda activate barrnap_env   # must contain Infernal + HMMER

# --- Paths ---
MODELS=~/Thesis/code/database_pipeline/intermediate/ssu_models
ROOT=~/Thesis/code/database_pipeline/intermediate/count_copies_per_genome
OUT=~/Thesis/code/database_pipeline/intermediate/ssu_align
mkdir -p "$OUT"

# --- ARCHAEA alignment ---
echo "[INFO] Aligning Archaea 16S centroids..."
cmalign --dna --matchonly --mxsize 4096 \
  "$MODELS/archaea.cm" \
  "$ROOT/archaea_16S_centroids.fasta" \
  > "$OUT/archaea_16S_centroids.stk"

# Convert to aligned FASTA
esl-reformat afa "$OUT/archaea_16S_centroids.stk" \
  > "$OUT/archaea_16S_centroids_ssu_align.fna"

# --- BACTERIA alignment ---
echo "[INFO] Aligning Bacteria 16S centroids..."
cmalign --dna --matchonly --mxsize 8192 \
  "$MODELS/bacteria.cm" \
  "$ROOT/bacteria_16S_centroids.fasta" \
  > "$OUT/bacteria_16S_centroids.stk"

# Convert to aligned FASTA
esl-reformat afa "$OUT/bacteria_16S_centroids.stk" \
  > "$OUT/bacteria_16S_centroids_ssu_align.fna"

echo "[DONE] Alignments saved in: $OUT"
