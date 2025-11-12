#!/usr/bin/bash -l
#SBATCH --job-name=raxml_check_16S
#SBATCH --output=logs/raxml_check_%j.out
#SBATCH --error=logs/raxml_check_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail
mkdir -p logs raxml

# Activate env that has raxml-ng
source ~/.bashrc
conda activate barrnap_env   # or another env with raxml-ng installed

# --- Paths ---
#ARC_MSA=~/Thesis/code/database_pipeline/intermediate/ssu_align/archaea_16S_centroids_ssu_align_best.fna
BAC_MSA=~/Thesis/code/database_pipeline/intermediate/ssu_align/bacteria_16S_centroids_ssu_align_best.fna

# --- Run RAxML-NG --check (produces *.raxml.reduced.phy) ---
# Archaea
#raxml-ng --check \
#  --msa "$ARC_MSA" \
#  --model GTR+G \
#  --prefix ~/Thesis/code/database_pipeline/intermediate/raxml/archaea_raxml-check \
#  --threads ${SLURM_CPUS_PER_TASK:-8}

# Bacteria
raxml-ng --check \
  --msa "$BAC_MSA" \
  --model GTR+G \
  --prefix ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check \
  --threads ${SLURM_CPUS_PER_TASK:-8}
