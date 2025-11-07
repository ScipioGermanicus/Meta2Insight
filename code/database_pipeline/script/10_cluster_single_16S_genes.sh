#!/usr/bin/bash -l
#SBATCH --job-name=cluster_16S
#SBATCH --output=cluster_16S_%j.out
#SBATCH --error=cluster_16S_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
set -x

# Activate env with vsearch
source /usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh
conda activate pyenv

ROOT=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome

cluster_id1 () {
  local kingdom="$1"
  local in_fa="$ROOT/${kingdom}_16S_genes.fasta"
  local out_cent="$ROOT/${kingdom}_16S_centroids.fasta"
  local out_uc="$ROOT/${kingdom}_16S_clusters.uc"

  echo "[INFO] Clustering ${kingdom} at 100% identity"
  vsearch --cluster_fast "$in_fa" \
          --id 1.0 \
          --centroids "$out_cent" \
          --uc "$out_uc" \
          --threads "${SLURM_CPUS_PER_TASK:-4}"
}

cluster_id1 archaea
cluster_id1 bacteria

echo "[DONE] Wrote *_16S_centroids.fasta and *_16S_clusters.uc"
