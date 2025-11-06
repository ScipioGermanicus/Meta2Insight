#!/bin/bash
#SBATCH --job-name=pipeline02_gtdbtk
#SBATCH --output=pipeline02_gtdbtk_%j.out
#SBATCH --error=pipeline02_gtdbtk_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00

cd ~/Thesis/code/database_pipeline/intermediate/GTDB-Tk
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate gtdbtk

export GTDBTK_DATA_PATH=/databases/GTDB/gtdbtk_packages/GTDB-TK_release226_2025_04_11

export TMPDIR=/scratch/$USER/$SLURM_JOB_ID
mkdir -p "$TMPDIR"

GENOMES=~/Thesis/code/database_pipeline/intermediate/MAGs_formatted
OUT=~/Thesis/code/database_pipeline/GTDB-Tk/

gtdbtk classify_wf \
  --skip_ani_screen \
  --genome_dir "$GENOMES" \
  --extension fa \
  --out_dir "$OUT" \
  --cpus 8 \
  --pplacer_cpus 2

