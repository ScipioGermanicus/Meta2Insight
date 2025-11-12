#!/usr/bin/bash -l
#SBATCH --job-name=iqtree_bac
#SBATCH --output=logs/iqtree_bac_%j.out
#SBATCH --error=logs/iqtree_bac_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G

set -euo pipefail

MSA=~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced.phy
OUTDIR=~/Thesis/code/database_pipeline/intermediate/iqtree
PREFIX="$OUTDIR/bacteria_16S_iqtree"

mkdir -p "$OUTDIR" logs

# Optional sanity checks
ls -l "$MSA"
conda run -n barrnap_env iqtree --version

# Build the tree
conda run -n barrnap_env iqtree \
  -s "$MSA" \
  -m GTR+G \
  -bb 1000 \
  -nt "${SLURM_CPUS_PER_TASK:-8}" \
  -seed 12345 \
  -redo \
  -pre "$PREFIX"
