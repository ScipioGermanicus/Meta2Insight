#!/usr/bin/bash -l
#SBATCH --job-name=raxml_eval_bac
#SBATCH --output=logs/raxml_eval_bac_%j.out
#SBATCH --error=logs/raxml_eval_bac_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --hint=nomultithread

set -euo pipefail
set -x

# Cap all implicit threading to 1 to avoid oversubscription
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export GOTO_NUM_THREADS=1
export RAYON_NUM_THREADS=1

BAC_PHYLIP=~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced.phy
BAC_TREE=~/Thesis/code/database_pipeline/intermediate/iqtree/bacteria_16S_iqtree.treefile
BAC_PREFIX=~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml

ls -l "$BAC_PHYLIP" "$BAC_TREE"

# Run in your conda env without relying on activation state
conda run -n barrnap_env raxml-ng --version
conda run -n barrnap_env raxml-ng --evaluate \
  --msa "$BAC_PHYLIP" \
  --tree "$BAC_TREE" \
  --model GTR+G \
  --prefix "$BAC_PREFIX" \
  --threads "${SLURM_CPUS_PER_TASK:-8}"
