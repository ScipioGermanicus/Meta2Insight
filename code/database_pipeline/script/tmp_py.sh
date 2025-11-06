#!/usr/bin/bash -l

cd ~/Thesis/code/database_pipeline/
set -euo pipefail

# Robust conda init (works even if 'conda init' wasn’t run)
if command -v conda >/dev/null 2>&1; then
  # Preferred: initialise shell hooks
  eval "$(conda shell.bash hook)" || source "$(conda info --base)/etc/profile.d/conda.sh"
else
  # Fallback: try common install paths (adjust if needed)
  source "$HOME/miniconda3/etc/profile.d/conda.sh" || \
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

conda activate pyenv
python script/05_a_quality_filtering.py
