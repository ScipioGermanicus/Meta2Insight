#!/usr/bin/bash -l
#SBATCH --job-name=copy_picrust_ref
#SBATCH --output=logs/copy_picrust_%j.out
#SBATCH --error=logs/copy_picrust_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
set -euo pipefail

BASE=~/Thesis/code/database_pipeline/intermediate
RAXML=$BASE/raxml
IQTREE=$BASE/iqtree

# Output directory for PICRUSt2 reference DB
OUT=~/Thesis/code/database_pipeline/intermediate/gtdb_r220_picrust_ref
BAC_REF=$OUT/bac_ref
ARC_REF=$OUT/arc_ref

mkdir -p "$BAC_REF" "$ARC_REF" logs

echo "[INFO] Copying bacteria reference files..."

# --- Bacteria files ---
# 1. 16S alignment used for HMM
cp "$RAXML/bacteria_raxml-check.raxml.reduced_dna.fna" \
   "$BAC_REF/bac_ref.fna"

# 2. HMM built from that alignment
cp "$RAXML/bacteria_raxml-check.raxml.reduced_dna.hmm" \
   "$BAC_REF/bac_ref.hmm"

# 3. Tree (your IQ-TREE tree or bestTree)
cp "$IQTREE/bacteria_16S_iqtree.treefile" \
   "$BAC_REF/bac_ref.tre"

# 4. Best model file from RAxML-NG
cp "$RAXML/bacteria_raxml.raxml.bestModel" \
   "$BAC_REF/bac_ref.model"

# 5. raxml_info (hand-constructed)
cp "$RAXML/bacteria_raxml.raxml_info" \
   "$BAC_REF/bac_ref.raxml_info"

echo "[OK] Bacteria reference copied to $BAC_REF"

# --- Archaea block (ENABLE LATER) ---
# echo "[INFO] Copying archaeal reference files..."
# cp "$RAXML/archaea_raxml-check.raxml.reduced_dna.fna"  "$ARC_REF/arc_ref.fna"
# cp "$RAXML/archaea_raxml-check.raxml.reduced_dna.hmm"  "$ARC_REF/arc_ref.hmm"
# cp "$IQTREE/archaea_16S_iqtree.treefile"                "$ARC_REF/arc_ref.tre"
# cp "$RAXML/archaea_raxml.raxml.bestModel"              "$ARC_REF/arc_ref.model"
# cp "$RAXML/archaea_raxml.raxml_info"                   "$ARC_REF/arc_ref.raxml_info"
# echo "[OK] Archaea reference copied to $ARC_REF"

echo "[DONE] PICRUSt2 reference folder structure created."
