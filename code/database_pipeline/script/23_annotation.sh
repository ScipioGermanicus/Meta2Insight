#!/usr/bin/bash -l
#SBATCH --job-name=emapper_mags
#SBATCH --output=logs/emapper_mags_%j.out
#SBATCH --error=logs/emapper_mags_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

ENV_NAME="eggnog_legacy"  # env with python2.7, emapper.py, prodigal

DATA_DIR=/databases/eggnog/eggnog_2020-03-31
GENOME_DIR=~/Thesis/code/database_pipeline/MAGs
PROT_DIR=~/Thesis/code/database_pipeline/intermediate/eggnog_proteins
OUTDIR=~/Thesis/code/database_pipeline/intermediate/eggnog_out
ID_MAP=~/Thesis/code/database_pipeline/intermediate/MAGs_formatted/id_map.tsv

mkdir -p "$OUTDIR" "$PROT_DIR" logs

# --- Load id_map: original_filename -> MAG0001 etc. ---
declare -A IDMAP
while IFS=$'\t' read -r orig new_id; do
    [[ "$orig" == "original_filename" ]] && continue
    [[ -z "$orig" ]] && continue
    IDMAP["$orig"]="$new_id"
done < "$ID_MAP"

echo "[INFO] Loaded ${#IDMAP[@]} mappings from $ID_MAP"

# --- Loop over MAGs ---
for f in "$GENOME_DIR"/*.fa; do
    [[ -e "$f" ]] || { echo "[ERROR] No .fa MAGs found in $GENOME_DIR"; exit 1; }

    base=$(basename "$f")          # e.g. NNF_AD15_mmlongv091.bin.1.1.fa
    base_noext="${base%.fa}"

    if [[ -n "${IDMAP[$base]:-}" ]]; then
        sample="${IDMAP[$base]}"   # e.g. MAG0001
    else
        sample="$base_noext"
        echo "[WARN] No mapping for $base in $ID_MAP, using $sample" >&2
    fi

    prot_faa="$PROT_DIR/${sample}.faa"

    echo "[INFO] Predicting proteins for $base -> $prot_faa"
    conda run -n "$ENV_NAME" prodigal \
      -i "$f" \
      -a "$prot_faa" \
      -p single

    echo "[INFO] Annotating $sample with eggnog-mapper (protein input)"
    conda run -n "$ENV_NAME" emapper.py \
      -m diamond \
      --data_dir "$DATA_DIR" \
      -i "$prot_faa" \
      -o "$sample" \
      --output_dir "$OUTDIR" \
      --cpu "${SLURM_CPUS_PER_TASK:-16}"
done

echo "[DONE] All MAGs processed with Prodigal + eggnog-mapper."
