#!/usr/bin/bash -l
#SBATCH --job-name=extract_16S
#SBATCH --output=extract_16S_%j.out
#SBATCH --error=extract_16S_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yr42on@student.aau.dk

set -euo pipefail

# --- Activate conda environment ---
source /usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh
conda activate barrnap_env

# --- Run the 16S extraction step ---
BASE=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/qc/bins
GFF_DIR=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/barrnap
OUT=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome
mkdir -p "$OUT"/{bacteria,archaea}

shopt -s nullglob
for dom in bacteria archaea; do
  echo "[INFO] Domain: $dom"
  for gff in "$GFF_DIR/$dom"/*.gff; do
    raw_id=$(basename "$gff" .gff)
    id=${raw_id%_genomic.fna}

    genome_gz="$BASE/$dom/${id}_genomic.fna.gz"
    genome_plain="$BASE/$dom/${id}_genomic.fna"
    bed="$OUT/$dom/${id}.16S.bed"
    outfa="$OUT/$dom/${id}_16S.fna"

    awk 'BEGIN{OFS="\t"} $0!~/^#/ && $3 ~ /rRNA/ && $0 ~ /16S/ {print $1, $4-1, $5, "16S_rRNA","0",$7}' "$gff" > "$bed"
    if [[ ! -s "$bed" ]]; then rm -f "$bed"; continue; fi

    tmpfa="$OUT/$dom/${id}.tmp.fa"
    if [[ -f "$genome_gz" ]]; then
      gzip -cd "$genome_gz" > "$tmpfa"
    elif [[ -f "$genome_plain" ]]; then
      cp "$genome_plain" "$tmpfa"
    else
      echo "[WARN] Missing genome for $id in $BASE/$dom"; rm -f "$bed"; continue
    fi

    bedtools getfasta -fi "$tmpfa" -bed "$bed" -s -name -fo "$outfa"
    rm -f "$tmpfa" "$bed"
  done
done
echo "[DONE] Wrote *_16S.fna where 16S was found."

# --- Run the counting script (Python) ---
python /home/student.aau.dk/yr42on/Thesis/code/database_pipeline/script/07_b_count_copies_per_genome.py
