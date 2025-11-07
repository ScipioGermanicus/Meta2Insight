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
set -x

# conda env with bedtools
source /usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh
conda activate barrnap_env

# --- paths (adjust if yours differ) ---
BASE=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/qc/bins
GFF_DIR=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/barrnap
OUT=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome
mkdir -p "$OUT"/{bacteria,archaea}

find_genome() {
  # $1 = domain (bacteria/archaea), $2 = bare MAG id (e.g. MAG0052)
  local dom="$1"; local id="$2"
  for cand in \
    "$BASE/$dom/${id}_genomic.fna.gz" \
    "$BASE/$dom/${id}_genomic.fna" \
    "$BASE/$dom/${id}.fa.gz" \
    "$BASE/$dom/${id}.fa" \
    "$BASE/$dom/${id}.fna.gz" \
    "$BASE/$dom/${id}.fna"
  do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
}

shopt -s nullglob
for dom in bacteria archaea; do
  echo "[INFO] Domain: $dom"
  gffs=( "$GFF_DIR/$dom"/*.gff )
  echo "[INFO] Found ${#gffs[@]} GFF files under $GFF_DIR/$dom"

  for gff in "${gffs[@]}"; do
    raw_id=$(basename "$gff" .gff)     # could be MAG0052 or MAG0052_genomic.fna
    id=${raw_id%_genomic.fna}          # normalise to MAG0052

    genome=$(find_genome "$dom" "$id" || true)
    if [[ -z "${genome:-}" ]]; then
      echo "[WARN] Missing genome for $id in $BASE/$dom"; continue
    fi

    bed="$OUT/$dom/${id}.16S.bed"
    outfa="$OUT/$dom/${id}_16S.fna"
    tmpfa="$OUT/$dom/${id}.tmp.fa"

    # 16S selection: accept “16S” OR “SSU” anywhere on rRNA features (case-insensitive)
    awk 'BEGIN{IGNORECASE=1; OFS="\t"} $0!~/^#/ && ($3=="rRNA" || $3=="RNA") && (/16S/ || /SSU/) {print $1, $4-1, $5, "16S_rRNA","0",$7}' "$gff" > "$bed"
    if [[ ! -s "$bed" ]]; then
      echo "[INFO] No 16S/SSU rows in $(basename "$gff"); skipping."
      rm -f "$bed"
      continue
    fi

    # real temp FASTA to avoid /dev/fd index warnings
    case "$genome" in
      *.gz) gzip -cd "$genome" > "$tmpfa" ;;
      *)     cp "$genome"     "$tmpfa" ;;
    esac
    rm -f "$tmpfa.fai"

    # extract sequences; -s keeps strand, -name uses BED name in header
    bedtools getfasta -fi "$tmpfa" -bed "$bed" -s -name -fo "$outfa"
    nseq=$(grep -c '^>' "$outfa" || true)
    echo "[OK] $id → $nseq 16S sequences → $(basename "$outfa")"

    rm -f "$tmpfa" "$bed"
  done
done

echo "[DONE] Wrote *_16S.fna where 16S/SSU was found."

# --- Run the counting script (Python) ---

conda activate pyenv
python3 /home/student.aau.dk/yr42on/Thesis/code/database_pipeline/script/07_b_count_copies_per_genome.py
