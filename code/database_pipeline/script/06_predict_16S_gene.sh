#!/usr/bin/bash -l
#SBATCH --job-name=barrnap
#SBATCH --output=barrnap_%j.out
#SBATCH --error=barrnap_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yr42on@student.aau.dk

set -euo pipefail
set -x  # debug: print commands

# --- Conda env ---
source /usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh
conda activate barrnap   # ensure this env has: barrnap bedtools samtools

# --- Inputs/outputs ---
BASE=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/qc/bins
OUT_BASE=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/barrnap
mkdir -p "$OUT_BASE"/{bacteria,archaea}

# Helper: stream fasta whether gzipped or not
cat_fasta() {
  case "$1" in
    *.gz)  gzip -cd "$1" ;;
    *)     cat "$1" ;;
  esac
}

run_barrnap() {
  local DOMAIN=$1   # bacteria | archaea
  local KINGDOM=$2  # bac | arc
  local INDIR="$BASE/$DOMAIN"
  local OUTDIR="$OUT_BASE/$DOMAIN"
  mkdir -p "$OUTDIR"
  shopt -s nullglob

  # sanity: ensure there are genomes to process
  files=( "$INDIR"/*_genomic.fna "$INDIR"/*_genomic.fna.gz )
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "[WARN] No genomes found in $INDIR matching *_genomic.fna(.gz)"; return 0
  fi

  echo "[INFO] Running Barrnap for $DOMAIN from $INDIR -> $OUTDIR"
  for f in "${files[@]}"; do
    [[ -e "$f" ]] || continue
    fname=$(basename "$f")
    id=${fname%_genomic.fna}; id=${id%.gz}

    # 1) Barrnap to GFF
    if [[ "$f" == *.gz ]]; then
      gzip -cd "$f" | barrnap --kingdom "$KINGDOM" --threads "${SLURM_CPUS_PER_TASK:-8}" > "$OUTDIR/${id}.gff"
    else
      barrnap --kingdom "$KINGDOM" --threads "${SLURM_CPUS_PER_TASK:-8}" < "$f" > "$OUTDIR/${id}.gff"
    fi

    # 2) Extract sequences using a real, indexed FASTA to avoid /dev/fd/*.fai issues
    tmpfa="$OUTDIR/${id}.tmp.fa"
    if [[ "$f" == *.gz ]]; then gzip -cd "$f" > "$tmpfa"; else cp "$f" "$tmpfa"; fi
    bedtools getfasta -fi "$tmpfa" -bed "$OUTDIR/${id}.gff" -fo "$OUTDIR/${id}_rRNA.fna"
    rm -f "$tmpfa"
  done
  echo "[INFO] Finished $DOMAIN"
}

# --- CALL THE FUNCTION(S)! ---
run_barrnap bacteria bac
run_barrnap archaea arc

echo "[DONE] Barrnap finished for all genomes."
