#!/usr/bin/env bash
#SBATCH --job-name=checkm_qc
#SBATCH --output=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/checkm_%x-%j.out
#SBATCH --error=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/checkm_%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yr42on@student.aau.dk
#SBATCH --chdir=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline

set -euo pipefail
set -x

# --- Resolve project root robustly (never use /var/spool/...):
ROOT="${SLURM_SUBMIT_DIR:-/home/student.aau.dk/yr42on/Thesis/code/database_pipeline}"
# If this file lives in ROOT/script/, you can still cd there for convenience:
cd "$ROOT"

# --- Tooling & data
CHECKM_BIN="/home/student.aau.dk/yr42on/.conda/envs/checkm/bin/checkm"
ENV_BIN="/home/student.aau.dk/yr42on/.conda/envs/checkm/bin"
export PATH="$ENV_BIN:$PATH"

# CheckM data root (we verified this earlier)
export CHECKM_DATA_PATH="$HOME/checkm_data"
test -f "$CHECKM_DATA_PATH/hmms/phylo.hmm" || { echo "Bad CHECKM_DATA_PATH: $CHECKM_DATA_PATH"; exit 2; }

# --- Inputs/outputs (ALWAYS under $ROOT, not relative to spool)
SRC_BAC_DIR="$ROOT/intermediate/MAGs_formatted/genomes_to_search_barrnap/bacteria"
SRC_ARC_DIR="$ROOT/intermediate/MAGs_formatted/genomes_to_search_barrnap/archaea"
OUT_ROOT="$ROOT/intermediate/CheckM"
EXT="fna"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Create outputs
mkdir -p "${OUT_ROOT}"/{bacteria,archaea,merged,unzipped/bacteria,unzipped/archaea}

# Sanity
[ -d "$SRC_BAC_DIR" ] || { echo "Missing $SRC_BAC_DIR"; exit 2; }
[ -d "$SRC_ARC_DIR" ] || { echo "Missing $SRC_ARC_DIR"; exit 2; }

# --- Unzip if needed to standardise on .fna
unzip_if_needed () {
  local src="$1" dest="$2" ext="$3"
  shopt -s nullglob
  local any=no
  for f in "$src"/*.${ext}.gz; do any=yes; gunzip -c "$f" > "$dest/$(basename "$f" .gz)"; done
  for f in "$src"/*.${ext};    do any=yes; cp -f "$f" "$dest/"; done
  [[ $any = yes ]] || { echo "No *.$ext(.gz) in $src"; exit 2; }
}
WORK_BAC_DIR="${OUT_ROOT}/unzipped/bacteria"
WORK_ARC_DIR="${OUT_ROOT}/unzipped/archaea"
unzip_if_needed "$SRC_BAC_DIR" "$WORK_BAC_DIR" "$EXT"
unzip_if_needed "$SRC_ARC_DIR" "$WORK_ARC_DIR" "$EXT"
ls -1 "$WORK_BAC_DIR"/*."$EXT" >/dev/null
ls -1 "$WORK_ARC_DIR"/*."$EXT" >/dev/null

# --- Run CheckM
"$CHECKM_BIN" lineage_wf -x "$EXT" -t "$THREADS" "$WORK_BAC_DIR" "${OUT_ROOT}/bacteria"
"$CHECKM_BIN" lineage_wf -x "$EXT" -t "$THREADS" "$WORK_ARC_DIR" "${OUT_ROOT}/archaea"

"$CHECKM_BIN" qa "${OUT_ROOT}/bacteria/lineage.ms" "${OUT_ROOT}/bacteria" -o 2 -f "${OUT_ROOT}/bacteria/checkm_qa.tsv"
"$CHECKM_BIN" qa "${OUT_ROOT}/archaea/lineage.ms"  "${OUT_ROOT}/archaea"  -o 2 -f "${OUT_ROOT}/archaea/checkm_qa.tsv"

# --- Minimal TSVs
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {for(i=1;i<=NF;i++) h[$i]=i; print "genome_id","checkm_completeness","checkm_contamination"; next}
            {print $h["Bin Id"], $h["Completeness"], $h["Contamination"]}' \
  "${OUT_ROOT}/bacteria/checkm_qa.tsv" > "${OUT_ROOT}/bacteria/checkm_results.min.tsv"

awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {for(i=1;i<=NF;i++) h[$i]=i; print "genome_id","checkm_completeness","checkm_contamination"; next}
            {print $h["Bin Id"], $h["Completeness"], $h["Contamination"]}' \
  "${OUT_ROOT}/archaea/checkm_qa.tsv" > "${OUT_ROOT}/archaea/checkm_results.min.tsv"

# --- Merge
{
  echo -e "genome_id\tcheckm_completeness\tcheckm_contamination"
  tail -n +2 "${OUT_ROOT}/bacteria/checkm_results.min.tsv"
  tail -n +2 "${OUT_ROOT}/archaea/checkm_results.min.tsv"
} > "${OUT_ROOT}/merged/checkm_results.min.tsv"

echo "Done:"
echo "  ${OUT_ROOT}/bacteria/checkm_results.min.tsv"
echo "  ${OUT_ROOT}/archaea/checkm_results.min.tsv"
echo "  ${OUT_ROOT}/merged/checkm_results.min.tsv"
