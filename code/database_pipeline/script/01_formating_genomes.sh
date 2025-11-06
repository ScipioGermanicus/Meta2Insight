#!/usr/bin/env bash

cd ~/Thesis/code/database_pipeline/MAGs

set -euo pipefail

# Config
SRC_DIR=~/Thesis/code/database_pipeline/MAGs
DST_DIR=~/Thesis/code/database_pipeline/intermediate/MAGs_formatted
OUT_PREFIX="MAG"          # will produce MAG0001, MAG0002, ...
PAD=4                     # zero padding width

mkdir -p "$DST_DIR"
cd "$DST_DIR"

# Mapping file: original filename -> new_id
MAP_FILE=id_map.tsv
: > "$MAP_FILE"           # truncate on each run (remove ":" to append instead)
printf "original_filename\tnew_id\n" >> "$MAP_FILE"

# Collect input files deterministically and robustly
# Accept .fa, .fna, and their .gz variants
mapfile -t FILES < <(find "$SRC_DIR" -maxdepth 1 -type f \
  \( -iname '*.fa' -o -iname '*.fna' -o -iname '*.fa.gz' -o -iname '*.fna.gz' \) \
  -printf '%f\n' | LC_ALL=C sort)

if (( ${#FILES[@]} == 0 )); then
  echo "No MAG FASTA files found in $SRC_DIR" >&2
  exit 1
fi

i=1
for base in "${FILES[@]}"; do
  new_id=$(printf "%s%0${PAD}d" "$OUT_PREFIX" "$i")
  tgt="${new_id}_genomic.fna.gz"

  src="$SRC_DIR/$base"

  # Always produce gzipped output with the exact suffix
  case "$base" in
    *.gz)
      # Input already gzipped: copy and rename to target name
      # Use -f to overwrite on re-runs, remove -f if you prefer to fail instead.
      cp -f "$src" "$tgt"
      ;;
    *)
      # Input not gzipped: copy then gzip
      cp -f "$src" "${new_id}_genomic.fna"
      gzip -f "${new_id}_genomic.fna"
      ;;
  esac

  # Record mapping (original basename, not full path)
  printf "%s\t%s\n" "$base" "$new_id" >> "$MAP_FILE"

  ((i++))
done

echo "Renamed $((${i}-1)) genomes into $(pwd) with pattern ${OUT_PREFIX}%0${PAD}d_genomic.fna.gz"
echo "Mapping written to: $MAP_FILE"

