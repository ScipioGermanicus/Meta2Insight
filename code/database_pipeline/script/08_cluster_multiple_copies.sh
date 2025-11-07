#!/usr/bin/bash -l
#SBATCH --job-name=cluster_16S
#SBATCH --output=cluster_16S_%j.out
#SBATCH --error=cluster_16S_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
set -x

# env with vsearch
source /usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh
conda activate pyenv

ROOT=/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome

cluster_domain () {
  dom="$1"   # bacteria | archaea
  in_multi="$ROOT/${dom}_16S_multiple"
  in_single="$ROOT/${dom}_16S_single"
  out_dir="$ROOT/${dom}_16S_clustered"
  mkdir -p "$out_dir"

  echo "[INFO] Domain: $dom"
  echo "[INFO] in_multi:  $in_multi"
  echo "[INFO] in_single: $in_single"
  echo "[INFO] out_dir:   $out_dir"

  mapfile -t multi_files < <(find "$in_multi" -maxdepth 1 -type f \( -name '*.fna' -o -name '*.fa' \) | sort)
  mapfile -t single_files < <(find "$in_single" -maxdepth 1 -type f \( -name '*.fna' -o -name '*.fa' \) | sort)
  echo "[INFO] Found ${#multi_files[@]} multi and ${#single_files[@]} single files"

  # Cluster each multi-copy file independently at 90% id; write centroid as <MAG>.fna
  for f in "${multi_files[@]}"; do
    base=$(basename "$f")
    echo "[CLUSTER] $f -> $out_dir/$base"
    vsearch --cluster_fast "$f" --id 0.90 --centroids "$out_dir/$base" --threads "${SLURM_CPUS_PER_TASK:-4}"
  done

  # Copy single-copy genomes as-is
  if ((${#single_files[@]} > 0)); then
    cp -n "${single_files[@]}" "$out_dir"/
  fi

  echo "[INFO] ${dom}: total MAG FASTAs in clustered dir: $(ls -1 "$out_dir"/*.f*a 2>/dev/null | wc -l)"
}

cluster_domain bacteria
cluster_domain archaea

echo "[DONE] 16S clustering."
