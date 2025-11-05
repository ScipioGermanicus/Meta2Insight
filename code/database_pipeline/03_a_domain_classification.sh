#!/usr/bin/env bash

# Domain classification



# specify working directory


# Example after producing domain map "domain_map.tsv" with two columns: genome_id \t domain
# and files named genome_id_genomic.fna.gz
while IFS=$'\t' read -r gid dom; do
  if [[ "$dom" == "Bacteria" ]]; then
    mv "${gid}_genomic.fna.gz" genomes_to_search_barrnap/bacteria/ || true
  else
    mv "${gid}_genomic.fna.gz" genomes_to_search_barrnap/archaea/ || true
  fi
done < domain_map.tsv
