#!/usr/bin/env barrnap_env
import sys
import gzip

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <ko_txt_gz_in> <id_map.tsv> <ko_txt_gz_out>")
    sys.exit(1)

ko_in = sys.argv[1]       # ko.txt.gz (old IDs)
id_map_file = sys.argv[2] # original_filename  new_id
ko_out = sys.argv[3]      # corrected ko.txt.gz

# Read ID map
id_map = {}
with open(id_map_file) as f:
    next(f)  # skip header
    for line in f:
        orig, new = line.strip().split("\t")
        # remove .fa extension so it matches genome_id in KO table
        orig_no_ext = orig.replace(".fa", "")
        id_map[orig_no_ext] = new

# Rewrite KO table
with gzip.open(ko_in, "rt") as fin, gzip.open(ko_out, "wt") as fout:
    header = next(fin).rstrip("\n")
    fout.write(header + "\n")

    for line in fin:
        genome_id, ko, count = line.rstrip("\n").split("\t")

        if genome_id not in id_map:
            raise ValueError(
                f"Genome ID '{genome_id}' not found in ID map. "
                "Check mapping file or naming consistency."
            )

        new_id = id_map[genome_id]
        fout.write(f"{new_id}\t{ko}\t{count}\n")

print(f"Corrected KO table written to: {ko_out}")
