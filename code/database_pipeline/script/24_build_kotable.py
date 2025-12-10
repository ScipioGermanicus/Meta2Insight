#!/usr/bin/env barrnap_env
"""
Build PICRUSt2-style KO trait tables (ko.txt.gz) from eggNOG-mapper annotations.

- Always requires bacterial annotations.
- Optionally takes archaeal annotations (for later when you have them).

Usage examples:

  # Bacteria only (what you can do right now)
  python make_ko_tables_optional_archaea.py \
      --bac-dir /path/to/bacteria_annotations \
      --bac-out default_files/bacteria/ko.txt.gz

  # Bacteria + archaea (later, when you have archaeal genomes)
  python make_ko_tables_optional_archaea.py \
      --bac-dir /path/to/bacteria_annotations \
      --bac-out default_files/bacteria/ko.txt.gz \
      --arc-dir /path/to/archaea_annotations \
      --arc-out default_files/archaea/ko.txt.gz

Assumptions:
  - Each genome has one eggNOG file named like:
        GENOME_ID.emapper.annotations
    (GENOME_ID is what you want PICRUSt2 to recognise as the genome ID.)
  - Files have a header line containing a "KEGG_ko" column (not necessarily
    the first line; there may be several comment lines before).
  - KEGG_ko field may be "-", empty, a single KO, or multiple KOs separated
    by commas or pipes, sometimes with "ko:" prefixes.

Output format (per domain):
  - A wide PICRUSt2-style matrix:

        assembly    ko:K00003    ko:K00008    ...
        MAG0001     1            3            ...
        MAG0002     0            2            ...

    where:
      * rows = genomes
      * columns = KOs prefixed with "ko:"
      * entries = copy numbers (gene counts)
"""

import sys
import os
import glob
import gzip
import argparse
from collections import Counter, defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build KO trait tables (ko.txt.gz) from eggNOG-mapper annotations."
    )
    parser.add_argument(
        "--bac-dir", required=True,
        help="Directory containing bacterial *.emapper.annotations files."
    )
    parser.add_argument(
        "--bac-out", required=True,
        help="Output path for bacterial ko.txt.gz."
    )
    parser.add_argument(
        "--arc-dir", required=False,
        help="Directory containing archaeal *.emapper.annotations files (optional)."
    )
    parser.add_argument(
        "--arc-out", required=False,
        help="Output path for archaeal ko.txt.gz (optional; required if --arc-dir is used)."
    )

    args = parser.parse_args()

    # Basic checks
    if not os.path.isdir(args.bac_dir):
        sys.stderr.write(f"Error: --bac-dir not found or not a directory: {args.bac_dir}\n")
        sys.exit(1)

    if args.arc_dir:
        if not os.path.isdir(args.arc_dir):
            sys.stderr.write(f"Error: --arc-dir not found or not a directory: {args.arc_dir}\n")
            sys.exit(1)
        if not args.arc_out:
            sys.stderr.write("Error: --arc-out is required when --arc-dir is provided.\n")
            sys.exit(1)

    return args


def extract_genome_id_from_filename(filename: str) -> str:
    """
    Given a filename like 'NNF_AD15_mmlongv091.bin.1.1.emapper.annotations',
    return the genome ID used downstream in PICRUSt2.

    Default behaviour:
        - take everything before '.emapper' as the ID;
        - if '.emapper' is missing, take everything before the first dot.

    Adjust this if your naming scheme is different.
    """
    base = os.path.basename(filename)
    if ".emapper" in base:
        return base.split(".emapper")[0]
    else:
        return base.split(".")[0]


def parse_eggnog_dir(ann_dir: str) -> dict:
    """
    Parse all *.emapper.annotations* files in ann_dir, and return:

        genome_id -> Counter(KO -> copy_number)

    i.e. for each genome (file) and each KO, how many times it appears.
    """
    pattern = os.path.join(ann_dir, "*.emapper.annotations*")
    ann_files = sorted(glob.glob(pattern))

    if not ann_files:
        sys.stderr.write(f"Warning: no *.emapper.annotations files found in {ann_dir}\n")

    ko_counts = defaultdict(Counter)

    for ann in ann_files:
        genome_id = extract_genome_id_from_filename(ann)
        sys.stderr.write(f"Processing {ann} -> genome_id {genome_id}\n")

        with open(ann) as f:
            ko_idx = None
            header_found = False

            for line in f:
                if line.startswith("#"):
                    # Comment or header line. We only treat it as the header
                    # if it actually contains the 'KEGG_ko' column.
                    text = line[1:].rstrip("\n")
                    cols = text.split("\t")
                    if "KEGG_ko" in cols:
                        ko_idx = cols.index("KEGG_ko")
                        header_found = True
                    # Else: just another comment line (version, command, time, etc.)
                    continue

                # From here on we are in data lines.
                if not header_found or ko_idx is None:
                    sys.stderr.write(
                        f"Error: no header line containing 'KEGG_ko' found in {ann}\n"
                    )
                    sys.exit(1)

                cols = line.rstrip("\n").split("\t")
                if len(cols) <= ko_idx:
                    # malformed line
                    continue

                kfield = cols[ko_idx].strip()
                if not kfield or kfield == "-" or kfield.upper() == "NA":
                    continue

                # Try both comma and pipe as separators for multiple KOs
                if "," in kfield:
                    raw_kos = kfield.split(",")
                elif "|" in kfield:
                    raw_kos = kfield.split("|")
                else:
                    raw_kos = [kfield]

                for ko in raw_kos:
                    ko = ko.strip()
                    if not ko or ko == "-" or ko.upper() == "NA":
                        continue
                    if ko.startswith("ko:"):
                        ko = ko[3:]

                    # Keep only canonical KEGG KO IDs like K01623
                    if not (ko.startswith("K") and len(ko) == 6 and ko[1:].isdigit()):
                        continue

                    ko_counts[genome_id][ko] += 1

    return ko_counts


def write_ko_table(ko_counts: dict, out_path: str):
    """
    Write a PICRUSt2-style wide KO table (gzipped) with format:

        assembly    ko:K00003    ko:K00008    ...
        MAG0001     1            3            ...
        MAG0002     0            2            ...

    where:
      - rows = genomes (assembly IDs),
      - columns = KOs prefixed with 'ko:',
      - entries = copy numbers (int).
    """
    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Collect the full KO universe across all genomes
    all_kos = set()
    for counter in ko_counts.values():
        all_kos.update(counter.keys())
    all_kos = sorted(all_kos)

    with gzip.open(out_path, "wt") as out:
        # Header: 'assembly' + all KOs (with 'ko:' prefix)
        header = ["assembly"] + [f"ko:{ko}" for ko in all_kos]
        out.write("\t".join(header) + "\n")

        # Rows: one per genome, with counts in the same KO order
        for genome_id in sorted(ko_counts.keys()):
            counter = ko_counts[genome_id]
            row = [genome_id]
            for ko in all_kos:
                n = counter.get(ko, 0)
                row.append(str(n))
            out.write("\t".join(row) + "\n")

    sys.stderr.write(
        f"Wrote KO matrix with {len(ko_counts)} genomes and {len(all_kos)} KOs to {out_path}\n"
    )


def main():
    args = parse_args()

    # Bacteria (always)
    sys.stderr.write(f"=== Bacteria: parsing {args.bac_dir} ===\n")
    bac_counts = parse_eggnog_dir(args.bac_dir)
    write_ko_table(bac_counts, args.bac_out)

    # Archaea (optional)
    if args.arc_dir:
        sys.stderr.write(f"=== Archaea: parsing {args.arc_dir} ===\n")
        arc_counts = parse_eggnog_dir(args.arc_dir)
        write_ko_table(arc_counts, args.arc_out)
    else:
        sys.stderr.write("No --arc-dir provided, skipping archaeal KO table.\n")


if __name__ == "__main__":
    main()
