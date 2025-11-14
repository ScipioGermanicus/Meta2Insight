#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import os

# =======================
# CONFIG (hard-coded paths)
# =======================

BASE = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate')

# Input
BAC_COPIES_IN = os.path.join(BASE, 'count_copies_per_genome/bacteria_16S_copies.txt')
BAC_FASTA     = os.path.join(BASE, 'gtdb_r220_picrust_ref/bac_ref/bac_ref.fna')

# Output
BAC_COPIES_OUT = os.path.join(BASE, 'gtdb_r220_picrust_ref/bac_ref/bacteria_16S_copies.txt')

# For later:
# ARC_COPIES_IN  = os.path.join(BASE, 'count_copies_per_genome/archaea_16S_copies.txt')
# ARC_FASTA      = os.path.join(BASE, 'gtdb_r220_picrust_ref/arc_ref/arc_ref.fna')
# ARC_COPIES_OUT = os.path.join(BASE, 'gtdb_r220_picrust_ref/arc_ref/archaea_16S_copies.txt')
# =======================

def filter_copy_table(fasta, copy_file, out_file):
    print(f"[INFO] Filtering: {copy_file}")
    print(f"[INFO] Using genomes from: {fasta}")

    # 1. Collect genome IDs from the final reference FASTA
    included = [rec.id for rec in SeqIO.parse(fasta, 'fasta')]
    included_set = set(included)
    print(f"[INFO] Genomes in reference: {len(included_set)}")

    # 2. Load table of copy numbers
    df = pd.read_csv(copy_file, index_col=0, header=None, sep="\t")
    print(f"[INFO] Input copy table entries: {df.shape[0]}")

    # 3. Filter to only genomes in bac_ref.fna
    df = df.loc[df.index.intersection(included_set), :]

    # 4. Reset index and name columns
    df = df.reset_index()
    df.columns = ['assembly', '16S_rRNA_Count']

    # 5. Cap counts above 10
    df.loc[df['16S_rRNA_Count'] > 10, '16S_rRNA_Count'] = 10

    # 6. Write output
    df.to_csv(out_file, index=False, sep="\t")
    print(f"[OK] Wrote filtered table: {out_file}")
    print(f"[INFO] Final rows: {df.shape[0]}")

def main():
    # Bacteria
    filter_copy_table(BAC_FASTA, BAC_COPIES_IN, BAC_COPIES_OUT)

    # Archaea (enable later)
    # filter_copy_table(ARC_FASTA, ARC_COPIES_IN, ARC_COPIES_OUT)

if __name__ == "__main__":
    main()
