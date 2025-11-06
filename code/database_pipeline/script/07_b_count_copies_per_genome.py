#!/usr/bin/env python3
from Bio import SeqIO
import os, sys
import pandas as pd
import numpy as np

# Adjust if you changed these
ROOT = "/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome"
DOMAINS = ["bacteria", "archaea"]  # run both in one go

def is_16S(record):
    h = (record.id + " " + (record.description or "")).lower()
    # Barrnap usually tags "16S_rRNA" and description contains "16S ribosomal RNA"
    return ("16s" in h) and ("rrna" in h)

all_stats = []
for kingdom in DOMAINS:
    indir = os.path.join(ROOT, kingdom)
    out_single = os.path.join(ROOT, f"{kingdom}_16S_single")
    out_multiple = os.path.join(ROOT, f"{kingdom}_16S_multiple")
    os.makedirs(out_single, exist_ok=True)
    os.makedirs(out_multiple, exist_ok=True)

    # Collect per-genome 16S counts and write FASTAs
    copies = []
    # Expect files named like MAG0001_rRNA.fna from your Barrnap step
    fasta_files = [f for f in os.listdir(indir) if f.endswith("_rRNA.fna")]
    fasta_files.sort()
    for fn in fasta_files:
        genome_id = fn.replace("_rRNA.fna", "")               # MAG0001_genomic.fna.tmp → we used MAG0001_genomic.fna → id becomes MAG0001_genomic.fna? No: we set id earlier to MAG0001, so this yields MAG0001
        genome_id = genome_id.replace("_genomic.fna","")      # normalise: MAG0001
        path = os.path.join(indir, fn)

        seqs_16S = [rec for rec in SeqIO.parse(path, "fasta") if is_16S(rec)]
        copies.append([genome_id, len(seqs_16S)])

        if len(seqs_16S) == 1:
            SeqIO.write(seqs_16S, os.path.join(out_single, f"{genome_id}.fna"), "fasta")
        elif len(seqs_16S) > 1:
            SeqIO.write(seqs_16S, os.path.join(out_multiple, f"{genome_id}.fna"), "fasta")
        # if 0, skip (no 16S found)

    # Save copy counts and print quick stats
    copies_tsv = os.path.join(ROOT, f"{kingdom}_16S_copies.txt")
    with open(copies_tsv, "w") as fh:
        for gid, cnt in copies:
            fh.write(f"{gid}\t{cnt}\n")

    df = pd.DataFrame(copies, columns=["genome","copies"]).set_index("genome")
    df_pos = df[df["copies"] > 0]
    if not df_pos.empty:
        mean_c = float(np.mean(df_pos["copies"].values))
        med_c = float(np.median(df_pos["copies"].values))
        max_c = int(np.max(df_pos["copies"].values))
        print(f"{kingdom}: n(genomes with ≥1 16S)={df_pos.shape[0]}  mean={mean_c:.3f}  median={med_c:.1f}  max={max_c}")
    else:
        print(f"{kingdom}: no genomes with detected 16S")

print("Done.")
