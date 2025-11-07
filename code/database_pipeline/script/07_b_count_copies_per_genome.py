#!/usr/bin/env python3
from Bio import SeqIO
import os, sys
import pandas as pd
import numpy as np

ROOT = "/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome"
DOMAINS = ["bacteria","archaea"]

def list_files(d):
    try:
        return sorted([f for f in os.listdir(d) if f.endswith("_16S.fna") or f.endswith("_16S.fa")])
    except FileNotFoundError:
        return []

for kingdom in DOMAINS:
    indir = os.path.join(ROOT, kingdom)
    out_single = os.path.join(ROOT, f"{kingdom}_16S_single")
    out_multiple = os.path.join(ROOT, f"{kingdom}_16S_multiple")
    os.makedirs(out_single, exist_ok=True)
    os.makedirs(out_multiple, exist_ok=True)

    files = list_files(indir)
    print(f"[INFO] {kingdom}: searching in {indir}")
    print(f"[INFO] {kingdom}: found {len(files)} *_16S.(fa|fna) files")
    if not files:
        # helpful debug: show what's actually in the folder
        try:
            print(f"[DEBUG] {kingdom} contents:", os.listdir(indir))
        except Exception as e:
            print(f"[DEBUG] cannot list {indir}: {e}")
        continue

    copies = []
    for fn in files:
        gid = fn.replace("_16S.fna","").replace("_16S.fa","").replace("_genomic.fna","")
        path = os.path.join(indir, fn)
        seqs = list(SeqIO.parse(path, "fasta"))
        n = len(seqs)
        copies.append([gid, n])
        if n == 1:
            SeqIO.write(seqs, os.path.join(out_single, f"{gid}.fna"), "fasta")
        elif n > 1:
            SeqIO.write(seqs, os.path.join(out_multiple, f"{gid}.fna"), "fasta")

    # write summary and print stats
    tsv = os.path.join(ROOT, f"{kingdom}_16S_copies.txt")
    with open(tsv, "w") as fh:
        for gid, n in copies:
            fh.write(f"{gid}\t{n}\n")

    df = pd.DataFrame(copies, columns=["genome","copies"]).set_index("genome")
    pos = df[df["copies"] > 0]
    if len(pos):
        print(f"{kingdom}: n={len(pos)} mean={pos['copies'].mean():.3f} median={pos['copies'].median():.1f} max={int(pos['copies'].max())}")
    else:
        print(f"{kingdom}: no genomes with 16S")
print("Done.")
