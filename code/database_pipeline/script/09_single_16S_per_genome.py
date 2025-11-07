#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os, pandas as pd

ROOT = "/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/count_copies_per_genome"
DOMAINS = ["bacteria","archaea"]

# --- OPTIONAL: enforce CheckM allow-list (set to True to enable) ---
USE_CHECKM_FILTER = False
CHECKM_TSV = "/home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/qc/checkm_filtered.tsv"

def norm_gid(fn: str) -> str:
    # normalize filenames to bare MAG id, e.g. MAG0001[...].(fa|fna) -> MAG0001
    x = fn
    for suf in ("_16S.fna","_16S.fa",".fna",".fa","_genomic.fna"):
        x = x.replace(suf,"")
    return x

def pick_longest(records):
    return sorted(records, key=lambda r: len(r.seq), reverse=True)[0] if records else None

# Load allow-list if requested
allow = None
if USE_CHECKM_FILTER:
    df = pd.read_csv(CHECKM_TSV, sep="\t")
    # your genome_id looks like MAG0078_genomic; drop trailing "_genomic" to match our MAG IDs
    df["MAG"] = df["genome_id"].str.replace(r"_genomic$", "", regex=True)
    allow = set(df["MAG"].astype(str))

for kingdom in DOMAINS:
    clustered = os.path.join(ROOT, f"{kingdom}_16S_clustered")
    singles   = os.path.join(ROOT, f"{kingdom}_16S_single")
    out_fa    = os.path.join(ROOT, f"{kingdom}_16S_genes.fasta")
    out_map   = os.path.join(ROOT, f"{kingdom}_16S_genes.map.tsv")
    os.makedirs(clustered, exist_ok=True); os.makedirs(singles, exist_ok=True)

    picked, source = {}, {}

    # 1) clustered: one file per MAG, may contain >1 centroid → pick longest
    for fn in sorted(os.listdir(clustered)):
        if not (fn.endswith(".fna") or fn.endswith(".fa")): continue
        gid = norm_gid(fn)
        if allow is not None and gid not in allow: continue
        recs = list(SeqIO.parse(os.path.join(clustered, fn), "fasta"))
        if recs:
            rec = pick_longest(recs)
            picked[gid] = SeqRecord(rec.seq, id=gid, description="")
            source[gid] = "clustered"

    # 2) singles: add if not already picked
    for fn in sorted(os.listdir(singles)):
        if not (fn.endswith(".fna") or fn.endswith(".fa")): continue
        gid = norm_gid(fn)
        if allow is not None and gid not in allow: continue
        if gid in picked: continue
        recs = list(SeqIO.parse(os.path.join(singles, fn), "fasta"))
        if recs:
            rec = pick_longest(recs)
            picked[gid] = SeqRecord(rec.seq, id=gid, description="")
            source[gid] = "single"

    # write outputs
    SeqIO.write([picked[k] for k in sorted(picked)], out_fa, "fasta")
    with open(out_map, "w") as fh:
        fh.write("genome_id\tsource\n")
        for k in sorted(picked):
            fh.write(f"{k}\t{source[k]}\n")

    print(f"{kingdom}: wrote {len(picked)} sequences → {out_fa}")
