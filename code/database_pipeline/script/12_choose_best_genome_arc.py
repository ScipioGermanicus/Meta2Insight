#!/usr/bin/env python3
"""
Choose the best genome per 16S cluster (ARCHAEA) with hard-coded paths only.
This simplified version avoids any duplicate path definitions and is ready to run directly.
Edit the CONFIG section with your correct paths, then execute:
    python 12_choose_best_genome.py
"""

import os
import sys
import csv
from collections import defaultdict
import pandas as pd
from Bio import SeqIO

# ==========================
# CONFIG — EDIT THESE PATHS
# ==========================
DOMAIN       = 'archaea'  # 'archaea' or 'bacteria' (used for metadata filtering + filenames)
CLUSTERS     = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/count_copies_per_genome/archaea_16S_clusters.uc')
ALIGNED_FASTA= os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/ssu_align/archaea_16S_centroids_ssu_align.fna')
METADATA     = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/qc/checkm_filtered.tsv')
ID_MAP       = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/MAGs_formatted/id_map.tsv')
OUTDIR       = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/choose_best_genome/')

# ==========================

def read_id_map(path: str) -> dict:
    mp = {}
    with open(path, newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        cols = {c.lower(): c for c in r.fieldnames}
        of = cols.get('original_filename') or cols.get('original') or r.fieldnames[0]
        ni = cols.get('new_id') or cols.get('id') or r.fieldnames[1]
        for row in r:
            mp[row[of]] = row[ni]
    return mp

def parse_clusters(path: str, genes_in_alignment: set) -> dict:
    clusters = defaultdict(list)
    with open(path, 'r') as f:
        for raw in f:
            if not raw.strip() or raw.startswith('#'):
                continue
            line = raw.rstrip('\n')
            parts = line.split('\t')
            rec_type = parts[0] if parts else ''

            # vsearch .uc cases
            if rec_type in {'S', 'C'}:
                q = parts[8] if len(parts) > 8 else None
                if q and q in genes_in_alignment:
                    clusters.setdefault(q, [])
                continue
            if rec_type == 'H':
                q = parts[8] if len(parts) > 8 else None
                t = parts[9] if len(parts) > 9 else None
                if t in genes_in_alignment:
                    clusters.setdefault(t, [])
                    if q is not None:
                        clusters[t].append(q)
                continue

            # fallback: last two columns [centroid, member]
            if len(parts) >= 2 and parts[-1] != '*':
                centroid, member = parts[-2], parts[-1]
                if centroid in genes_in_alignment:
                    clusters.setdefault(centroid, []).append(member)

    for gid in genes_in_alignment:
        clusters.setdefault(gid, [])
    return clusters

def choose_best(md: pd.DataFrame, candidates: list) -> str:
    avail = [g for g in candidates if g in md.index]
    if not avail:
        return sorted(candidates)[0]
    sub = md.loc[avail, :].copy()
    comp = pd.to_numeric(sub['checkm_completeness'], errors='coerce')
    cont = pd.to_numeric(sub['checkm_contamination'], errors='coerce')
    sub['_comp'] = comp.fillna(float('-inf'))
    sub['_cont'] = cont.fillna(float('inf'))
    sub = sub.sort_values(by=['_comp','_cont'], ascending=[False, True])
    return sub.index[0]

def main():
    outdir = OUTDIR or os.path.dirname(os.path.abspath(ALIGNED_FASTA))
    os.makedirs(outdir, exist_ok=True)

    md = pd.read_csv(METADATA, sep='\t', header=0, index_col=0)

    # Filter metadata to domain if column exists
    if 'domain' in md.columns:
        md = md[md['domain'].astype(str).str.lower() == DOMAIN.lower()]

    # Optional renaming map
    if ID_MAP:
        id_map = read_id_map(ID_MAP)
        md.index = [id_map.get(idx, idx) for idx in md.index]

    # Harmonise IDs
    md.index = [i.replace('_genomic', '') for i in md.index]

    records = list(SeqIO.parse(ALIGNED_FASTA, 'fasta'))
    for r in records:
        if r.id.endswith('_genomic'):
            r.id = r.id[:-8]
    genes_16S = [r.id for r in records]
    genes_set = set(genes_16S)

    clmap = parse_clusters(CLUSTERS, genes_set)

    best_map = {}
    processed_rows = []
    for centroid, members in clmap.items():
        candidates = [centroid] + [m for m in members if m != centroid]
        ordered = []
        for c in candidates:
            if c not in ordered:
                ordered.append(c)
        best = choose_best(md, ordered)
        best_map[centroid] = best
        processed_rows.append((centroid, best, ','.join(ordered)))

    out_fa = os.path.splitext(ALIGNED_FASTA)[0] + '_best.fna'
    new_records = []
    out_ids = []
    for r in records:
        new_id = best_map.get(r.id, r.id)
        r.id = new_id
        r.description = ''
        new_records.append(r)
        out_ids.append(new_id)
    SeqIO.write(new_records, out_fa, 'fasta')

    proc_path = os.path.join(outdir, f"{DOMAIN}_16S_clusters_processed.txt")
    with open(proc_path, 'w') as f:
        f.write('Centroid\tBest\tAll genomes\n')
        for c, b, al in processed_rows:
            f.write(f"{c}\t{b}\t{al}\n")

    md_reduced = md.reindex(sorted(set(out_ids)))
    md_out = os.path.join(outdir, f"{DOMAIN}_metadata_clusters_ssu_align_centroids.csv")
    md_reduced.to_csv(md_out)

    print(f"[OK] Wrote: {out_fa}")
    print(f"[OK] Wrote: {proc_path}")
    print(f"[OK] Wrote: {md_out}")

if __name__ == '__main__':
    main()