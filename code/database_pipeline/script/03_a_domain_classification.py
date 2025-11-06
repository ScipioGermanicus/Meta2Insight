#!/usr/bin/env python3
import os, shutil, pandas as pd

# --- Resolve paths relative to this script ---
BASE = os.path.dirname(os.path.abspath(__file__))          # .../code/database_pipeline/script
ROOT = os.path.abspath(os.path.join(BASE, ".."))           # .../code/database_pipeline

# GTDB-Tk summaries
summaries = [
    os.path.join(ROOT, "intermediate", "GTDB-Tk", "gtdbtk.bac120.summary.tsv"),
    os.path.join(ROOT, "intermediate", "GTDB-Tk", "gtdbtk.ar53.summary.tsv"),
]

# Formatted genomes location (required)
FORMATTED_DIR = os.path.join(ROOT, "intermediate", "MAGs_formatted")
GENOME_DIR = FORMATTED_DIR
BAC_DIR = os.path.join(GENOME_DIR, "genomes_to_search_barrnap", "bacteria")
ARC_DIR = os.path.join(GENOME_DIR, "genomes_to_search_barrnap", "archaea")
for d in (GENOME_DIR, BAC_DIR, ARC_DIR):
    os.makedirs(d, exist_ok=True)

# --- Load rename map: user_genome -> formatted new_id (without suffix) ---
RENAME_MAP_PATH = os.path.join(FORMATTED_DIR, "id_map.tsv")
if not os.path.exists(RENAME_MAP_PATH):
    raise SystemExit(f"Expected rename map at {RENAME_MAP_PATH}. "
                     "It must map raw user_genome IDs to formatted basenames (column names like original_filename/new_id).")

rm = pd.read_csv(RENAME_MAP_PATH, sep="\t", dtype=str)
rm.columns = [c.strip().lower() for c in rm.columns]
old_candidates = [c for c in rm.columns if c in ("original_filename","original","raw_id","old_id","user_genome")]
new_candidates = [c for c in rm.columns if c in ("new_id","formatted_id","mag_id")]
if not old_candidates or not new_candidates:
    raise ValueError("id_map.tsv must have columns like original_filename/new_id (or raw_id/new_id).")
old_col, new_col = old_candidates[0], new_candidates[0]

def stem_no_suffix(x: str) -> str:
    b = os.path.basename(str(x))
    # strip common genome suffixes/extensions, keep bare stem
    for suf in ("_genomic.fna.gz","_genomic.fna",".fa.gz",".fna.gz",".fasta.gz",".fa",".fna",".fasta"):
        if b.endswith(suf):
            return b[:-len(suf)]
    return b.rsplit(".", 1)[0] if "." in b else b

rename_map = {stem_no_suffix(rm.at[i, old_col]): str(rm.at[i, new_col]) for i in rm.index}

# --- Read GTDB-Tk summaries and infer domain ---
dfs = []
for s in summaries:
    if os.path.exists(s):
        df = pd.read_csv(s, sep="\t")
        if {"user_genome","classification"}.issubset(df.columns):
            dfs.append(df[["user_genome","classification"]].copy())
        else:
            raise ValueError(f"{s} missing required columns: user_genome/classification")
if not dfs:
    raise SystemExit("No GTDB-Tk summary files found at expected paths.")

md = pd.concat(dfs, ignore_index=True).drop_duplicates("user_genome")

def infer_domain(cl):
    if isinstance(cl, str):
        if cl.startswith("d__Bacteria"): return "Bacteria"
        if cl.startswith("d__Archaea"):  return "Archaea"
    return "Unknown"

md["domain"] = md["classification"].map(infer_domain)
md = md[md["domain"] != "Unknown"].copy()

# --- Write lists next to formatted genomes ---
md[["user_genome","domain"]].to_csv(os.path.join(GENOME_DIR, "domain_map.tsv"), sep="\t", index=False)
md.query("domain=='Bacteria'")["user_genome"].to_csv(os.path.join(GENOME_DIR, "bacteria.txt"), index=False, header=False)
md.query("domain=='Archaea'")["user_genome"].to_csv(os.path.join(GENOME_DIR, "archaea.txt"), index=False, header=False)

# --- Move formatted files by domain (formatted-only, fail if missing) ---
SUFFIXES = ["_genomic.fna.gz","_genomic.fna",".fa.gz",".fa",".fna.gz",".fna",".fasta.gz",".fasta"]

def find_formatted_path(user_genome: str):
    # map to new_id if present, else assume user_genome is already the formatted basename
    bases = [rename_map.get(user_genome, user_genome)]
    for base in dict.fromkeys(bases):  # dedupe, preserve order
        for suf in SUFFIXES:
            p = os.path.join(FORMATTED_DIR, base + suf)
            if os.path.exists(p):
                return p
    return None

moved, missing = 0, []
for gid, dom in md[["user_genome","domain"]].itertuples(index=False):
    src = find_formatted_path(gid)
    if not src:
        missing.append(gid); continue
    dst_dir = BAC_DIR if dom == "Bacteria" else ARC_DIR
    shutil.move(src, os.path.join(dst_dir, os.path.basename(src)))
    moved += 1

print(f"Moved {moved} genomes.")
if missing:
    raise SystemExit(
        f"{len(missing)} genomes were not found in the formatted directory. "
        f"Check id_map.tsv and filenames. Examples: {missing[:10]}"
    )

def strip_suffix(name: str) -> str:
    for suf in ("_genomic.fna.gz","_genomic.fna",".fa.gz",".fna.gz",".fasta.gz",".fa",".fna",".fasta"):
        if name.endswith(suf):
            return name[:-len(suf)]
    return os.path.splitext(name)[0]

# regenerate lists from moved files to guarantee consistency
bac_files = sorted(os.listdir(BAC_DIR))
arc_files = sorted(os.listdir(ARC_DIR))

with open(os.path.join(GENOME_DIR, "bacteria.txt"), "w") as f:
    for fn in bac_files:
        if not fn.startswith("."):
            f.write(strip_suffix(fn) + "\n")

with open(os.path.join(GENOME_DIR, "archaea.txt"), "w") as f:
    for fn in arc_files:
        if not fn.startswith("."):
            f.write(strip_suffix(fn) + "\n")

