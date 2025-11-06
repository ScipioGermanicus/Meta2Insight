#!/usr/bin/env python3
import os, re, sys, glob, shutil
import pandas as pd
from datetime import datetime
from typing import Optional

# ============================= CONFIG ========================================
INFILE = "intermediate/CheckM/merged/checkm_results.min.tsv"
OUTDIR = "intermediate/qc"
COMPLETENESS_MIN = 90.0
CONTAMINATION_MAX = 10.0
APPEND_GENOMIC_SUFFIX = False        # set True if you want IDs to gain "_genomic"
PLACE_MODE = "symlink"               # "symlink", "copy", or "hardlink"

# If you moved your MAGs, add/adjust candidates below; script picks the first that exists
# These are resolved relative to .../code/database_pipeline/
MAG_DIR_CANDIDATES = [
    os.path.join("intermediate", "genomes_to_search_barrnap"),
    os.path.join("intermediate", "MAGs_formatted", "genomes_to_search_barrnap"),
]

# File endings to consider when locating MAG files
EXTS = (".fa", ".fna", ".fasta")
GZ = ("", ".gz")

# ============================= LOGGING =======================================
os.makedirs(OUTDIR, exist_ok=True)
LOGFILE = os.path.join(OUTDIR, "log.txt")
def log(msg: str) -> None:
    ts = datetime.now().strftime("[%H:%M:%S]")
    print(f"{ts} {msg}")
    with open(LOGFILE, "a") as lf:
        lf.write(f"{ts} {msg}\n")

# ============================= PATH SETUP ====================================
BASE = os.path.dirname(os.path.abspath(__file__))                 # .../code/database_pipeline/script
ROOT = os.path.abspath(os.path.join(BASE, ".."))                  # .../code/database_pipeline
abs_in = os.path.join(ROOT, INFILE) if not os.path.isabs(INFILE) else INFILE
if not os.path.exists(abs_in):
    log(f"ERROR: INFILE not found: {abs_in}")
    sys.exit(1)

# pick genomes dir
GENOMES_DIR = None
for rel in MAG_DIR_CANDIDATES:
    cand = os.path.join(ROOT, rel)
    if os.path.isdir(cand):
        GENOMES_DIR = cand
        break
if GENOMES_DIR is None:
    log("ERROR: Could not find genomes directory in any of:")
    for rel in MAG_DIR_CANDIDATES:
        log("  - " + os.path.join(ROOT, rel))
    sys.exit(1)

# ============================= START =========================================
log("Starting CheckM quality filtering")
log(f"INFILE: {abs_in} (size: {os.path.getsize(abs_in)} bytes)")

# ============================= HELPERS =======================================
ARCHAEA_HINTS = re.compile(
    r"Archaea|Euryarchaeota|Thermoproteota|Crenarchaeota|Halobacteriota|Nanoarchaeota|Micrarchaeota|DPANN",
    re.I,
)
def guess_domain(marker_lineage: str) -> Optional[str]:
    if not isinstance(marker_lineage, str):
        return None
    return "Archaea" if ARCHAEA_HINTS.search(marker_lineage) else "Bacteria"

def normalize_id(gid: str) -> str:
    gid = gid.strip()
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", gid, flags=re.I)

def maybe_add_genomic_suffix(gid: str) -> str:
    if APPEND_GENOMIC_SUFFIX and not gid.endswith("_genomic"):
        gid = f"{gid}_genomic"
    return gid

# ============================= PARSE CHECKM ==================================
# Pass 1: try minimal 3-column block (genome_id, completeness, contamination)
min_rows = []
with open(abs_in, "r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.strip() or line.lstrip().startswith("-"):
            continue
        parts = [p.strip() for p in line.rstrip("\n").split("\t") if p.strip() != ""]
        if len(parts) == 3 and parts[0].lower() != "genome_id":
            try:
                cpl, cnt = float(parts[1]), float(parts[2])
                min_rows.append((normalize_id(parts[0]), cpl, cnt))
            except ValueError:
                pass
df_min = pd.DataFrame(min_rows, columns=["genome_id","checkm_completeness","checkm_contamination"]).drop_duplicates()
log(f"Parsed {len(df_min)} minimal entries")

# Pass 2: parse wide lines (duplicated panels); take first panel’s Completeness/Contamination
wide_rows = []
if df_min.empty:
    lead = re.compile(r"^\s*([A-Za-z0-9._+-]+)\s+([a-z]__.+?)\s{2,}", re.I)
    with open(abs_in, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith(("Bin Id","genome_id","-")):
                continue
            m = lead.match(raw)
            if not m:
                continue
            gid = normalize_id(m.group(1))
            lineage = (m.group(2) or "").strip()
            tail = raw[m.end():]
            tokens = [t for t in re.split(r"\s{2,}|\t", tail.strip()) if t]
            nums = []
            for t in tokens:
                tt = t.replace(",", "")
                try:
                    nums.append(float(tt))
                except ValueError:
                    break
            if len(nums) >= 5:
                completeness, contamination = nums[3], nums[4]
                wide_rows.append((gid, lineage, completeness, contamination))
    df_wide = pd.DataFrame(wide_rows, columns=["genome_id","marker_lineage","checkm_completeness","checkm_contamination"]) \
                .drop_duplicates(subset=["genome_id"])
    log(f"Extracted {len(df_wide)} wide entries")
else:
    df_wide = pd.DataFrame(columns=["genome_id","marker_lineage","checkm_completeness","checkm_contamination"])

if df_min.empty and df_wide.empty:
    log("No parsable entries found. Double-check file format.")
    sys.exit(1)

# Merge preference: use minimal for numbers if present, add lineage from wide
if not df_min.empty:
    df = df_min.merge(df_wide[["genome_id","marker_lineage"]], on="genome_id", how="left")
else:
    df = df_wide.copy()

df["domain"] = df["marker_lineage"].map(guess_domain).fillna("Bacteria")
df["checkm_completeness"] = pd.to_numeric(df["checkm_completeness"], errors="coerce")
df["checkm_contamination"] = pd.to_numeric(df["checkm_contamination"], errors="coerce")
df = df.dropna(subset=["genome_id","checkm_completeness","checkm_contamination"])

keep = df[
    (df["checkm_completeness"] >= COMPLETENESS_MIN) &
    (df["checkm_contamination"] <= CONTAMINATION_MAX)
].copy()

log(f"Kept {len(keep)} genomes ≥{COMPLETENESS_MIN}% completeness and ≤{CONTAMINATION_MAX}% contamination")

# ============================= WRITE TABLES/LISTS =============================
df.to_csv(os.path.join(OUTDIR, "checkm_clean_all.tsv"), sep="\t", index=False)
keep.to_csv(os.path.join(OUTDIR, "checkm_filtered.tsv"), sep="\t", index=False)

bac_ids = sorted(maybe_add_genomic_suffix(g) for g in keep.loc[keep["domain"]=="Bacteria","genome_id"].unique())
arc_ids = sorted(maybe_add_genomic_suffix(g) for g in keep.loc[keep["domain"]=="Archaea","genome_id"].unique())

with open(os.path.join(OUTDIR, "bacteria.txt"), "w") as f:
    f.write("\n".join(bac_ids) + ("\n" if bac_ids else ""))
with open(os.path.join(OUTDIR, "archaea.txt"), "w") as f:
    f.write("\n".join(arc_ids) + ("\n" if arc_ids else ""))

log(f"Wrote {len(bac_ids)} bacteria and {len(arc_ids)} archaea genome IDs")
log("Done.")

# ============================= GENOME DIR RESOLUTION ==========================
# Resolve GENOMES_DIR relative to repo root
for i, rel in enumerate(MAG_DIR_CANDIDATES):
    MAG_DIR_CANDIDATES[i] = os.path.join(ROOT, rel)
GENOMES_DIR = next((p for p in MAG_DIR_CANDIDATES if os.path.isdir(p)), None)
if GENOMES_DIR is None:
    log("ERROR: Could not find genomes directory after resolution.")
    sys.exit(1)

log(f"GENOMES_DIR: {os.path.abspath(GENOMES_DIR)}")
try:
    top = sorted(os.listdir(GENOMES_DIR))[:10]
    log(f"Top-level entries in GENOMES_DIR: {top}")
except Exception as e:
    log(f"WARNING: cannot list GENOMES_DIR: {e}")

# quick visibility test
test_hits = []
for ext in EXTS:
    for gz in GZ:
        test_hits += glob.glob(os.path.join(GENOMES_DIR, "**", f"*{ext}{gz}"), recursive=True)
log(f"Recursive visibility test: found {len(test_hits)} MAG files under GENOMES_DIR")

# ============================= FILE FINDER & PLACEMENT =======================
def find_genome_file(genome_id: str) -> Optional[str]:
    # Build candidate prefixes (with/without '_genomic')
    cands = {genome_id, re.sub(r"_genomic$", "", genome_id)}
    cands.add(genome_id if genome_id.endswith("_genomic") else genome_id + "_genomic")
    # Try exact-start and permissive-start patterns
    for prefix in cands:
        for ext in EXTS:
            for gz in GZ:
                pat = os.path.join(GENOMES_DIR, "**", f"{prefix}{ext}{gz}")
                hits = glob.glob(pat, recursive=True)
                if hits:
                    return hits[0]
                pat = os.path.join(GENOMES_DIR, "**", f"{prefix}*{ext}{gz}")
                hits = glob.glob(pat, recursive=True)
                if hits:
                    return hits[0]
    return None

# Placement destinations
bins_root = os.path.join(OUTDIR, "bins")
dst_bac = os.path.join(bins_root, "bacteria")
dst_arc = os.path.join(bins_root, "archaea")
os.makedirs(dst_bac, exist_ok=True)
os.makedirs(dst_arc, exist_ok=True)

log(f"Filtered IDs: bacteria={len(bac_ids)}, archaea={len(arc_ids)}")
probe = (bac_ids[:3] + arc_ids[:3])[:5]
for gid in probe:
    p = find_genome_file(gid)
    log(f"Probe resolve: {gid} -> {p}")

def place_file(src: str, dst: str) -> None:
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        return
    try:
        if PLACE_MODE == "symlink":
            os.symlink(src, dst)
        elif PLACE_MODE == "hardlink":
            os.link(src, dst)
        elif PLACE_MODE == "copy":
            shutil.copy2(src, dst)
        else:
            raise ValueError(f"Unknown PLACE_MODE: {PLACE_MODE}")
    except OSError as e:
        if PLACE_MODE in ("symlink", "hardlink"):
            log(f"Link failed ({e}); falling back to copy for {dst}")
            shutil.copy2(src, dst)
        else:
            raise

missing = []

def place_ids(id_list, dst_dir) -> int:
    placed = 0
    for gid in id_list:
        src = find_genome_file(gid)
        if not src:
            missing.append(gid)
            continue
        dst = os.path.join(dst_dir, os.path.basename(src))
        place_file(src, dst)
        placed += 1
    return placed

placed_bac = place_ids(bac_ids, dst_bac)
placed_arc = place_ids(arc_ids, dst_arc)

log(f"Placed {placed_bac} bacterial MAGs into {dst_bac} ({PLACE_MODE})")
log(f"Placed {placed_arc} archaeal MAGs into {dst_arc} ({PLACE_MODE})")
if missing:
    log(f"WARNING: {len(missing)} IDs had no matching file in GENOMES_DIR; first few: {missing[:5]}")

try:
    log(f"Example files in bacteria/: {sorted(os.listdir(dst_bac))[:3]}")
    log(f"Example files in archaea/: {sorted(os.listdir(dst_arc))[:3]}")
except Exception as e:
    log(f"Could not list destination dirs: {e}")
