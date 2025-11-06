#!/usr/bin/env python3
import re, os, sys, pandas as pd
from datetime import datetime
from typing import Optional

# --- CONFIG ---
INFILE = "intermediate/CheckM/merged/checkm_results.min.tsv"  # <-- singular if that's the real name
OUTDIR = "intermediate/qc"
COMPLETENESS_MIN = 90.0
CONTAMINATION_MAX = 10.0
APPEND_GENOMIC_SUFFIX = False  # set True if you want "<id>_genomic"

# --- IO setup ---
os.makedirs(OUTDIR, exist_ok=True)
LOGFILE = os.path.join(OUTDIR, "log.txt")
def log(msg):
    ts = datetime.now().strftime("[%H:%M:%S]")
    print(f"{ts} {msg}")
    with open(LOGFILE, "a") as lf: lf.write(f"{ts} {msg}\n")

abs_in = os.path.abspath(INFILE)
if not os.path.exists(INFILE):
    print(f"[ERROR] INFILE not found: {abs_in}")
    sys.exit(1)
log(f"Starting CheckM quality filtering")
log(f"INFILE: {abs_in} (size: {os.path.getsize(INFILE)} bytes)")

# --- Helpers ---
ARCHAEA_HINTS = re.compile(
    r"Archaea|Euryarchaeota|Thermoproteota|Crenarchaeota|Halobacteriota|Nanoarchaeota|Micrarchaeota|DPANN",
    re.I,
)
def guess_domain(marker_lineage: str) -> Optional[str]:
    if not isinstance(marker_lineage, str): return None
    return "Archaea" if ARCHAEA_HINTS.search(marker_lineage) else "Bacteria"

def normalize_id(gid: str) -> str:
    gid = gid.strip()
    gid = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", gid, flags=re.I)
    return gid

def maybe_add_genomic_suffix(gid: str) -> str:
    if APPEND_GENOMIC_SUFFIX and not gid.endswith("_genomic"):
        gid = f"{gid}_genomic"
    return gid

# --- Pass 1: try minimal 3-col table ---
min_rows = []
with open(INFILE, "r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.strip() or line.lstrip().startswith("-"): continue
        parts = [p.strip() for p in line.rstrip("\n").split("\t") if p.strip() != ""]
        if len(parts) == 3 and parts[0].lower() != "genome_id":
            try:
                cpl, cnt = float(parts[1]), float(parts[2])
                min_rows.append((normalize_id(parts[0]), cpl, cnt))
            except ValueError:
                pass
df_min = pd.DataFrame(min_rows, columns=["genome_id","checkm_completeness","checkm_contamination"]).drop_duplicates()
log(f"Parsed {len(df_min)} minimal entries")

# --- Pass 2: parse wide lines (handles space-only repeats, digits in lineage) ---
wide_rows = []
if df_min.empty:
    # capture: <gid> <lineage> then the rest; allow digits/() in lineage; stop at 2+ spaces
    lead = re.compile(r"^\s*([A-Za-z0-9._+-]+)\s+([a-z]__.+?)\s{2,}", re.I)
    with open(INFILE, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            s = raw.strip()
            if not s: continue
            if s.startswith("Bin Id") or s.startswith("genome_id") or s.startswith("-"): continue
            m = lead.match(raw)
            if not m: continue
            gid = normalize_id(m.group(1))
            lineage = (m.group(2) or "").strip()
            tail = raw[m.end():]

            # Take the first "panel": up to the point where the header would repeat.
            # Heuristic: split by 2+ spaces, keep numeric-looking tokens; take first ~20 tokens.
            tokens = [t for t in re.split(r"\s{2,}|\t", tail.strip()) if t]
            # The first panel should start with numbers: #genomes, #markers, #marker_sets, Completeness, Contamination, ...
            # Convert piecewise to floats, skipping non-numerics
            nums = []
            for t in tokens:
                tt = t.replace(",", "")  # be safe with thousands separators
                try:
                    nums.append(float(tt))
                except ValueError:
                    # first non-numeric in panel probably indicates column headers repeating; stop here
                    break

            if len(nums) >= 5:
                completeness = nums[3]
                contamination = nums[4]
                wide_rows.append((gid, lineage, completeness, contamination))

    df_wide = pd.DataFrame(
        wide_rows, columns=["genome_id","marker_lineage","checkm_completeness","checkm_contamination"]
    ).drop_duplicates(subset=["genome_id"])
    log(f"Extracted {len(df_wide)} wide entries")
else:
    df_wide = pd.DataFrame(columns=["genome_id","marker_lineage","checkm_completeness","checkm_contamination"])

# --- Merge and filter ---
if df_min.empty and df_wide.empty:
    log("No parsable entries found. Double-check filename, working directory, and file format.")
    sys.exit(0)

df = df_min.copy()
if df.empty:
    df = df_wide.copy()
else:
    df = df_min.merge(df_wide[["genome_id","marker_lineage"]], on="genome_id", how="left")

df["domain"] = df["marker_lineage"].map(guess_domain).fillna("Bacteria")
df["checkm_completeness"] = pd.to_numeric(df["checkm_completeness"], errors="coerce")
df["checkm_contamination"] = pd.to_numeric(df["checkm_contamination"], errors="coerce")
df = df.dropna(subset=["genome_id","checkm_completeness","checkm_contamination"])

keep = df[
    (df["checkm_completeness"] >= COMPLETENESS_MIN) &
    (df["checkm_contamination"] <= CONTAMINATION_MAX)
].copy()

log(f"Kept {len(keep)} genomes ≥{COMPLETENESS_MIN}% completeness and ≤{CONTAMINATION_MAX}% contamination")

# --- Save outputs to OUTDIR ---
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
