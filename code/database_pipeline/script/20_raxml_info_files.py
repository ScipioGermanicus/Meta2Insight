#!/usr/bin/env python3
"""
Create a RAxML_info-style file for PICRUSt2/SEPP from a RAxML-NG --evaluate log.

Domain: bacteria only (archaea block is commented out for future use).

Inputs (hard-coded):
  - bacteria_raxml.raxml.log   (from raxml-ng --evaluate)
Outputs:
  - bacteria_raxml.raxml_info  (RAxML 7.x-style info file)
"""

import os
import re

# ======================
# CONFIG — EDIT IF NEEDED
# ======================
RAXML_DIR   = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/raxml')
BAC_PREFIX  = os.path.join(RAXML_DIR, 'bacteria_raxml')
BAC_LOG     = BAC_PREFIX + '.raxml.log'
BAC_INFO    = BAC_PREFIX + '.raxml_info'

# For later:
# ARC_PREFIX  = os.path.join(RAXML_DIR, 'archaea_raxml')
# ARC_LOG     = ARC_PREFIX + '.raxml.log'
# ARC_INFO    = ARC_PREFIX + '.raxml_info'
# ======================


def parse_raxml_ng_log(log_path: str):
    """Parse key values out of a RAxML-NG --evaluate log."""
    patterns = None
    base_freqs = None
    subs_rates = None
    final_ll = None
    elapsed = None
    called_line = None
    cmd_line = None

    with open(log_path, 'r') as f:
        lines = f.readlines()

    # Find "RAxML-NG was called at ..." and the next non-empty line (the command)
    for i, line in enumerate(lines):
        if "RAxML-NG was called at" in line:
            called_line = line.rstrip("\n")
            # next non-empty line
            for j in range(i + 1, len(lines)):
                cmd = lines[j].strip()
                if cmd:
                    cmd_line = cmd
                    break
            break

    for line in lines:
        line_stripped = line.strip()

        # Alignment patterns
        m = re.search(r"Alignment comprises\s+(\d+)\s+partitions\s+and\s+(\d+)\s+patterns", line_stripped)
        if m:
            patterns = int(m.group(2))
            continue

        # Base frequencies
        if line_stripped.startswith("Base frequencies (ML):"):
            # after colon: four floats
            parts = line_stripped.split(":", 1)[1].strip().split()
            base_freqs = parts  # ['0.241626', '0.217110', '0.322555', '0.218709']
            continue

        # Substitution rates
        if line_stripped.startswith("Substitution rates (ML):"):
            parts = line_stripped.split(":", 1)[1].strip().split()
            subs_rates = parts  # six floats
            continue

        # Final log-likelihood
        m2 = re.search(r"Final LogLikelihood:\s+(-?[0-9.eE+-]+)", line_stripped)
        if m2:
            final_ll = m2.group(1)
            continue

        # Elapsed time
        m3 = re.search(r"Elapsed time:\s+([0-9.eE+-]+)\s+seconds", line_stripped)
        if m3:
            elapsed = m3.group(1)
            continue

    if patterns is None:
        raise RuntimeError(f"Could not parse alignment patterns from {log_path}")
    if base_freqs is None:
        raise RuntimeError(f"Could not parse base frequencies from {log_path}")
    if subs_rates is None:
        raise RuntimeError(f"Could not parse substitution rates from {log_path}")
    if final_ll is None:
        raise RuntimeError(f"Could not parse final log-likelihood from {log_path}")
    if elapsed is None:
        # not strictly required, but nice to have; default to 0 if missing
        elapsed = "0.0"
    if called_line is None or cmd_line is None:
        raise RuntimeError(f"Could not find 'RAxML-NG was called at' section in {log_path}")

    return {
        "patterns": patterns,
        "base_freqs": base_freqs,
        "subs_rates": subs_rates,
        "final_ll": final_ll,
        "elapsed": elapsed,
        "called_line": called_line,
        "cmd_line": cmd_line,
    }


def make_raxml_info(log_path: str, info_path: str):
    """Build the RAxML_info-style file from a RAxML-NG log."""
    parsed = parse_raxml_ng_log(log_path)

    patterns   = parsed["patterns"]
    base_freqs = parsed["base_freqs"]
    subs_rates = parsed["subs_rates"]
    final_ll   = parsed["final_ll"]
    elapsed    = parsed["elapsed"]
    called     = parsed["called_line"]
    cmd        = parsed["cmd_line"]

    # NOTE:
    # - We hard-code alpha[0] = 1.000000 like the reference hack.
    # - pplacer only seems to care that the file parses, not the exact time.
    # - Base frequencies and rates are already in AC/AG/AT/CG/CT/GT order
    #   according to the RAxML documentation and the PICRUSt2 dev notes.

    txt = []
    txt.append("This is RAxML version 7.7.2 released by Alexandros Stamatakis on July 31 2013.")
    txt.append("")
    txt.append("This is a RAxML_info file from an --evaluate run, manually reformatted")
    txt.append("")
    txt.append("Partition: 0")
    txt.append(f"Alignment Patterns: {patterns}")
    txt.append("Name: No Name Provided")
    txt.append("DataType: DNA")
    txt.append("Substitution Matrix: GTR")
    txt.append("")
    txt.append(called)
    txt.append("")
    txt.append(cmd)
    txt.append("")
    txt.append("Base frequencies: " + " ".join(base_freqs))
    txt.append("")
    txt.append(f"Inference[0]: Time {elapsed} CAT-based likelihood -0000, best rearrangement setting 5")
    txt.append("alpha[0]: 1.000000 rates[0] ac ag at cg ct gt: " + " ".join(subs_rates))
    txt.append("")
    txt.append("")
    txt.append("NOT conducting any final model optimizations on all 1 trees under CAT-based")
    txt.append("model ....")
    txt.append("")
    txt.append(f"Final GAMMA  likelihood: {final_ll}")
    txt.append("")

    with open(info_path, "w") as out:
        out.write("\n".join(txt))

    print(f"[OK] Wrote RAxML_info: {info_path}")


def main():
    # Bacteria
    if not os.path.exists(BAC_LOG):
        raise SystemExit(f"[ERROR] Log file not found: {BAC_LOG}")
    make_raxml_info(BAC_LOG, BAC_INFO)

    # Archaea (enable later)
    # if os.path.exists(ARC_LOG):
    #     make_raxml_info(ARC_LOG, ARC_INFO)
    # else:
    #     print(f\"[WARN] Archaea log not found, skipping: {ARC_LOG}\")


if __name__ == "__main__":
    main()
