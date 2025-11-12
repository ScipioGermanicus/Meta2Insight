#!/usr/bin/env python3
# save as: 20_phylip_to_fasta_bac.py
# Converts RAxML-NG *.raxml.reduced.phy -> FASTA for HMM prep.
# Archaea lines included but commented for later use.

import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ====== CONFIG ======
BAC_PHYLIP = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced.phy')
BAC_FASTA  = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced.fna')

# ARC_PHYLIP = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/raxml/archaea_raxml-check.raxml.reduced.phy')
# ARC_FASTA  = os.path.expanduser('~/Thesis/code/database_pipeline/intermediate/raxml/archaea_raxml-check.raxml.reduced.fna')
# ==========================================

def phylip_relaxed_to_fasta(phylip_path: str, fasta_path: str) -> None:
    seq_records = []
    with open(phylip_path, 'r') as f:
        header = f.readline()  # e.g. "N M"
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            # relaxed PHYLIP from raxml-ng: "<name><space><sequence>"
            # names have no spaces; sequence is everything after first space
            parts = line.split(' ', 1)
            if len(parts) < 2:
                continue  # skip malformed lines
            name, seq = parts[0], parts[1].replace(' ', '')
            # make a SeqRecord (no need for Bio.Alphabet)
            seq_records.append(SeqRecord(Seq(seq), id=name, description=''))
    msa = MultipleSeqAlignment(seq_records)
    AlignIO.write(msa, fasta_path, 'fasta')
    print(f"[OK] Wrote FASTA: {fasta_path}  (seqs: {len(seq_records)})")

def main():
    # Bacteria
    phylip_relaxed_to_fasta(BAC_PHYLIP, BAC_FASTA)
    # Archaea (enable later)
    # phylip_relaxed_to_fasta(ARC_PHYLIP, ARC_FASTA)

if __name__ == '__main__':
    main()
