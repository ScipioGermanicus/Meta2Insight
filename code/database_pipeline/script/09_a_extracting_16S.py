#!/usr/bin/env python3

# Python script for filtering: only include genomes with >90% completion and <10% redundancy



# tweak file names and paths


import pandas as pd
import os

barrnap_bacteria = os.listdir('intermediate/qc/bins/bacteria')
barrnap_archaea = os.listdir('intermediate/qc/bins/archaea')

bac_md = pd.read_csv('bac120_metadata_r220.tsv', index_col=0, header=0, sep='\t')
arc_md = pd.read_csv('ar53_metadata_r220.tsv', index_col=0, header=0, sep='\t')
md = pd.concat([bac_md, arc_md])
md = md[md['checkm_completeness'] >= 90]
md = md[md['checkm_contamination'] <= 10]
genomes = [g.replace('RS_', '').replace('GB_', '') for g in md.index.values]
barrnap_bacteria = [g.replace('_genomic.fna.gz', '') for g in barrnap_bacteria]
barrnap_archaea = [g.replace('_genomic.fna.gz', '') for g in barrnap_archaea]
genomes = set(genomes)
barrnap_bacteria = set(barrnap_bacteria)
barrnap_archaea = set(barrnap_archaea)
bac_count = 0
for genome in barrnap_bacteria:
  if genome not in genomes:
    m = os.system('mv genomes_to_search_barrnap/bacteria/'+genome+'_genomic.fna.gz gtdb_genomes/')
  else:
    bac_count += 1


arc_count = 0
for genome in barrnap_archaea:
  if genome not in genomes:
    m = os.system('mv genomes_to_search_barrnap/archaea/'+genome+'_genomic.fna.gz gtdb_genomes/')
  else:
    arc_count += 1
