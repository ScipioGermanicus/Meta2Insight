#!/usr/bin/env python3


# Counts copies per genome and gets files with multiple 16S copies only

# check file names and file paths


from Bio import SeqIO
import os
import pandas as pd
import numpy as np

kingdom = 'archaea' #note that this needs to be run twice - once with 'archaea' and once with 'bacteria'
genomes = os.listdir('barrnap_'+kingdom)

copies = []
for genome in genomes:
  count = 0
  sequences = []
  for record in SeqIO.parse('barrnap_'+kingdom+'/'+genome, "fasta"):
    if '16S' in record.id:
      sequences.append(record)
      count += 1
  copies.append([genome, count])
  if sequences != []:
    if len(sequences) == 1:
      w = SeqIO.write(sequences, 'barrnap_'+kingdom+'_16S_single/'+genome, "fasta")
    else:
      w = SeqIO.write(sequences, 'barrnap_'+kingdom+'_16S_multiple/'+genome, "fasta")

with open(kingdom+'_16S_copies.txt', 'w') as f:
  for copy in copies:
    w = f.write(copy[0].replace('_genomic.fna', '')+'\t'+str(copy[1])+'\n')
    
k08 = pd.read_csv(kingdom+'_16S_copies.txt', sep='\t', header=None, index_col=0)
k08 = k08[k08.max(axis=1) > 0]
print(kingdom)
print('80%: ', k08.shape[0], np.mean(k08.iloc[:, 0].values), np.median(k08.iloc[:, 0].values), np.max(k08.iloc[:, 0].values))
