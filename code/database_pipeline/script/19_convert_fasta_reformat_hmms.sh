#!/usr/bin/env bash

# FASTA -> DNA AFA (explicit nucleotides) and Stockholm
esl-reformat -d -o ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced_dna.fna afa \
  ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced.fna

esl-reformat -o ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced_dna.sto stockholm \
  ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced_dna.fna

# hmmbuild the 16S HMM
hmmbuild --cpu ${SLURM_CPUS_PER_TASK:-8} \
  ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced_dna.hmm \
  ~/Thesis/code/database_pipeline/intermediate/raxml/bacteria_raxml-check.raxml.reduced_dna.sto

# (Archaea later)
# esl-reformat -d -o .../archaea_raxml-check.raxml.reduced_dna.fna afa .../archaea_raxml-check.raxml.reduced.fna
# esl-reformat -o .../archaea_raxml-check.raxml.reduced_dna.sto stockholm .../archaea_raxml-check.raxml.reduced_dna.fna
# hmmbuild --cpu ... .../archaea_raxml-check.raxml.reduced_dna.hmm .../archaea_raxml-check.raxml.reduced_dna.sto
