#!/usr/bin/env bash
cd ~/Thesis/code/database_pipeline/


# 00_base_pipeline.sh
    # Sets up the database pipeline
    # Calls the other scripts in the correct order


# 01_formating.sh (FINISHED)
    # Fixes naming of files and compresses them
    # Splits them into bacteria and archaea folders based on GTDB-Tk output
bash script/01_a_formating_genomes.sh


# 02_taxonomy.sh (ALMOST FINISHED)
    # Fix issue: Perform on formatted MAGs, not unformatted ones
    # Overall done - needs specification of paths and testing
    # Runs GTDB-Tk to assign taxonomy
    # Outputs are stored in 
J1=$(sbatch --parsable script/02_a_infer_taxonomy.sh)


# 03_domain_classification.py (ALMOST FINISHED)
    # Needs only specification of paths based on GTDB-Tk output directories
    # Domain classification based on GTDB-Tk output
    # Splits genomes into bacteria and archaea folders
source ~/.bashrc
conda activate pyenv
python script/03_a_domain_classification.py


# 04_quality_control.sh (ALMOST FINISHED)
    # Evaluate after CheckM is done
    # Runs CheckM for quality metadata
J2=$(sbatch --parsable script/04_a_quality_control.sh)


# 05_quality_filtering.sh (FINISHED)
    # Check paths and filenames once all is in place
    # Uses CheckM1 output to filter genomes based on quality thresholds
source ~/.bashrc
conda activate pyenv
python script/05_a_quality_filtering.py

# 06_predict_16S.sh (FINISHED)
    # Check paths once all is in place
    # Predicts 16S rRNA gene using Barrnap
    # Outputs are stored in separate folders for bacteria and archaea
    # >90% completion and <10% redundancy
J3=$(sbatch --parsable script/06_a_predict_16S.sh)


# 07_count_copies_per_genome (FINISHED)
    # Check paths once all is in place
    # Counts copies per genome and gets files with multiple 16S copies only
J4=$(sbatch --parsable script/07_a_count_copies_per_genome.sh)


# 08_cluster_multiple_copies.py (FINISHED)
    # Check paths once all is in place
    # For each MAG in "*_16S_multiple/", it clusters the sequences within that file at --id 0.90 and writes centroids as <MAG>.fna into "*_16S_clustered/""
    # For "*_16S_single/", it just copies the single 16S FASTA into the same clustered folder
J5=$(sbatch --parsable script/08_cluster_multiple_copies.sh)


# 09_single_16S_per_genome.py (FINISHED)
    # Check paths and file names
    # Picks the longest 16S (if multiple), creates one combined FASTA per domain, writes sequence into that FASTA
source ~/.bashrc
conda activate pyenv
python script/09_single_16S_per_genome.py

# 10_cluster_single_16S_genes (FINISHED)
    # Takes final per-MAG 16S seqs and removes exact duplicates
J5=$(sbatch --parsable script/10_cluster_single_16S_genes.sh)


# 11_align_sequences

# ...







