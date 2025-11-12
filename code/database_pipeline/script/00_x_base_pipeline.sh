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
python3 script/12_choose_best_genome.py \
  --clusters ~/Thesis/code/database_pipeline/intermediate/count_copies_per_genome/archaea_16S_clusters.uc \
  --aligned-fasta ~/Thesis/code/database_pipeline/intermediate/ssu_align/archaea_16S_centroids_ssu_align.fna \
  --metadata ~/Thesis/code/database_pipeline/intermediate/CheckM/merged/checkm_results_wide.tsv \
  --domain archaea


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
J6=$(sbatch --parsable script/10_cluster_single_16S_genes.sh)


# 11_align_sequences)
    # Aligns the 16S sequences 
    # Uses Infernal (cmalign) # Rfam models instead of ssualign (not available)
J7=$(sbatch --parsable script/11_align_16S.sh)


# 12_choose_best_genome_arc.py
    # Selects best genome per cluster using CheckM completeness->contamination
    # Rewrites aligned centroid FASTA IDs to the chosen best
source ~/.bashrc
conda activate barrnap_env
python script/12_choose_best_genome_arc.py


# 13_choose_best_genome_bac
    # Same thing as step 12, but for bacteria
source ~/.bashrc
conda activate barrnap_env
python script/13_choose_best_genome_bac.py
   

# 14_raxmlng_check.s (INCOMPLETE)
#--># Skip archaea for now - requires at least four genomes
    # Auto-detects alignment length, sequence problems, and drops duplicates, writing a cleaned PHYLIP alignment


# 15_convert_archaea_alignment (MISSING)
# COMPLETELY MISSING (ARCHAEAL-SPECIFIC)


# 16_build_trees.py (INCOMPLETE)
    # HAS TO BE EDITED FOR ARCHAEA TO BE INCLUDED
    # Builds trees using IQ-Tree based on given genomes
    # One for bacteria, one for archaea


# 17_run_raxml
    # UNCOMMENT ARCHAEA PART LATER


# 18_convert_phylip_to_fasta.py
    # Archaea
    # Converts the cleaned PHYLIP alignment from RAxML-NG back to FASTA format for downstream analyses
source ~/.bashrc
conda activate barrnap_env
python script/18_convert_phylip_to_fasta.py
   
# 19_convert_fasta_reformat_hmms.sh
    # ARCHAEA
    # Converts FASTA to DNA FASTA and Stockholm with HMMER/Easel tools, builds HMM












