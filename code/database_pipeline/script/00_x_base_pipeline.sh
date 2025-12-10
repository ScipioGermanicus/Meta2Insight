#!/usr/bin/env bash
cd ~/Thesis/code/database_pipeline/


# 00_base_pipeline.sh
    # Calls scripts in correct order


# 01_formating.sh
    # Renames all MAG FASTA files to a consistent PICRUSt2-compatible pattern (e.g. MAG0001_genomic.fna.gz)
    # Compresses all genomes to gzipped .fna.gz and stores them in a clean output directory
    # Creates a mapping file linking original filenames to new MAG IDs for downstream tracking
bash script/01_a_formating_genomes.sh


# 02_taxonomy.sh
    ## Fix issue: Perform on formatted MAGs, not unformatted ones
    ## Overall done - needs specification of paths and testing
    # Runs GTDB-Tk taxonomy assignment on all formatted MAGs to classify them into bacterial or archaeal lineages
    # Uses the official GTDB-Tk reference package and writes full classification outputs to the GTDB-Tk directory
    # Produces summary files later used to split genomes by domain and continue the PICRUSt2 custom-DB pipeline
J1=$(sbatch --parsable script/02_a_infer_taxonomy.sh)


# 03_domain_classification.py
    ## Needs only specification of paths based on GTDB-Tk output directories
    # Uses GTDB-Tk summary files to assign each formatted MAG to Bacteria or Archaea based on its GTDB classification
    # Maps original MAG names to their new formatted IDs, moves the corresponding genome files into domain-specific folders, and writes updated bacteria/archaea ID lists
    # Produces a domain_map.tsv plus bacteria.txt and archaea.txt that define the genome sets used in the subsequent 16S gene prediction (barrnap) step
source ~/.bashrc
conda activate pyenv
python3 script/12_choose_best_genome.py \
  --clusters ~/Thesis/code/database_pipeline/intermediate/count_copies_per_genome/archaea_16S_clusters.uc \
  --aligned-fasta ~/Thesis/code/database_pipeline/intermediate/ssu_align/archaea_16S_centroids_ssu_align.fna \
  --metadata ~/Thesis/code/database_pipeline/intermediate/CheckM/merged/checkm_results_wide.tsv \
  --domain archaea


# 04_quality_control.sh
    ## Evaluate after CheckM is done
    # Runs CheckM lineage_wf on the bacteria and archaea MAG sets to estimate completeness and contamination for each genome
    # Normalises inputs to uncompressed .fna files, then generates per-domain QA tables and minimal TSVs with just genome_id, completeness, and contamination
    # Merges bacterial and archaeal results into a single summary file that will be used later to quality-filter genomes for the custom PICRUSt2 database
J2=$(sbatch --parsable script/04_a_quality_control.sh)


# 05_quality_filtering.sh
    ## Check paths and filenames once all is in place
    # Parses merged CheckM results, filters MAGs by completeness ≥90% and contamination ≤10%, and infers domain (Bacteria/Archaea) where possible
    # Writes full and filtered QC tables plus bacteria/archaea ID lists, then locates the corresponding MAG files on disk
    # Places only high-quality bacterial and archaeal MAGs into dedicated bins folders (via symlinks/copies), defining the QC-passed genome set for the custom PICRUSt2 database
source ~/.bashrc
conda activate pyenv
python script/05_a_quality_filtering.py

# 06_predict_16S.sh
    ## Check paths once all is in place
    # Runs Barrnap on all QC-passed bacterial and archaeal MAGs to predict rRNA (16S/SSU) loci using the appropriate kingdom models
    # For each genome, writes a GFF of rRNA hits and extracts the corresponding rRNA sequences into *_rRNA.fna files using bedtools
    # Produces per-domain Barrnap outputs that will feed into the subsequent 16S alignment and copy-number estimation steps of the custom PICRUSt2 pipeline
J3=$(sbatch --parsable script/06_a_predict_16S.sh)


# 07_count_copies_per_genome
    # Check paths once all is in place
    # Counts copies per genome and gets files with multiple 16S copies only
J4=$(sbatch --parsable script/07_a_count_copies_per_genome.sh)


# 08_cluster_multiple_copies.py
    # Check paths once all is in place
    # For each MAG in "*_16S_multiple/", it clusters the sequences within that file at --id 0.90 and writes centroids as <MAG>.fna into "*_16S_clustered/""
    # For "*_16S_single/", it just copies the single 16S FASTA into the same clustered folder
J5=$(sbatch --parsable script/08_cluster_multiple_copies.sh)


# 09_single_16S_per_genome.py
    # Check paths and file names
    # Picks the longest 16S (if multiple), creates one combined FASTA per domain, writes sequence into that FASTA
source ~/.bashrc
conda activate pyenv
python script/09_single_16S_per_genome.py

# 10_cluster_single_16S_genes
    # Takes final per-MAG 16S seqs and removes exact duplicates
J6=$(sbatch --parsable script/10_cluster_single_16S_genes.sh)


# 11_align_16S.sh
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
    # ARCHAEA later
    # Converts FASTA to DNA FASTA and Stockholm with HMMER/Easel tools, builds HMM


# 20_raxml_info_files.py
    # remember archaea later


# 21_formatting_raxml_info_files.sh
    # remember archaea later

# 22_filter_16S_copies_bac.py


# 23_annotation.sh
    # remember archaea later

# 24_build_kotable.py
python 24_build_kotable.py \
  --bac-dir /home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/eggnog_out \
  --bac-out /home/student.aau.dk/yr42on/Thesis/code/database_pipeline/intermediate/ko.txt.gz



# This is what it should look like when also using archaea (adjust paths)
source ~/.bashrc
conda activate barrnap_env
python 24_build_kotable.py \
  --bac-dir /home/student.aau.dk/yr42on/Thesis/data/database/eggnog/bacteria \
  --bac-out /home/student.aau.dk/yr42on/Thesis/picrust_db/default_files/bacteria/ko.txt.gz \
  --arc-dir /home/student.aau.dk/yr42on/Thesis/data/database/eggnog/archaea \
  --arc-out /home/student.aau.dk/yr42on/Thesis/picrust_db/default_files/archaea/ko.txt.gz







