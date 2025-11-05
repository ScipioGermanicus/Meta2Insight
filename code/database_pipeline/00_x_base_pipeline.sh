#!/usr/bin/env bash

# 0_base_pipeline.sh
    # Sets up the database pipeline
    # Calls the other scripts in the correct order


# 1_formating.sh
    # Fixes naming of files and compresses them
    # Splits them into bacteria and archaea folders based on GTDB-Tk output

bash ~/Thesis/code/database_pipeline/01_a_formating_genomes.sh

# 2_taxonomy.sh
    # Runs GTDB-Tk to assign taxonomy
    # Outputs are stored in 


# 3_assigning.sh
    # Domain classification


# 4_qc.sh 
    # (Only if CheckM has not been applied on the data beforehand)
    # Runs CheckM to assess genome quality and obtain metadata


# 5_extracting.sh
    # Extracts 16S rRNA sequences using Barrnap
    # Outputs are stored in separate folders for bacteria and archaea


# 6_count_copies_per_genome
    # Counts copies per genome and gets files with multiple 16S copies only


# ...







