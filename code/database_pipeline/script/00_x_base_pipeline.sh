#!/usr/bin/env bash
cd ~/Thesis/code/database_pipeline/


# 0_base_pipeline.sh
    # Sets up the database pipeline
    # Calls the other scripts in the correct order


# 1_formating.sh (FINISHED)
    # Fixes naming of files and compresses them
    # Splits them into bacteria and archaea folders based on GTDB-Tk output
bash script/01_a_formating_genomes.sh


# 2_taxonomy.sh (ALMOST FINISHED)
    # Fix issue: Perform on formatted MAGs, not unformatted ones
    # Overall done - needs specification of paths and testing
    # Runs GTDB-Tk to assign taxonomy
    # Outputs are stored in 
J1=$(sbatch --parsable script/02_a_infer_taxonomy.sh)


# 3_domain_classification.py (ALMOST FINISHED)
    # Needs only specification of paths based on GTDB-Tk output directories
    # Domain classification based on GTDB-Tk output
    # Splits genomes into bacteria and archaea folders
source ~/.bashrc
conda activate pyenv
python script/03_a_domain_classification.py


# 4_quality_control.sh (ALMOST FINISHED)
    # Evaluate after CheckM is done
    # Runs CheckM for quality metadata
J2=$(sbatch --parsable script/04_a_quality_control.sh)


# 5_quality_filtering.sh (FINISHED)
    # Check paths and filenames once all is in place
    # Uses CheckM1 output to filter genomes based on quality thresholds
source ~/.bashrc
conda activate pyenv
python script/05_a_quality_filtering.py

# 6_predict_16S.sh (IN PROGRESS)
    # Predicts 16S rRNA gene using Barrnap
    # Outputs are stored in separate folders for bacteria and archaea
J3=$(sbatch --parsable script/06_a_predict_16S.sh)


# 7_count_copies_per_genome (FINISHED)
    # Check paths once all is in place
    # Counts copies per genome and gets files with multiple 16S copies only
J4=$(sbatch --parsable script/07_a_count_copies_per_genome.sh)




# ...







