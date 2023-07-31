#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Make a series of directories/subdirectories to conduct the analysis in.
# Move/copy files to directories.
#-------------------------------------------------------------------------------
# Inputs:
# The contents of the GitHub repository 
# https://github.com/davidprince84/Collins-et-al_NE-R000875-1_obj3_mRNA-seq_scripts.git
# The mRNA-seq read files from GEO GSE175623 downloaded using sra. 
#
# Outputs:
# A series of named directories/subdirectories containing the scripts from
# the GitHub repository. Gene lists in the appropriate folders.
#-------------------------------------------------------------------------------

# STEP 1: SET UP THE DIRECTORY STRUCTURE BY MAKING DIRECTORIES ----

# NOTE: This script assumes that the working directory is the downloaded GitHub repository.

# Making New Directories ----

# Make the initial directory for the analysis containing the project name.

mkdir ../NER0008751_obj3_exp3_dmel

# Change into the new directory to make further directories.

cd ../NER0008751_obj3_exp3_dmel

# Make directories to store the data, scripts and outputs.

mkdir 00_data
mkdir 01_scripts
mkdir 02_outputs

# Change into the data directory to make further directories.

cd 00_data

# Make directories to store different types of data.

mkdir 00_fastq
mkdir 01_fasta
mkdir 02_gtf
mkdir 03_bed
mkdir 04_sums
mkdir 05_gene_list_csv
mkdir 10_transcripts2genes

# Change into the outputs directory to make further directories.

cd ../02_outputs

# Make directories to store outputs from analyses on all data, ovary, head and fat body samples.

mkdir 00_all_tissues
mkdir 01_ovaries
mkdir 02_head
mkdir 03_fatbody

# STEP 2: MOVE THE GITHUB REPOSITORY INTO THE SCRIPTS FOLDER ----

# Copy the Repository ----

# Change directory to the GitHub repository.

cd ../../Collins-et-al_NE-R000875-1_obj3_mRNA-seq_scripts

# Copy all the scripts and directories to the scripts folder.
# Args:
# -r: copy directories recursively

cp -r * ../NER0008751_obj3_exp3_dmel/01_scripts

# Copy the .git directory to the scripts folder (as the previous command does not do this).

cp -r .git ../NER0008751_obj3_exp3_dmel/01_scripts

# STEP 3: DELETE THE EMPTY GITHUB REPOSITORY ----

# Change Directory ----

cd ../

# Remove the Empty Repository ----
# Args:
# -r: remove directories and their contents recursively
# -f: do not ask permission to delete the files

rm -rf Collins-et-al_NE-R000875-1_obj3_mRNA-seq_scripts

# STEP 4: MOVE GENE LISTS TO APPROPRIATE DIRECTORY ----

mv 01_scripts/chen_2014_table_S1.csv 00_data/05_gene_list_csv/chen_2014_table_S1.csv

mv 01_scripts/pacifico_2018_table_S1_female.csv 00_data/05_gene_list_csv/pacifico_2018_table_S1_female.csv

mv 01_scripts/korb_2021_table_s1.csv 00_data/05_gene_list_csv/korb_2021_table_s1.csv

# STEP 5: MOVE THE MRNA-SEQ READS ----

# Change directory to where the mRNA-seq reads were downloaded to.

# Use the following (commented out) code below to move the reads to the 
# 00_data/00_fastq folder where PATH/TO/ is the path from the current 
# directory to NER0008751_obj3_dmel.

# mv *fastq.gz PATH/TO/NER0008751_obj3_dmel/00_data/00_fastq  # PLEASE CHANGE.
