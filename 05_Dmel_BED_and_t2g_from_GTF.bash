#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Generate BED file from GTF for use with RSeQC. Isolate list of 
# transcripts and corresponding genes from GTF file.
#-------------------------------------------------------------------------------
# Inputs:
# Drosophila melanogaster (Dmel) GTF file.
#
# Outputs:
# Dmel BED file. Text document linking transcript and gene IDs.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module add python/anaconda/2019.10/3.7  # Where Mikado is installed.

# STEP 1: PRINT TRANSCRIPT AND GENE IDS FROM GTF TO A NEW FILE ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Change Directory ----

cd ../00_data/02_gtf

# Awk Columns from GTF ----

# Using awk to search the GTF.
# If the third column equals transcript then print the twelth column (transcript id) 
# and tenth column (gene id) separated by a tab.

awk '$3 == "transcript" {print $12"\t"$10}' dmel_BDGP6_22_97.gtf > transcripts2genes.txt

# STEP 2: EDIT THE FILE ----

# Removing Unwanted Characters ----

# Remove ; from the file using sed.

sed 's/;//g' transcripts2genes.txt > Dmel_6_22_97_transcripts2genes.txt

# STEP 3: MOVE THE FILE TO THE CORRECT SUBDIRECTORY ----

# Move File ----

mv Dmel_6_22_97_transcripts2genes.txt ../10_transcripts2genes

# Remove Intermediate File ----

rm -f transcripts2genes.txt

# STEP 4: GENERATING BED12 FILE FOR RSEQC ----

# Generating BED12 File ----

# BED12 file of annotations needed for RSeQC. Mikado is used to convert the Ensembl GTF file
# to a BED12 file.
# Args: 
# -of: specify output format {bed12, gtf or gff3}

mikado util convert -of bed12 dmel_BDGP6_22_97.gtf ../03_bed/Dmel_6_22_97.BED
