#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Generating indexes.
# NOTE: Indexes and file conversions for programmes that require them are 
# generated in alphabetical order based on the name of the programme.
#-------------------------------------------------------------------------------
# Inputs:
# Drosophila melanogaster (Dmel) genome and transcriptome fasta files, and GTF 
# file. 
#
# Outputs:
# Indexes for Kallisto and HISAT2.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add HISAT2/2.1.0
module add kallisto/0.46.1

# STEP 1: UNZIP THE GENOME AND TRANSCRIPTOME FILES ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Change Directory ----

cd ../00_data/01_fasta

# Unzip Fasta Files ----

# Unzip genome fasta file.
# Args:
# --stdout: Keep file unchanged and write output to standard output.

gunzip --stdout Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz > dmel_BDGP6_22_97_dna_toplevel.fasta

# Unzip transcriptome fasta file.

gunzip --stdout Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz > dmel_BDGP6_22_97_cdna_all.fasta

# STEP 2: UNZIP THE GTF FILE ----

# Change Directory ----

cd ../02_gtf

# Unzip GTF File ----

gunzip --stdout Drosophila_melanogaster.BDGP6.22.97.gtf.gz > dmel_BDGP6_22_97.gtf

# STEP 3: GENERATE HISAT2 INDEX OF GENOME ----

# Change Directory ----

cd ../01_fasta

# Build HISAT2 Index ----

hisat2-build dmel_BDGP6_22_97_dna_toplevel.fasta dmel_HISAT2_index

# STEP 4: GENERATE KALLISTO INDEX OF TRANSCRIPTOME ----

# Indexing command.
# Args:
# -i: fileName for index to be created.

kallisto index -i dmel_cDNA_index.idx dmel_BDGP6_22_97_cdna_all.fasta
