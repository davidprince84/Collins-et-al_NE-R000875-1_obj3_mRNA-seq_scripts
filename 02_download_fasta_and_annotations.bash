#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Download required files from Ensembl.
#-------------------------------------------------------------------------------
# Inputs:
# None
#
# Outputs:
# Drosophila melanogaster (Dmel) genome and transcriptome files as zipped 
# fasta files, Dmel annotations as zipped GTF file and checksums for all files.
#-------------------------------------------------------------------------------

# STEP 1: DOWNLOAD DMEL FASTA FILES ----
# NOTE: script assumes the working directory is 01_scripts.

# Change Directory ----

cd ../00_data/01_fasta

# Download Genome Fasta File and Checksums ----

wget ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz

# Args:
# -O: Name of output file.

wget ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/dna/CHECKSUMS -O ../04_sums/dmel_BDGP6_22_97_dna_toplevel_checksums.txt

# Download Transcriptome Fasta File and Checksums ----

wget ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz

wget ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/cdna/CHECKSUMS -O ../04_sums/dmel_BDGP6_22_97_cdna_all_checksums.txt

# STEP 2: DOWNLOAD DMEL ANNOTATION FILES ----

# Download GTF File and Checksums ----

# Change directory.

cd ../02_gtf

# Download files.

wget ftp://ftp.ensembl.org/pub/release-97/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.97.gtf.gz

wget ftp://ftp.ensembl.org/pub/release-97/gtf/drosophila_melanogaster/CHECKSUMS -O ../04_sums/dmel_BDGP6_22_97_gtf_checksums.txt
