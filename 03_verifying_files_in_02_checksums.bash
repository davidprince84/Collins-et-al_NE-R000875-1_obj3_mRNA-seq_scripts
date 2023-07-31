#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Select relevant entries from checksums files and compare with the
# checksums generated from the downloaded files.
#-------------------------------------------------------------------------------
# Inputs:
# Dmel fasta and annotation files, and associated checksum files, downloaded
# from Ensembl using script 02.
#
# Outputs:
# .txt file stating whether the checksums match between the downloaded files
# and the files present on Ensembl.
#-------------------------------------------------------------------------------

# STEP 1: SELECT RELEVANT CHECKSUMS ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Print Relevant Checksums to New Document for Comparison ----

# Change directory.

cd ../00_data/04_sums

# Print the relevant checksums.

# Genome checksums.

awk '/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz/' dmel_BDGP6_22_97_dna_toplevel_checksums.txt > original_ensembl_checksums.txt

# cDNA checksums. 

awk '/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz/' dmel_BDGP6_22_97_cdna_all_checksums.txt >> original_ensembl_checksums.txt

# GTF checksums.

awk '/Drosophila_melanogaster.BDGP6.22.97.gtf.gz/' dmel_BDGP6_22_97_gtf_checksums.txt >> original_ensembl_checksums.txt

# STEP 2: CONDUCT CHECKSUMS ON DOWNLOADED FILES ----
# Print the checksums results and the file name to a new document so that they can be compared.

# Genome Checksums ----

echo "$(sum ../01_fasta/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz) $(echo Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz)" > downloaded_ensembl_checksums.txt

# cDNA Checksums ----

echo "$(sum ../01_fasta/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz) $(echo Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz)" >> downloaded_ensembl_checksums.txt

# GTF Checksums ----

echo "$(sum ../02_gtf/Drosophila_melanogaster.BDGP6.22.97.gtf.gz) $(echo Drosophila_melanogaster.BDGP6.22.97.gtf.gz)" >> downloaded_ensembl_checksums.txt

# STEP 3: COMPARE THE ORIGINAL AND DOWNLOADED CHECKSUMS ----
# Check that the checksums of the downloaded file match those provided with the file
# by comparing the files.

# Make Empty File for Results Output ----

touch ../../02_outputs/00_all_tissues/01_checksums_ensembl_comparison_results.txt

# Compare Results ----

# grep arguments:
# -q: quiet, exits with zero status (mapped to true in the "if" statement) if any matches
# found, and exits with non-zero status (i.e. "false") if no matches found.
# -f fileName: obtains pattern(s) for comparison from fileName. 

while IFS= read -r FILE
do
if grep -q "$FILE" original_ensembl_checksums.txt
then
	RESULTS="$FILE \t Checksums match, therefore downloaded file is likely to be complete"
else 
	RESULTS="$FILE \t Checksums do not match, therefore downloaded file is likely to be incomplete"
fi
echo -e $RESULTS >> ../../02_outputs/00_all_tissues/01_checksums_ensembl_comparison_results.txt
done < "downloaded_ensembl_checksums.txt"
