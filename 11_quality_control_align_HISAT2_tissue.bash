#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Align mRNA-seq reads to Drosophila melanogaster (Dmel) genome using HISAT2.
#-------------------------------------------------------------------------------
# Inputs:
# HISAT2 index of the D. melanogaster genome, raw mRNA-seq read files.
#
# Outputs:
# .txt files summarizing the stats of each alignment command. Sorted .bam files
# of the reads. Indexes of .bam files. MultiQC report of HISAT2 alignment statistics.
#-------------------------------------------------------------------------------

# LOAD MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add HISAT2/2.1.0
module add python/anaconda/2019.10/3.7  # where MultiQC is installed.
module add samtools/1.10

# STEP 1: SET VARIABLES BASED ON ARGUMENTS ----

# $1 takes the value of the first argument after the bash script name when the script is executed.

if [ ${1} = ovaries ]
then
	NUMBER=01
	TISSUE=ovary
elif [ ${1} = head ]
then 
	NUMBER=02
	TISSUE=head
elif [ ${1} = fatbody ]
then 
	NUMBER=03
	TISSUE=fatbody
else
	echo "Arguments incorrect, must be ovaries, head or fatbody"
fi

# STEP 2: MAKE OUTPUT DIRECTORIES FOR RESULTS ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

mkdir ../02_outputs/${NUMBER}_${1}/10_hisat2_aligned_reads

mkdir ../02_outputs/${NUMBER}_${1}/11_hisat2_alignment_summaries

# STEP 3: ALIGN MRNA-SEQ TO DROSOPHILA MELANOGASTER GENOME ----

# Change Directory ----

cd ../00_data/00_fastq

# HISAT2 Alignment ----

# HISAT2 arguments: 
# -x: the basename of the index of the reference genome.
# --summary-file fileName and --new-summary: Print alignment summary in machine-friendly style to fileName.
# --dta: Tailor alignment report for Downstream Transcriptome Assembly, in this case with StringTie.
# -1: name of file containing read 1.
# -2: name of file containing read 2.

# Samtools arguments:
# -o fileName: write the final, sorted output to fileName.

for FILE in *${TISSUE}*R1.fastq.gz 
do
if [ ${1} = head ]
then
	OUTPUT=${FILE:24:21}
	FILEHANDLE=${FILE:0:45}
elif [ ${1} = fatbody ]
then 
	OUTPUT=${FILE:24:24}
	FILEHANDLE=${FILE:0:48}
else [ ${1} = ovaries ]
	OUTPUT=${FILE:24:21}
	FILEHANDLE=${FILE:0:46}
fi
hisat2 -x ../01_fasta/dmel_HISAT2_index --summary-file ../../02_outputs/${NUMBER}_${1}/11_hisat2_alignment_summaries/${OUTPUT}.txt --new-summary --dta -1 ${FILE} -2 ${FILEHANDLE}_R2.fastq.gz | samtools sort -o ../../02_outputs/${NUMBER}_${1}/10_hisat2_aligned_reads/${FILEHANDLE}_hisat2.bam
done

# STEP 4: GENERATE MULTIQC REPORT OF ALIGNMENT SUMMARIES ----

# Change Directory ----

cd ../../02_outputs/${NUMBER}_${1}/11_hisat2_alignment_summaries

# Generate Report ----
# Args:
# -n: rename report
# -o: save report in directory 

multiqc . -n 10_NER0008751_obj3_exp3_${1}_hisat2_multiqc.html -o ../02_quality_control_reports

# STEP 5: INDEX ALIGNMENT FILES ----

# Change Directory ----

cd ../10_hisat2_aligned_reads

# Generate Indexes ----
# Loop through all .bam files in the directory and generate an index.

for BAMFILE in *hisat2.bam
do
samtools index ${BAMFILE}
done
