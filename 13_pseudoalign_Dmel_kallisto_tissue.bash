#!/bin/bash

#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Pseudoalign mRNA-seq reads to Drosophila melanogaster (Dmel) transcriptome 
# using Kallisto. 
#-------------------------------------------------------------------------------
# Inputs:
# Raw mRNA-seq files. Dmel transcriptome index for Kallisto.
#
# Outputs:
# Directories containing results of pseudoalignment (abundances).
# MultiQC summary of pseudoalignments.
# Text file containing study design table.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add kallisto/0.46.1
module add python/anaconda/2019.10/3.7  # Where MultiQC is installed.

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

# STEP 2: MAKE DIRECTORIES FOR KALLISTO OUTPUTs ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

mkdir ../02_outputs/${NUMBER}_${1}/20_kallisto_pseudoalignment_abundances

mkdir ../02_outputs/${NUMBER}_${1}/21_kallisto_pseudoalignment_summaries

# STEP 3: PSEUDOALIGN MRNA-SEQ READS USING KALLISTO ----

# Change Directory ----

cd ../00_data/00_fastq

# Pseudo Align with Kallisto ----

# Looping over mRNA-seq files to kallisto pseudoalign.
# Args:
# -i: fileName for index to be used for quanitification.
# -o: directory to write output to.
# --rf-stranded: strand specific reads, first read reverse.
# NOTE: 2> captures standard error to file, which is where the pseudoalignment stats are sent.

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
	OUTPUT=${FILE:24:22}
	FILEHANDLE=${FILE:0:46}
fi
kallisto quant -i ../01_fasta/dmel_cDNA_index.idx -o ../../02_outputs/${NUMBER}_${1}/20_kallisto_pseudoalignment_abundances/${OUTPUT} --rf-stranded ${FILE} ${FILEHANDLE}_R2.fastq.gz 2> ../../02_outputs/${NUMBER}_${1}/21_kallisto_pseudoalignment_summaries/${OUTPUT}.txt
done

# STEP 4: SUMMARISE PSEUDOALIGNMENTS WITH MULTIQC ----

# Change Directory ----

cd ../../02_outputs/${NUMBER}_${1}/21_kallisto_pseudoalignment_summaries/

# Run MultiQC ----
# -n fileName: call report fileName rather than default of multiqc_report.html
# -o: save report in directory 

multiqc . -n 12_NER0008751_obj3_exp3_${1}_kallisto_multiqc_report.html -o ../02_quality_control_reports

# STEP 5: GENERATE STUDY DESIGN FILE ----

# Initiate the file.
# Args:
# -e: enable interpretation of backslash escapes.

echo -e "sample\ttreatment\ttime_point" > ${1}_study_design.txt

# Loop over file names to fill in the details.

	for NAME in DC*.txt 
	do
if [ ${1} = head ]
then
	SAMPLE=${NAME:10:11}
	TREATMENT=${NAME:10:1}
	TIMEPOINT=${NAME:12:4}
elif [ ${1} = fatbody ]
then 
	SAMPLE=${NAME:13:11}
	TREATMENT=${NAME:13:1}
	TIMEPOINT=${NAME:15:4}
else 
	SAMPLE=${NAME:11:11}
	TREATMENT=${NAME:11:1}
	TIMEPOINT=${NAME:13:4}
fi
echo -e "${SAMPLE}\t${TREATMENT}\t${TIMEPOINT}" >> ${1}_study_design.txt
done
