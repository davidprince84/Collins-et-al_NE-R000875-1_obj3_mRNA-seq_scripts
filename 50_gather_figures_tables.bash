#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq and comparative analysis with other data.
# Tasks: Make a new directory and copy all figure and table files across.
#-------------------------------------------------------------------------------
# Inputs:
# Completed analysis (scripts #01 - #45)
#
# Outputs:
# A directory containing all the figures and tables from the analysis.
#-------------------------------------------------------------------------------

# STEP 1: MAKE THE NEW DIRECTORY ----
# NOTE: script assumes the working directory is 01_scripts.

# Make Directory ----

mkdir ../02_outputs/20_all_figures_tables

# STEP 2: COPY RELEVANT FILES TO NEW DIRECTORY ----

# All Tissues ----

# Change directory.

cd ../02_outputs/00_all_tissues

# Copy figures.

cp *exp3_fig* ../20_all_figures_tables

# Copy tables.

cp *exp3_table* ../20_all_figures_tables

# Ovaries ----

# Change directory.

cd ../01_ovaries/02_quality_control_reports

# Copy file.

cp *exp3_supplementary_file* ../../20_all_figures_tables

# Copy figure.

cp *exp3_fig* ../../20_all_figures_tables

# Change directory.

cd ../31_DESeq2_DEG_lists

# Copy table. 

cp *exp3_table* ../../20_all_figures_tables

# Change directory.

cd ../32_DEG_top50_heatmaps

# Copy figure.

cp *exp3_fig* ../../20_all_figures_tables

# Head ----

# Change directory.

cd ../../02_head/02_quality_control_reports

# Copy file.

cp *exp3_supplementary_file* ../../20_all_figures_tables

# Copy figure.

cp *exp3_fig* ../../20_all_figures_tables

# Change directory.

cd ../31_DESeq2_DEG_lists

# Copy table. 

cp *exp3_table* ../../20_all_figures_tables

# Change directory.

cd ../32_DEG_top50_heatmaps

# Copy figure.

cp *exp3_fig* ../../20_all_figures_tables

# Fat body ----

# Change directory.

cd ../../03_fatbody/02_quality_control_reports

# Copy file.

cp *exp3_supplementary_file* ../../20_all_figures_tables

# Copy figure.

cp *exp3_fig* ../../20_all_figures_tables

# Change directory.

cd ../31_DESeq2_DEG_lists

# Copy table. 

cp *exp3_table* ../../20_all_figures_tables

# Change directory.

cd ../32_DEG_top50_heatmaps

# Copy figure.

cp *exp3_fig* ../../20_all_figures_tables

# STEP 3: REMOVE NUMBERS PRE-FIXING EACH FILE ----

# Change Directory ----

cd ../../20_all_figures_tables

# Remove Pre-fix Numbers ----

for FILE in *
do
	NEWFILE=${FILE:3}
mv $FILE $NEWFILE
done
