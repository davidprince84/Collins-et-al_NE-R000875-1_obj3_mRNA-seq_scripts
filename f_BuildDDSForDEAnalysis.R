#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 12.01.2021
# File last updated: 09.06.2023
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Tasks: Summarise transcript-level counts to gene level and create DESeq2
# data set for differential expression analysis.
#-------------------------------------------------------------------------------
# Inputs:
# Abundances from Kallisto pseudoalignment. Table describing the experimental
# design.

# Outputs:
# dds (DESeq2 data set) object.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(tximport)  # tximport()
library(DESeq2)  # DESeqDataSetFromTximport(), DESeq() 

# LOADING DATA ----
# Note:
# This function requires the "Dmel_6_22_97_transcripts2genes.txt" file to be 
# loaded in the environment.

# FUNCTION DEFINITION ----

BuildDDSForDEAnalysis <- function (x) {
  # Summarises transcript-level counts to gene level and create DESeq2
  # data set (dds).
  #
  # Args:
  #   x: string denoting which tissue is to be analysed (fatbody", "head" 
  #      or "ovaries").
  #
  # Returns:
  #   A DESeq2 object with the transcript level counts summarised to gene level.
  
  # Set variables based on argument.
  
  if (x == "head") {
    tissuePath = "02_head"
  } else if (x == "fatbody") {
    tissuePath = "03_fatbody"
  } else if (x == "ovaries") {
    tissuePath = "01_ovaries"
  } else {
    stop('Argument x must be "fatbody", "head" or "ovaries".')
  }
  
  # Set common internal variable.
  
  baseDir <- #"path/to/NER0008751_obj3_exp3_dmel/02_outputs/"
  # This string is specific to your computer, PLEASE CHANGE accordingly.
  
  abundanceDir <- "20_kallisto_pseudoalignment_abundances"
  summaryDir <- "21_kallisto_pseudoalignment_summaries"
  study <- "_study_design.txt"
  
  # Generate character vector with sample names.
  
  sample_id <- 
    read.table(file.path(baseDir, tissuePath, summaryDir, 
                         paste0(x, study)),
               header = TRUE)
  
  # Add column containing full sample ID, to match the folder names.
  
  sample_id$full_id <- paste0(dir(file.path(baseDir, tissuePath, abundanceDir)))
  
  # Add condition column to sample_id.
  
  sample_id$condition <- paste0(sample_id$treatment, "_", sample_id$time_point) 
  
  # Generate named vector with path to quantification files.
  
  files <- file.path(baseDir, tissuePath, abundanceDir, sample_id$full_id, "abundance.h5")
  
  names(files) <- paste0(sample_id$sample)
  
  # Import transcript-level estimates and summarise to gene level (Bter only)
  # producing "original counts and offsets".
  
  txi <- tximport(files, type = "kallisto", tx2gene = t2g)  

  # Initiate DESeq2 data set (dds) object.
  
  ddsKallisto <- DESeqDataSetFromTximport(txi, sample_id, ~condition)
  # Model has condition as the main effect.
  
  # Filter rows of dds object to remove any genes with less than 10 counts across 
  # all samples.
  
  ddsKallisto <- ddsKallisto[rowSums(counts(ddsKallisto)) >= 10, ]
  
  # Conduct differential expression analysis.
  
  ddsKallisto <- DESeq(ddsKallisto)
  
  # Return ddsKallisto
  
  return(ddsKallisto)
  
}
