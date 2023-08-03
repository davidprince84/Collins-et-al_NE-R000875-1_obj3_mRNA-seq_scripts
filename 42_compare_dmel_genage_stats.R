#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: Comparative analysis with other data.
# Tasks: Statistical comparison of overlap between current study differentially 
# expressed genes (DEGs) with Drosophila melanogaster (Dmel) GenAge data set.
#-------------------------------------------------------------------------------
# Inputs:
# GenAge Dmel gene list (from script #40). Lists of current study DEGs.

# Outputs:
# .csv file summarising the results of the statistical tests in a table.
# .csv file recording all the DEGs that overlap between current study and 
# GenAge. 
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# org.Dm.eg.db loaded via the script in the "Load data" section.

# LOAD DATA ----

# Load GenAge Data ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("40_compare_dmel_genage_prep.R")

# Remove unnecessary data.frames.

rm(chenDownDEGsFB, chenUpDEGsFB, pacificoUpDEGsFB, pacificoDownDEGsFB)

# FUNCTION DEFINITIONS ----

CompareDEGsAndDmelGenAge <- function (x, y) {
  # Compares the overlap between current study DEG lists current study and Dmel
  # genes in the GenAge database.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("fatbody", "head"  
  #      or "ovaries").
  #   y: string denoting the name of treatment from the current study to be 
  #      analysed ("M" or "H").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between the current 
  #      study and Dmel GenAge genes.
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Load current study lists.
  
  upDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                "up_regulated_LFC0.csv"))
  
  downDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                  "down_regulated_LFC0.csv"))
  
  backgroundList <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Reduce number of columns and rename column in background list.
  
  backgroundList <- as.data.frame(backgroundList[, 1])
  
  colnames(backgroundList) <- "Flybase_gene_ID"
  
  # Combine current study up- and down-regulated genes, simplify and rename column.
  
  currentStudyDEGs <- rbind(upDEGs, downDEGs)
  
  currentStudyDEGs <- as.data.frame(currentStudyDEGs[, c("symbol", "name")])
  
  colnames(currentStudyDEGs) <- c("Flybase_gene_ID", "Dmel_gene_name")
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Treatment" = y,
                          "Comparison_list" = "Dmel GenAge",
                          "Number_current_study_DEGs" = length(currentStudyDEGs$Flybase_gene_ID),
                          "Number_GenAge_genes" = length(genAgeGenes$symbol),
                          "Number_overlapping_DEGs" = NA,
                          "Number_DEGs_only_in_current_study" = NA,
                          "Number_genes_only_in_GenAge" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = (0.05/6),
                          "Percentage_of_current_study_DEGs_overlapping" = NA) 
  
  # Calculate overlap between current study and GenAge.
  
  # Determine overlapping genes.
  
  overlappingGenes <- currentStudyDEGs[currentStudyDEGs$Flybase_gene_ID %in% genAgeGenes$Flybase_gene_ID, ]
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$Flybase_gene_ID)
  
  resultsDF$Percentage_of_current_study_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_current_study_DEGs)*100, digits = 1)
  
  # Write results to data.frame, with extra blank space so that 
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes
  
  colnames(genesToOutput) <- c(paste0(x, "_", y, "_Flybase_gene_ID"), 
                               paste0(x, "_", y, "_Dmel_gene_name"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:300, ] <- ""
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_DEGs_only_in_current_study <- 
    resultsDF$Number_current_study_DEGs - resultsDF$Number_overlapping_DEGs
  
  resultsDF$Number_genes_only_in_GenAge <-
    resultsDF$Number_GenAge_genes - resultsDF$Number_overlapping_DEGs
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundList$Flybase_gene_ID) - resultsDF$Number_overlapping_DEGs - 
    resultsDF$Number_DEGs_only_in_current_study - resultsDF$Number_genes_only_in_GenAge
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both DEG lists, only 
  # M DEG list, only H DEG list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_DEGs, 
                               resultsDF$Number_DEGs_only_in_current_study,
                               resultsDF$Number_genes_only_in_GenAge, 
                               resultsDF$Number_genes_in_neither))
  
  # Change the dimensions to 2 rows and 2 columns.
  
  dim(contingencyTable) <- c(2,2)
  
  # Conduct two-tailed Fisher's Exact Test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance.
  
  fisherResults <- fisher.test(contingencyTable)
  
  # Add results of statistical tests to resultsDF.
  
  resultsDF$p_value <- fisherResults$p.value
  
  resultsDF$Odds_ratio <- fisherResults$estimate[[1]]
  
  # Return list of results.
  
  listToReturn <- list("genes" = genesToOutput, 
                       "stats" = resultsDF)
  
  return(listToReturn)
  
}

# EXECUTED STATEMENTS ----

# Head Comparison ----

# Set working directory.

setwd("../02_outputs/02_head/31_DESeq2_DEG_lists/")

# Compare gene lists.

headMResults <- CompareDEGsAndDmelGenAge("head", "M")

headHResults <- CompareDEGsAndDmelGenAge("head", "H")

# Fat Body Comparison ----

# Set working directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare gene lists.

fatBodyMResults <- CompareDEGsAndDmelGenAge("fatbody", "M")

fatBodyHResults <- CompareDEGsAndDmelGenAge("fatbody", "H")

# Ovaries Comparison ----

# Set working directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists/")

# Compare gene lists.

ovariesMResults <- CompareDEGsAndDmelGenAge("ovaries", "M")

ovariesHResults <- CompareDEGsAndDmelGenAge("ovaries", "H")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(headMResults$stats, headHResults$stats,
                       fatBodyMResults$stats, fatBodyHResults$stats,
                       ovariesMResults$stats, ovariesHResults$stats)

# Combine genes results.

combinedGenes <- cbind(headMResults$genes, headHResults$genes,
                       fatBodyMResults$genes, fatBodyHResults$genes,
                       ovariesMResults$genes, ovariesHResults$genes)

# Save Results ----

setwd("../../00_all_tissues/")

write.csv(combinedStats, "21_NER0008751_obj3_exp3_table_S9_GenAge_stats.csv",
          row.names = FALSE)

write.csv(combinedGenes, "22_NER0008751_obj3_exp3_table_D7_GenAge_overlapping_genes.csv",
          row.names = FALSE)
