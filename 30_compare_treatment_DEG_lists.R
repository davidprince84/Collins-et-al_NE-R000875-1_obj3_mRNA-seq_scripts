#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Compare the differentially expressed genes (DEGs) between treatments 
# within tissues.
#-------------------------------------------------------------------------------
# Inputs:
# Lists of DEGs from the current study.

# Outputs:
# .csv file summarising the results of the statistical tests in a table.
# .csv file recording all the DEGs that overlap between the two treatments.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# Data loaded by CompareTreatmentDEGs function.

# FUNCTION DEFINITIONS ----

CompareTreatmentDEGs <- function(x, y) {
  # Compares the overlap in DEG lists between the M and H treatments
  # for a given tissue, and direction of differential expression.
  #
  # Args:
  #   x: string denoting from which tissue to compare DEG lists ("fatbody",
  #      "head" or "ovaries").
  #   y: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in RNA_2 compared to RNA_1) or 
  #      down-regulated genes (genes more highly expressed in RNA_1 compared to 
  #      RNA_2) ("up" or "down").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between the treatments. 
  # NOTE: The working directory needs to be set to the location of the DEG 
  # lists for the specified tissue.
  
  # Set variables based on arguments.
  
  mDEGs <- read.csv(paste0(x, "_results_M_treatment_", y, "_regulated_LFC0.csv"))
  
  hDEGs <- read.csv(paste0(x, "_results_H_treatment_", y, "_regulated_LFC0.csv"))
  
  backgroundList <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Direction_of_expression" = paste0(y, "-regulated"),
                          "Number_M_DEGs" = length(mDEGs$symbol),
                          "Number_H_DEGs" = length(hDEGs$symbol),
                          "Number_overlapping_DEGs" = NA,
                          "Number_DEGs_only_in_M" = NA,
                          "Number_DEGs_only_in_H" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = (0.05/3),
                          "Percentage_of_M_DEGs_overlapping" = NA) 
  
  # Calculate overlap between M and H treatments.
  
  # Determine overlapping genes.
  
  overlappingGenes <- mDEGs[mDEGs$symbol %in% hDEGs$symbol, ]
 
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$symbol)
  
  resultsDF$Percentage_of_M_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_M_DEGs)*100, digits = 1)
  
  # Prepare overlapping gene IDs for output.
  # Write results to data.frame, with extra blank space so that
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes[, c("symbol", "name")]
  
  colnames(genesToOutput) <- c(paste(x, y, "regulated_Dmel_gene_symbol", sep = "_"),
                               paste(x, y, "regulated_Dmel_gene_name", sep = "_"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:2000, ] <- ""
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_DEGs_only_in_M <- 
    resultsDF$Number_M_DEGs - resultsDF$Number_overlapping_DEGs
  
  resultsDF$Number_DEGs_only_in_H <-
    resultsDF$Number_H_DEGs - resultsDF$Number_overlapping_DEGs
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundList$gene_symbol) - resultsDF$Number_overlapping_DEGs - 
    resultsDF$Number_DEGs_only_in_M - resultsDF$Number_DEGs_only_in_H
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both DEG lists, only 
  # control DEG list, only removal DEG list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_DEGs, 
                               resultsDF$Number_DEGs_only_in_M,
                               resultsDF$Number_DEGs_only_in_H, 
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
  
  # Make and return the list of outputs.
  
  listToReturn <- list("genes" = genesToOutput, 
                       "stats" = resultsDF)
  
  return(listToReturn)
  
}

# EXECUTED STATEMENTS ----

# Head DEG Overlap ----

# Change working directory.
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

setwd("../02_outputs/02_head/31_DESeq2_DEG_lists/")

# Compare DEG lists.

headUpRegOverlap <- CompareTreatmentDEGs("head", "up")

headDownRegOverlap <- CompareTreatmentDEGs("head", "down")

# Fat Body DEG Overlap ----

# Change directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare DEG lists.

fatBodyUpRegOverlap <- CompareTreatmentDEGs("fatbody", "up")

fatBodyDownRegOverlap <- CompareTreatmentDEGs("fatbody", "down")

# Ovaries DEG Overlap ----

# Change directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists/")

# Compare DEG lists.

ovariesUpRegOverlap <- CompareTreatmentDEGs("ovaries", "up")

ovariesDownRegOverlap <- CompareTreatmentDEGs("ovaries", "down")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(headUpRegOverlap$stats,
                       headDownRegOverlap$stats,
                       fatBodyUpRegOverlap$stats,
                       fatBodyDownRegOverlap$stats,
                       ovariesUpRegOverlap$stats,
                       ovariesDownRegOverlap$stats)

# Combine overlapping genes results.

combinedGenes <- cbind(headUpRegOverlap$genes,
                       headDownRegOverlap$genes,
                       fatBodyUpRegOverlap$genes,
                       fatBodyDownRegOverlap$genes,
                       ovariesUpRegOverlap$genes,
                       ovariesDownRegOverlap$genes)

# Save Results ----

setwd("../../00_all_tissues/")

write.csv(combinedStats, "11_NER0008751_obj3_exp3_table_S6_comparison_stats.csv", 
          row.names = FALSE)

write.csv(combinedGenes, "12_NER0008751_obj3_exp3_table_D4_overlapping_genes.csv", 
          row.names = FALSE)
