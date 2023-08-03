#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: Comparative analysis with other data.
# Tasks: Comparison of current study differentially expressed genes (DEGs) with 
# Drosophila melanogaster (Dmel) differentially expressed genes from Chen et al.
# (2014) and Pacifico et al. (2018).
#-------------------------------------------------------------------------------
# Inputs:
# Lists of Dmel DEGs from the current study. Lists of Dmel DEGs from Chen et al.
# (2014) and Pacifico et al. (2018). 

# Outputs:
# .csv file summarising the results of the statistical tests in a table.
# .csv file recording all the DEGs that overlap between the current study and 
# Chen et al. (2014) or Pacifico et al. (2018).
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# org.Dm.eg.db loaded via the script in the "Load data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("40_compare_dmel_genage_prep.R")

# FUNCTION DEFINITIONS ----

CompareDEGs <- function (x, y, z) {
  # Compares the overlap in DEG lists between the current study and Dmel data 
  # from previous studies for a given tissue, and direction of differential 
  # expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("head" or
  #      "fatbody").
  #   y: string denoting the name of treatment to be analysed. 
  #       ("M" or "H").
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in RNA_2 compared to RNA_1) or down-
  #      regulated genes (genes more highly expressed in RNA_1 than RNA_2) 
  #      ("up" or "down").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap. 
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Set variables based on arguments.
  
  if (x == "fatbody" && z == "up") {
    dmelList <- chenUpDEGsFB
  } else if (x == "fatbody" && z == "down") {
    dmelList <- chenDownDEGsFB 
  } else if (x == "head" && z == "up") {
    dmelList <- pacificoUpDEGsFB
  } else if (x == "head" && z == "down") {
    dmelList <- pacificoDownDEGsFB
  } else {
    stop ("Argument x or z is incorrect.")
  }
  
  # Load current study lists.
  
  currentStudyDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                      z , "_regulated_LFC0.csv"))
  
  backgroundList <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Reduce number of columns and rename.
  
  backgroundList <- as.data.frame(backgroundList[, 1])
  
  colnames(backgroundList) <- "symbol"
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Treatment" = y,
                          "Direction_of_expression" = paste0(z, "-regulated"),
                          "Number_current_study_DEGs" = length(currentStudyDEGs[, 1]),
                          "Number_previous_Dmel_study_DEGs" = length(dmelList[, 1]),
                          "Number_overlapping_DEGs" = NA,
                          "Number_DEGs_only_in_current_study" = NA,
                          "Number_DEGs_only_in_previous_Dmel_study" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = (0.05/4),
                          "Percentage_of_current_study_DEGs_overlapping" = NA) 
  
  # Calculate overlap between current study and previous Dmel study.
  
  # Determine overlapping genes.
  
  overlappingGenes <- currentStudyDEGs[currentStudyDEGs$symbol %in% dmelList$Dmel_gene_ID, ]
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$symbol)
  
  resultsDF$Percentage_of_current_study_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_current_study_DEGs)*100, digits = 1)
  
  # Prepare overlapping gene IDs for output.
  # Write results to data.frame, with extra blank space so that
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes[, c("symbol", "name")]
  
  colnames(genesToOutput) <- c(paste(x, y, z, "regulated_Flybase_gene_ID", sep = "_"),
                               paste(x, y, z, "regulated_Dmel_gene_name", sep = "_"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:2000, ] <- ""
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_DEGs_only_in_current_study <- 
    resultsDF$Number_current_study_DEGs - resultsDF$Number_overlapping_DEGs
  
  resultsDF$Number_DEGs_only_in_previous_Dmel_study <-
    resultsDF$Number_previous_Dmel_study_DEGs - resultsDF$Number_overlapping_DEGs
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundList$symbol) - resultsDF$Number_overlapping_DEGs - 
    resultsDF$Number_DEGs_only_in_current_study - resultsDF$Number_DEGs_only_in_previous_Dmel_study
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both DEG lists, only 
  # M DEG list, only H DEG list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_DEGs, 
                               resultsDF$Number_DEGs_only_in_current_study,
                               resultsDF$Number_DEGs_only_in_previous_Dmel_study, 
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

setwd("../02_outputs/02_head/31_DESeq2_DEG_lists/")

# Compare DEG lists.

headMUpOverlap <- CompareDEGs("head", "M", "up")

headMDownOverlap <- CompareDEGs("head", "M", "down")

headHUpOverlap <- CompareDEGs("head", "H", "up")

headHDownOverlap <- CompareDEGs("head", "H", "down")

# Fat Body DEG Overlap ----

# Change directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare DEG lists.

fatBodyMUpOverlap <- CompareDEGs("fatbody", "M", "up")

fatBodyMDownOverlap <- CompareDEGs("fatbody", "M", "down")

fatBodyHUpOverlap <- CompareDEGs("fatbody", "H", "up")

fatBodyHDownOverlap <- CompareDEGs("fatbody", "H", "down")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(headMUpOverlap$stats,
                       headMDownOverlap$stats,
                       headHUpOverlap$stats,
                       headHDownOverlap$stats,
                       fatBodyMUpOverlap$stats,
                       fatBodyMDownOverlap$stats,
                       fatBodyHUpOverlap$stats,
                       fatBodyHDownOverlap$stats)

# Combine genes results.

combinedGenes <- cbind(headMUpOverlap$genes,
                       headMDownOverlap$genes,
                       headHUpOverlap$genes,
                       headHDownOverlap$genes,
                       fatBodyMUpOverlap$genes,
                       fatBodyMDownOverlap$genes,
                       fatBodyHUpOverlap$genes,
                       fatBodyHDownOverlap$genes)

# Save Results ----

setwd("../../00_all_tissues/")

write.csv(combinedStats, "23_NER0008751_obj3_exp3_table_S8_Dmel_studies_stats_results.csv", 
          row.names = FALSE)

write.csv(combinedGenes, "24_NER0008751_obj3_exp3_table_D6_Dmel_studies_overlapping_genes.csv", 
          row.names = FALSE)
