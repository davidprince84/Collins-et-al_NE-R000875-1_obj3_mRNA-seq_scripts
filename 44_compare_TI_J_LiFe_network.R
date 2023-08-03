#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: Comparative analysis with other data.
# Tasks: Compare current study differentially expressed genes (DEGs) with 
# Drosophila melanogaster (Dmel) TI-J-LiFe network genes, as defined by Korb 
# et al. (2021).
#-------------------------------------------------------------------------------
# Inputs:
# Dmel DEGs from the current study. List of Dmel TI-J-LiFe genes from Korb et 
# al. (2021).

# Outputs:
# .csv file recording all the results of the statistical tests.
# .csv file recording all the genes that overlap between the TI-J-LiFe genes and 
# the DEGs from the current study.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Load Korb et al. (2021) Data ----

# Set working directory.

setwd("../00_data/05_gene_list_csv")

# Load data.

TIJLiFeGenes <- read.csv("korb_2021_table_s1.csv",
                         stringsAsFactors = FALSE)

# FUNCTION DEFINITIONS ----

CompareDEGsAndDmelGenes <- function (x, y, z = "all") {
  # Compares the overlap between Dmel DEG lists from the current study and Dmel 
  # genes in the TI-J-LiFe list.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("head", 
  #      "fatbody" or "ovaries").
  #   y: string denoting the name of treatment to be analysed. 
  #       ("M" or "H").
  #   z: number denoting the number of DEGs from the current study to compare 
  #      with the TI-J-LiFe genes. Number will select the top z most 
  #      positive and negative expressed DEGs based on log fold change, 
  #      therefore overall number of DEGs is 2 x z. (Default = all DEGs).
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between the Dmel genes 
  #      and the DEGs from the current study. 
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
  
  colnames(backgroundList) <- "Dmel_gene_ID"
  
  # Combine up- and down-regulated DEGs.
  
  currentStudyDEGs <- rbind(upDEGs, downDEGs)  
  
  # Filter DEGs, as specified by argument z.
  
  if (z == "all") {
    filteredDEGs <- currentStudyDEGs  
  } else if ((2 * z) > length(currentStudyDEGs[, 1])) {
    print('Argument z is too large, there are not suffient DEGs, using all DEGs instead.')
    filteredDEGs <- currentStudyDEGs
    z <- "all, as z too large"
  } else {
    # Make all log2FoldChange values absolute, so that they can be ordered by
    # change irrespective of direction.
    
    currentStudyDEGs$log2FoldChange <- abs(currentStudyDEGs$log2FoldChange)
    
    # Sort data.frame so that rows are sorted by logFoldChange.
    
    sortedDEGs <- 
      currentStudyDEGs[order(-(currentStudyDEGs$log2FoldChange)), ]
    
    # Index the top 2 x z genes.
    
    filteredDEGs <- sortedDEGs[(1:(2 * z)), ]
    
  }
  
  # Simplify to just the gene symbols and rename column.
  
  filteredDEGs <- as.data.frame(filteredDEGs[, c("symbol", "name")])
  
  colnames(filteredDEGs) <- c("Flybase_gene_ID", "name")
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Treatment" = y,
                          "Number_top_positive_and_negative_DEGs" = z,
                          "Comparison_list" = "TI-J-LiFe",
                          "Number_DEGs" = length(filteredDEGs[, 1]),
                          "Number_TI_J_LiFe_genes" = length(TIJLiFeGenes[, 1]),
                          "Number_overlapping_genes" = NA,
                          "Number_genes_only_in_DEGs" = NA,
                          "Number_genes_only_in_TI_J_LiFe" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = NA,
                          "Percentage_of_TI_J_LiFe_overlapping" = NA) 
  
  # Calculate overlap between DEGs and TI-J-LiFe.
  
  # Determine overlapping genes and rename columns.
  
  overlappingGenes <- 
    filteredDEGs[filteredDEGs$Flybase_gene_ID %in% TIJLiFeGenes$Flybase_gene_ID, ]
  
  valuesOfArgs <- paste(x, y, z, sep = "_")
  
  colnames(overlappingGenes) <- c(paste0(valuesOfArgs, "_Flybase_gene_ID"),
                                  paste0(valuesOfArgs, "_gene_name"))
                                 
  # Write results to data.frame, with extra blank space so that 
  # all results can be combined later.
  
  overlappingGenesDF <- overlappingGenes
  
  if (is.na(overlappingGenesDF[1, 1])) {
    overlappingGenesDF[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(overlappingGenesDF[, 1]) + 1
  
  overlappingGenesDF[fillStart:(length(TIJLiFeGenes[, 1])), ] <- ""
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_genes <- length(overlappingGenes[, 1])
  
  resultsDF$Percentage_of_TI_J_LiFe_overlapping <- 
    round((resultsDF$Number_overlapping_genes/resultsDF$Number_TI_J_LiFe_genes)*100, digits = 1)
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_genes_only_in_DEGs <- 
    resultsDF$Number_DEGs - resultsDF$Number_overlapping_genes
  
  resultsDF$Number_genes_only_in_TI_J_LiFe <-
    resultsDF$Number_TI_J_LiFe_genes - resultsDF$Number_overlapping_genes
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundList[, 1]) - resultsDF$Number_overlapping_genes - 
    resultsDF$Number_genes_only_in_DEGs - resultsDF$Number_genes_only_in_TI_J_LiFe
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of overlapping genes, only 
  # present in DEG list, only present in TIJLife list, and present in neither 
  # list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_genes, 
                               resultsDF$Number_genes_only_in_DEGs,
                               resultsDF$Number_genes_only_in_TI_J_LiFe, 
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
  
  listToReturn <- list("genes" = overlappingGenesDF, 
                       "stats" = resultsDF)
  
  return(listToReturn)
  
}

# EXECUTED STATEMENTS ----

# Head Comparison ----

# Set working directory.

setwd("../../02_outputs/02_head/31_DESeq2_DEG_lists/")

# Compare gene lists.

headMAllResults <- CompareDEGsAndDmelGenes("head", "M")

headMTop50Results <- CompareDEGsAndDmelGenes("head", "M", 50)

headMTop100Results <- CompareDEGsAndDmelGenes("head", "M", 100)

headMTop200Results <- CompareDEGsAndDmelGenes("head", "M", 200)

headMTop300Results <- CompareDEGsAndDmelGenes("head", "M", 300)

headMTop500Results <- CompareDEGsAndDmelGenes("head", "M", 500)

headHAllResults <- CompareDEGsAndDmelGenes("head", "H")

headHTop50Results <- CompareDEGsAndDmelGenes("head", "H", 50)

headHTop100Results <- CompareDEGsAndDmelGenes("head", "H", 100)

headHTop200Results <- CompareDEGsAndDmelGenes("head", "H", 200)

headHTop300Results <- CompareDEGsAndDmelGenes("head", "H", 300)

headHTop500Results <- CompareDEGsAndDmelGenes("head", "H", 500)

# Fat Body Comparison ----

# Set working directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare gene lists.

fatBodyMAllResults <- CompareDEGsAndDmelGenes("fatbody", "M")

fatBodyMTop50Results <- CompareDEGsAndDmelGenes("fatbody", "M", 50)

fatBodyMTop100Results <- CompareDEGsAndDmelGenes("fatbody", "M", 100)

fatBodyMTop200Results <- CompareDEGsAndDmelGenes("fatbody", "M", 200)

fatBodyMTop300Results <- CompareDEGsAndDmelGenes("fatbody", "M", 300)

fatBodyMTop500Results <- CompareDEGsAndDmelGenes("fatbody", "M", 500)

fatBodyHAllResults <- CompareDEGsAndDmelGenes("fatbody", "H")

fatBodyHTop50Results <- CompareDEGsAndDmelGenes("fatbody", "H", 50)

fatBodyHTop100Results <- CompareDEGsAndDmelGenes("fatbody", "H", 100)

fatBodyHTop200Results <- CompareDEGsAndDmelGenes("fatbody", "H", 200)

fatBodyHTop300Results <- CompareDEGsAndDmelGenes("fatbody", "H", 300)

fatBodyHTop500Results <- CompareDEGsAndDmelGenes("fatbody", "H", 500)

# Ovaries Comparison ----

# Set working directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists/")

# Compare gene lists.

ovariesMAllResults <- CompareDEGsAndDmelGenes("ovaries", "M")

ovariesMTop50Results <- CompareDEGsAndDmelGenes("ovaries", "M", 50)

ovariesMTop100Results <- CompareDEGsAndDmelGenes("ovaries", "M", 100)

ovariesMTop200Results <- CompareDEGsAndDmelGenes("ovaries", "M", 200)

ovariesMTop300Results <- CompareDEGsAndDmelGenes("ovaries", "M", 300)

ovariesMTop500Results <- CompareDEGsAndDmelGenes("ovaries", "M", 500)

ovariesHAllResults <- CompareDEGsAndDmelGenes("ovaries", "H")

ovariesHTop50Results <- CompareDEGsAndDmelGenes("ovaries", "H", 50)

ovariesHTop100Results <- CompareDEGsAndDmelGenes("ovaries", "H", 100)

ovariesHTop200Results <- CompareDEGsAndDmelGenes("ovaries", "H", 200)

ovariesHTop300Results <- CompareDEGsAndDmelGenes("ovaries", "H", 300)

ovariesHTop500Results <- CompareDEGsAndDmelGenes("ovaries", "H", 500)

# Combine Results and Save ----

# Combine results for stats.

combinedStats <- rbind(headMTop50Results$stats, headMTop100Results$stats,
                       headMTop200Results$stats, headMTop300Results$stats,
                       headMTop500Results$stats, headMAllResults$stats, 
                       headHTop50Results$stats, headHTop100Results$stats,
                       headHTop200Results$stats, headHTop300Results$stats,
                       headHTop500Results$stats, headHAllResults$stats,
                       fatBodyMTop50Results$stats, fatBodyMTop100Results$stats,
                       fatBodyMTop200Results$stats, fatBodyMTop300Results$stats,
                       fatBodyMTop500Results$stats, fatBodyMAllResults$stats,
                       fatBodyHTop50Results$stats, fatBodyHTop100Results$stats, 
                       fatBodyHTop200Results$stats, fatBodyHTop300Results$stats,
                       fatBodyHTop500Results$stats, fatBodyHAllResults$stats,
                       ovariesMTop50Results$stats, ovariesMTop100Results$stats,
                       ovariesMTop200Results$stats, ovariesMTop300Results$stats,
                       ovariesMTop500Results$stats, ovariesMAllResults$stats,
                       ovariesHTop50Results$stats, ovariesHTop100Results$stats,
                       ovariesHTop200Results$stats, ovariesHTop300Results$stats,
                       ovariesHTop500Results$stats, ovariesHAllResults$stats)

# Combine results for overlapping genes.

combinedGenes <- cbind(headMTop50Results$genes, headMTop100Results$genes,
                       headMTop200Results$genes, headMTop300Results$genes,
                       headMTop500Results$genes, headMAllResults$genes, 
                       headHTop50Results$genes, headHTop100Results$genes,
                       headHTop200Results$genes, headHTop300Results$genes,
                       headHTop500Results$genes, headHAllResults$genes,
                       fatBodyMTop50Results$genes, fatBodyMTop100Results$genes,
                       fatBodyMTop200Results$genes, fatBodyMTop300Results$genes,
                       fatBodyMTop500Results$genes, fatBodyMAllResults$genes,
                       fatBodyHTop50Results$genes, fatBodyHTop100Results$genes, 
                       fatBodyHTop200Results$genes, fatBodyHTop300Results$genes,
                       fatBodyHTop500Results$genes, fatBodyHAllResults$genes,
                       ovariesMTop50Results$genes, ovariesMTop100Results$genes,
                       ovariesMTop200Results$genes, ovariesMTop300Results$genes,
                       ovariesMTop500Results$genes, ovariesMAllResults$genes,
                       ovariesHTop50Results$genes, ovariesHTop100Results$genes,
                       ovariesHTop200Results$genes, ovariesHTop300Results$genes,
                       ovariesHTop500Results$genes, ovariesHAllResults$genes)

# Save results.

setwd("../../00_all_tissues/")

write.csv(combinedStats, "30_NER0008751_obj3_exp3_table_S10_TI-J-LiFe_stats_results.csv",
          row.names = FALSE)

write.csv(combinedGenes, "31_NER0008751_obj3_exp3_table_D8_TI-J-LiFe_overlapping_genes.csv",
          row.names = FALSE)
