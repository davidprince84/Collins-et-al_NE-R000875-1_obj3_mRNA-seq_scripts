#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Conduct gene ontology (GO) enrichment analysis on lists of 
# differentially expressed genes (DEGs) from the current study.
#-------------------------------------------------------------------------------
# Inputs:
# Overall gene lists for each tissue from the current study.

# Outputs:
# .txt files stating the significantly enriched Biological Processes GO terms.
# NOTE: .txt files used rather than .csv files to prevent Excel from automatically
# formatting GeneRatio and BgRatio columns, which can result in dates being 
# inserted rather than fractions.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(clusterProfiler)  # enrichGO(), simplify().
library(org.Dm.eg.db)

# LOAD DATA ----

# No data loaded.

# FUNCTION DEFINITIONS ----

CalculateEnrichedGOTerms <- function(x, y, z) {
  # Calculates the enriched GO terms for a list of Drosophila melanogaster (Dmel)
  # DEGs using GO annotations for the Dmel. Shows which GO terms are redundant.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("fatbody", "head" 
  #      or "ovaries").
  #   y: string denoting the name of treatment to be analysed. 
  #       ("M" or "H").
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in RNA_2 compared to RNA_1) or down-
  #      regulated genes (genes more highly expressed in RNA_1 than RNA_2) 
  #      ("up" or "down").
  #
  # Returns:
  #   Nothing. A .txt file is saved to the working directory with the results.
  #
  # NOTE: The working directory needs to be set to the location where the results
  # are to be saved.
  
  # Set variables based on arguments.
  
  if (y == "M") {
    columnNameToIndex1 <- "DEG_in_M"
    columnNameToIndex2 <- "expression_with_age_in_M"
  } else if (y == "H") {
    columnNameToIndex1 <- "DEG_in_H"
    columnNameToIndex2 <- "expression_with_age_in_H"
  } else {
    stop('Argument y must equal "M" or "H".')
  }
  
  if (z == "up") {
    directionOfDEGs <- "up-regulated"
  } else if (z == "down") {
    directionOfDEGs <- "down-regulated"
  } else {
    stop('Argument z must equal "up" or "down".')
  }
  
  # Load data.frame of all expressed genes for tissue x.
  
  setwd("../31_DESeq2_DEG_lists")
  
  allGenes <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  setwd("../40_GO_enrichment_analysis")
  
  # Remove columns no longer needed.
  
  allGenes <- allGenes[, c(1, 20, 21, 28, 29)]

  # Index DEGs based on arguments y and z.
  
  yDEGs <- 
    allGenes[allGenes[, columnNameToIndex1] == "YES" ,]
  
  yzDEGs <- 
    yDEGs[yDEGs[, columnNameToIndex2] == directionOfDEGs ,]
  
  # Conduct GO analysis.
  
  goEnrichResults <- enrichGO(yzDEGs$gene_symbol,
                              "org.Dm.eg.db",
                              keyType = "FLYBASE", 
                              ont = "BP",
                              universe = allGenes$gene_symbol,
                              qvalueCutoff = 0.05, minGSSize = 5)
  
  # Make data.frame of the enrichment results (if any).
  
  if (!is.null(goEnrichResults)) {
    goEnrichDF <- as.data.frame(goEnrichResults, stringsAsFactors = FALSE)
  }
  
  # Set file name.
  
  fileName <- paste(x, y, z, "regulated_DEGs_enriched_GO_results.txt", sep = "_") 
  
  # Save file with "No enriched GO terms found." if necessary, else
  # simplify GO results to remove redundancy.
  
  if (is.null(goEnrichResults)) {
    # Write a file stating that no enriched GO terms were returned.
    finalTable <- data.frame("Description" = "No enriched GO terms found.")
  } else if (length(goEnrichDF$ID) == 0) {
    # Write a file stating that no enriched GO terms were returned.
    finalTable <- data.frame("Description" = "No enriched GO terms found.")
  } else {
    # Simplify GO results.
    
    simpleGOResults <- as.data.frame(simplify(goEnrichResults),
                                     stringsAsFactors = FALSE)
    
    # Add column to specify terms are not redundant.
    
    simpleGOResults$redundant <- rep(0, times = length(simpleGOResults$ID))
    
    # Convert enrichResult object to data.frame.
    
    goEnrichDF <- as.data.frame(goEnrichResults, stringsAsFactors = FALSE)
    
    # Index the redundant terms from the results.
    
    redundantGODF <- goEnrichDF[!(goEnrichDF$ID %in% simpleGOResults$ID), ]
    
    # Add column specifying term is redundant.
    
    redundantGODF$redundant <- rep(1, times = length(redundantGODF$ID))
    
    # Bind the rows of the two data frames together.
    
    finalTable <- rbind(simpleGOResults, redundantGODF)
  } 
  
  # Save the results.
  
  write.table(finalTable, fileName, sep = "\t", row.names = FALSE)

} 

# EXECUTED STATEMENTS ----

# Calculate Enriched GO terms ----

# Head.

# Set working directory.
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

setwd("../02_outputs/02_head/")

# Create directory.

dir.create("40_GO_enrichment_analysis")
# Will produce a warning if directory already exists.

# Change directory.

setwd("40_GO_enrichment_analysis")

# Calculate enriched GO terms.

CalculateEnrichedGOTerms("head", "M", "up")

CalculateEnrichedGOTerms("head", "M", "down")

CalculateEnrichedGOTerms("head", "H", "up")

CalculateEnrichedGOTerms("head", "H", "down")

# Fat body.

# Create directory.

dir.create("../../03_fatbody/40_GO_enrichment_analysis")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../03_fatbody/40_GO_enrichment_analysis")

# Calculate enriched GO terms.

CalculateEnrichedGOTerms("fatbody", "M", "up")

CalculateEnrichedGOTerms("fatbody", "M", "down")

CalculateEnrichedGOTerms("fatbody", "H", "up")

CalculateEnrichedGOTerms("fatbody", "H", "down")

# Ovaries.

# Create directory.

dir.create("../../01_ovaries/40_GO_enrichment_analysis")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../01_ovaries/40_GO_enrichment_analysis")

# Calculate enriched GO terms.

CalculateEnrichedGOTerms("ovaries", "M", "up")

CalculateEnrichedGOTerms("ovaries", "M", "down")

CalculateEnrichedGOTerms("ovaries", "H", "up")

CalculateEnrichedGOTerms("ovaries", "H", "down")
