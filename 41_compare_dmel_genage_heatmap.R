#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: Comparative analysis with other data.
# Tasks: Comparison of differentially expressed genes (DEGs) from the current 
# study with Drosophila melanogaster (Dmel) GenAge data set.
#-------------------------------------------------------------------------------
# Inputs:
# GenAge Dmel gene list (from script #40). Lists of current study DEGs.

# Outputs:
# A heatmap of the relevant results.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# org.Dm.eg.db loaded via the script in the "Load data" section.
# pheatmap and ggplot2 loaded via the script in the "Custom functions" section.

# LOAD DATA ----

# Load GenAge Data ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("40_compare_dmel_genage_prep.R")

# Add a new column for current study expression to genAgeGenes.

genAgeGenes$Study_expression <- NA

# LOAD CUSTOM FUNCTIONS ----

source("f_PlotCustomHeatmap.R")

# FUNCTION DEFINITIONS ----

ReportPresenceOfGenAgeGenes <- function (x, y = "none", z = "none") {
  # Reports whether a list of DEGs contains Dmel GenAge genes.
  #
  # Args:
  #   x: string denoting the study that the gene list is taken from ("current", 
  #      "chen" or "pacifico").
  #   y: string denoting which tissue is to be analysed for Dmel genes from the
  #      current study ("none" (default), "fatbody", "head" or "ovaries").
  #   z: string denoting the name of treatment from the current study to be analysed. 
  #      ("none" (default), "M" or "H").
  #
  # Returns:
  #   genAgeGenes data.frame with presence of differentially expressed genes 
  #   for stated study added.
  #
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue of the current study.
  
  # Load DEGs list based on arguments.
  
  if (x == "current") {
    upregList <- read.csv(paste0(y, "_results_", z, 
                                 "_treatment_up_regulated_LFC0.csv"))
    downregList <- read.csv(paste0(y, "_results_", z, 
                                   "_treatment_down_regulated_LFC0.csv"))  
    columnName <- "symbol"
  } else if (x == "chen") {
    upregList <- chenUpDEGsFB
    downregList <- chenDownDEGsFB
    columnName <- "Dmel_gene_ID"
  } else if (x == "pacifico") {
    upregList <- pacificoUpDEGsFB
    downregList <- pacificoDownDEGsFB
    columnName <- "Dmel_gene_ID"
  } else {
    stop('Argument x must be "current", "chen" or "pacifico".')
  }
  
  # Designate internal copy of genAgeGenes.
  
  genAge <- genAgeGenes
  
  # Loop over GenAge data and add "1" to Study_expression column.
  
  for(i in (1:length(genAge[, 1]))) {
    if (genAge[i, "Flybase_gene_ID"] %in% upregList[, columnName]) {
      genAge[i, "Study_expression"] <- 1
    } else if (genAge[i, "Flybase_gene_ID"] %in% downregList[, columnName]) {
      genAge[i, "Study_expression"] <- -1
    } else {
      genAge[i, "Study_expression"] <- 0
    }
  }
  
  # Return genAge
  
  return(genAge)
  
}

# EXECUTED STATEMENTS ----

# Current study.

# Ovaries.

# Change working directory.

setwd("../02_outputs/01_ovaries/31_DESeq2_DEG_lists/")

# Run function.

ovariesMResults <- ReportPresenceOfGenAgeGenes("current", "ovaries", "M")

ovariesHResults <- ReportPresenceOfGenAgeGenes("current", "ovaries", "H")

# Fat body.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

fatBodyMResults <- ReportPresenceOfGenAgeGenes("current", "fatbody", "M")

fatBodyHResults <- ReportPresenceOfGenAgeGenes("current", "fatbody", "H")

# Head.

setwd("../../02_head/31_DESeq2_DEG_lists/")

headMResults <- ReportPresenceOfGenAgeGenes("current", "head", "M")

headHResults <- ReportPresenceOfGenAgeGenes("current", "head", "H")

# Dmel previous studies.

# Chen et al. (2014)

chenResults <- ReportPresenceOfGenAgeGenes("chen")

# Pacifico et al. (2018) 

pacificoResults <- ReportPresenceOfGenAgeGenes("pacifico")

# Combine results.

combinedResults <- genAgeGenes[, c(2, 3, 7, 8, 9)]

# Add Current and previous study results.

combinedResults$pacifico_dmel_head <- pacificoResults[, 10]

combinedResults$current_head_M <- headMResults[, 10]

combinedResults$current_head_H <- headHResults[, 10]

combinedResults$chen_dmel_fatbody <- chenResults[, 10]

combinedResults$current_fatbody_M <- fatBodyMResults[, 10]

combinedResults$current_fatbody_H <- fatBodyHResults[, 10]

combinedResults$current_ovaries_M <- ovariesMResults[, 10]

combinedResults$current_ovaries_H <- ovariesHResults[, 10]

# Plot and Save Heatmap ----

# Prepare data for heatmap.

# Heatmap showing genes which are differentially expressed in a Dmel study 
# (Chen or Pacifico) and at least 1 current study sample, or in current study
# head or fat body.

# Remove rows where only DEGs are in ovaries.

noZeroResults <- combinedResults[rowSums(combinedResults[, c(7, 8, 10, 11)] == 0, na.rm=TRUE) < 4, ]

# Remove ‘Juvenile hormone esterase binding protein 29 (DmP29)’ as it is annotated
# with both pro- and anti-longevity influences, and R won't allow
# this to be annotated on a heatmap (as the row names in the annotation
# data.frame would be the same).

noZeroResults <- noZeroResults[!(noZeroResults[, "name"] == "Juvenile hormone esterase binding protein 29 (DmP29)"), ]

# Format names where necessary.

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0033153"), "name"] <- 
  "Growth arrest and DNA damage-inducible 45"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0036813"), "name"] <- 
  "Aut1"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0010504"), "name"] <- 
  "kermit"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0028425"), "name"] <- 
  "Juvenile hormone Inducible-21"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0023215"), "name"] <- 
  "Mnt" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0031037"), "name"] <- 
  "CG14207" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0028717"), "name"] <- 
  "Lnk" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0042132"), "name"] <- 
  "CG18809" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0260990"), "name"] <- 
  "yata"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0020238"), "name"] <- 
  "14-3-3 epsilon" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0039044"), "name"] <- 
  "p53" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0023511"), "name"] <- 
  "ER degradation enhancer, mannosidase alpha-like 1"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0038475"), "name"] <- 
  "Keap1"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0010379"), "name"] <- 
  "Akt kinase"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0015279"), "name"] <- 
  "Phosphatidylinositol 3-kinase 92E"	

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0038788"), "name"] <- 
  "Sirtuin 2"	

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0024291"), "name"] <- 
  "Sirtuin 1"	

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0015765"), "name"] <- 
  "p38a MAP kinase"	

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0026379"), "name"] <- 
  "Phosphatase and tensin homolog"	

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0003028"), "name"] <- 
  "ovo" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0015805"), "name"] <- 
  "Histone deacetylase 1" 

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0261560"), "name"] <- 
  "Thor" 
  
noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0030718"), "name"] <- 
  "NADH dehydrogenase (ubiquinone) 20 kDa subunit" 

# Plotting matrix of +1/0/-1 values

plotMatrix <- as.matrix(noZeroResults[, c(6:13)])

# Add gene names to the matrix.

row.names(plotMatrix) <- noZeroResults[, "name"]

# Plot heatmap.

heatmapPlot <- PlotCustomHeatmap(41)

# Save heatmap.

setwd("../../00_all_tissues")

ggsave("20_NER0008751_obj3_exp3_fig_S13_GenAge_heatmap.svg", plot = heatmapPlot, 
       width = 42, height = 58, units = "cm")
