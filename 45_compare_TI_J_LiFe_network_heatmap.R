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
# Heatmap of TI-J-LiFe genes differentially expressed in the current study DEGs.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# pheatmap and ggplot2 loaded via the script in the "Custom functions" section.

# LOAD CUSTOM FUNCTIONS ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("f_PlotCustomHeatmap.R")

# LOAD DATA ----

# Load Korb et al. (2021) Data ----

# Set working directory.

setwd("../00_data/05_gene_list_csv")

# Load data.

TIJLiFeGenes <- read.csv("korb_2021_table_s1.csv",
                         stringsAsFactors = FALSE)

# FUNCTION DEFINITIONS ----

ReportPresenceOfTIJLiFeGenes <- function (x, y) {
  # Reports whether a list of DEGs contains Dmel TIJLiFe genes.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("head",
  #      "fatbody" or "ovaries").
  #   y: string denoting the name of treatment to be analysed. 
  #       ("M" or "H").
  #
  # Returns:
  #   Data.frame of results.
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Load DEGs list based on arguments.
  
  upregList <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                "up_regulated_LFC0.csv"))
  
  downregList <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                  "down_regulated_LFC0.csv"))
  
  # Loop over TIJLiFeGenes and add "1" to Study_expression column.
  
  for(i in (1:length(TIJLiFeGenes[, 1]))) {
    if (TIJLiFeGenes[i, "Flybase_gene_ID"] %in% upregList[, "symbol"]) {
      TIJLiFeGenes[i, "Study_expression"] <- 1
    } else if (TIJLiFeGenes[i, "Flybase_gene_ID"] %in% downregList[, "symbol"]) {
      TIJLiFeGenes[i, "Study_expression"] <- -1
    } else {
      TIJLiFeGenes[i, "Study_expression"] <- 0
    }
  }
  
  # Return TIJLiFeGenes.
  
  return(TIJLiFeGenes)
  
}

# EXECUTED STATEMENTS ----

# Current study.

# Ovaries.

# Change working directory to 02_outputs.

setwd("../../02_outputs/01_ovaries/31_DESeq2_DEG_lists/")

# Run function.

ovariesMResults <- ReportPresenceOfTIJLiFeGenes("ovaries", "M")

ovariesHResults <- ReportPresenceOfTIJLiFeGenes("ovaries", "H")

# Fat body.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

fatBodyMResults <- ReportPresenceOfTIJLiFeGenes("fatbody", "M")

fatBodyHResults <- ReportPresenceOfTIJLiFeGenes("fatbody", "H")

# Head.

setwd("../../02_head/31_DESeq2_DEG_lists/")

headMResults <- ReportPresenceOfTIJLiFeGenes("head", "M")

headHResults <- ReportPresenceOfTIJLiFeGenes("head", "H")

# Combine results.

combinedResults <- headMResults

colnames(combinedResults) <- 
  c("Flybase_gene_ID", "name", "head_M")

# Add current study results.

combinedResults$head_H <- headHResults[, 3]

combinedResults$fatbody_M <- fatBodyMResults[, 3]

combinedResults$fatbody_H <- fatBodyHResults[, 3]

combinedResults$ovaries_M <- ovariesMResults[, 3]

combinedResults$ovaries_H <- ovariesHResults[, 3]

# Plot and Save Heatmap ----

# Prepare data for heatmap.

# Heatmap showing genes which are differentially expressed in the current study.
# i.e. remove rows where all the DEG lists show 0 for expression.

noZeroResults <- combinedResults[rowSums(combinedResults[, 3:8] == 0, na.rm = TRUE) < 6, ]

# Format names where necessary.

noZeroResults$name

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0010548"), "name"] <- 
  "Aldehyde dehydrogenase type III FBgn0010548"

noZeroResults[which(noZeroResults$Flybase_gene_ID == "FBgn0012036"), "name"] <- 
  "Aldehyde dehydrogenase type III FBgn0012036"

# Plotting matrix of +1/0/-1 values.

plotMatrix <- as.matrix (noZeroResults[, c(3:8)])

# Add gene names to the matrix.

row.names(plotMatrix) <- noZeroResults[, "name"]

# Plot heatmap.

heatmapPlot <- PlotCustomHeatmap(45)

# Save heatmap.

setwd("../../00_all_tissues")

ggsave("32_NER0008751_obj3_exp3_fig_5_TIJLiFe_heatmap.svg", plot = heatmapPlot, 
       width = 42, height = 58, units = "cm")
