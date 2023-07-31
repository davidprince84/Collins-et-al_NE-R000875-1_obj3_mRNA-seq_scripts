#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce a heatmap of differentially expressed genes for each
# tissue.
#-------------------------------------------------------------------------------
# Inputs:
# Dmel_6_22_97_transcripts2genes.txt file, dmelAnnotations object and Kallisto 
# abundances.

# Outputs:
# .svg figure of a heatmap for each tissue.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ggplot2)  # ggsave().
library(pheatmap)  # pheatmap().
library(RColorBrewer)  # source of annotation colours.
# org.Dm.eg.db loaded via script in "Load data" section. DESeq2 and tximport 
# loaded via the script in the "Custom functions" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Load Transcript to Gene File ----

setwd("../00_data/10_transcripts2genes")

t2g <- read.table(file = "Dmel_6_22_97_transcripts2genes.txt",
                  header = FALSE,
                  col.names = c("TXNAME",
                                "GENEID"))

# Load DmelAnnotations Data.Frame ----

setwd("../../01_scripts")

source("21_generate_dmel_gene_names.R")

colnames(dmelAnnotations) <- c("symbol", "gene_name")

# LOAD CUSTOM FUNCTIONS ----

source("f_BuildDDSForDEAnalysis.R")

# FUNCTION DEFINITIONS ----

DetermineTop50DEGs <- function (x) {
  # Generate data.frame of the 50 most highly differentially expressed genes (DEGs) 
  # (or all DEGs if fewer than 50 in total) for a given tissue.
  # 
  # Args:
  #   x: string denoting which tissue is to be analysed ("fatbody", "head"  
  #      or "ovaries").
  #
  # Returns:
  #   A data.frame of the most highly differentially expressed genes.
  
  # Set variables based on argument.
  
  if (x == "head") {
    dds <- ddsHead
    reorderColumns <- c(1, 5, 8, 2, 6, 9, 3, 7, 10, 4, 11, 12)
  } else if (x == "fatbody") {
    dds <- ddsFatBody
    reorderColumns <- c(1, 5, 9, 2, 6, 10, 3, 7, 12, 4, 8, 11)
  } else if (x == "ovaries") {
    dds <- ddsOvaries
    reorderColumns <- c(1, 5, 9, 2, 6, 10, 3, 7, 12, 4, 8, 11)
  } else {
    stop('Argument x must be "fatbody", "head" or "ovaries".')
  }
  
  # Generate data.frame of log transformed counts for each gene.
  
  # rlog transform the counts.
  
  ddsTransform <- rlog(dds)
  
  # Extract transformed counts into a data.frame.
  
  rlogNormCounts <- as.data.frame(assay(ddsTransform))
  
  # Reorder columns in rlogNormCounts.
  
  rlogNormCounts <-rlogNormCounts[, reorderColumns]
  
  # Add a column with the gene symbol to the data.frame.
  
  rlogNormCounts$symbol <- row.names(rlogNormCounts)
  
  # Extract the results of the DEG analysis for each gene, for both the M
  # and H treatments.
  
  mResults <- as.data.frame(results(dds, contrast = c("condition", "M_RNA2", "M_RNA1"),
                            alpha = 0.05, lfcThreshold = 0)) # Stats for M
  
  hResults <- as.data.frame(results(dds, contrast = c("condition", "H_RNA2", "H_RNA1"),
                                          alpha = 0.05, lfcThreshold = 0)) # Stats for H
  
  # Remove NAs from padj column (represent genes excluded from analysis as all 
  # counts were 0 or it contained an extreme count outlier).
  
  mResults <- mResults[!is.na(mResults$padj), ]
  
  hResults <- hResults[!is.na(hResults$padj), ]
  
  # Add gene symbol as a column.
  
  mResults$symbol <- row.names(mResults)
  
  hResults$symbol <- row.names(hResults)
  
  # Index only the significant results (padj < 0.05).
  
  sigMResults <- mResults[mResults[, "padj"] < 0.05 , ]
  
  sigHResults <- hResults[hResults[, "padj"] < 0.05 , ]
  
  # Combine M and H results.
  
  allSigResults <- rbind(sigMResults, sigHResults)
  
  # Make all log2FoldChange values absolute, so that they can be ordered by
  # change irrespective of direction.
  
  allSigResults$log2FoldChange <- abs(allSigResults$log2FoldChange)
  
  # Sort data.frame so that rows are sorted by logFoldChange.
  
  sortedSigResults <- 
    allSigResults[order(-(allSigResults$log2FoldChange)), ]
  
  # Duplicated genes removed for data.frame, leaving the highest logFoldchange.
  
  nonDuplicateSigResults <- sortedSigResults[!duplicated(sortedSigResults$symbol), ]
  
  # Select top 50 DEGs by highest logFoldChange (or all DEGs if fewer).
  
  if (length(nonDuplicateSigResults$log2FoldChange) < 50) {
    numberOfGenes <- length(nonDuplicateSigResults$log2FoldChange)
  } else {
    numberOfGenes <- 50
  }
  
  # Index top DEGs.
  
  topDEGs <- nonDuplicateSigResults[1:numberOfGenes, ]
  
  # Combine topDEGs with the matrix of transformed counts.
  
  plottingData <- merge(topDEGs, rlogNormCounts, all.x = TRUE)
  
  # Remove columns no longer needed.
  
  plottingData <- plottingData[, c(1, 8:19)]
  
  # Add gene names to symbols.
  
  plottingData <- merge(plottingData, dmelAnnotations, all.x = TRUE)
  
  # Return plottingData.
  
  return(plottingData)
  
}

PlotHeatmap <- function (x) {
  # Takes a data.frame of expression data and generates a heatmap.
  #
  # Args:
  #   x: string denoting which tissue the data are for ("fatbody", "head" 
  #      or "ovaries").
  #
  # Returns:
  #   A heatmap showing expression of the genes in the data.frame.
  
  # Set variables based on argument.
  
  if (x == "head") {
    plottingData <- headData
  } else if (x == "fatbody") {
    plottingData <- fatBodyData
  } else if (x == "ovaries") {
    plottingData <- ovariesData
  } else {
    stop('Argument x must be "fatbody", "head" or "ovaries".')
  }
  
  # Transform data.frame into a matrix.
  
  plotMatrix <- as.matrix(plottingData[, 2:13])
  
  # Add gene names to the matrix.
  
  row.names(plotMatrix) <- plottingData[, "gene_name"]
  
  # Generate keys and annotations.
  
  colAnnotations <- 
    data.frame("Time point" = rep(c(rep("RNA_1", times = 3), rep("RNA_2", times = 3)), times = 2), 
               Treatment = c(rep("M", times = 6), rep("H", times = 6)),
               check.names = FALSE)
  
  row.names(colAnnotations) <- colnames(plotMatrix)
  
  annotationColours <- 
    list("Time point" = c("RNA_1" = "#984EA3", "RNA_2" = "#FF7F00"),
         Treatment = c("M" = "#0072B2", "H" = "#56B4E9"))
  
  # Plot heatmap.
  
  heatmapPlot <- pheatmap(plotMatrix,
                          scale = "row",
                          cluster_rows = TRUE,
                          cluster_cols = FALSE,
                          gaps_col = c(3, 6, 9),
                          annotation_col = colAnnotations,
                          annotation_colors = annotationColours,
                          legend = TRUE,
                          cellwidth=15,
                          cellheight=10,
                          fontsize = 7)
  
  # Return heatmap.
  
  return(heatmapPlot)
  
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

# Head.

ddsHead <- BuildDDSForDEAnalysis("head")

# Fat body.

ddsFatBody <- BuildDDSForDEAnalysis("fatbody")

# Ovaries.

ddsOvaries <- BuildDDSForDEAnalysis("ovaries")

# Make Heatmaps ----

# Data for heatmaps.

# Head.

headData <- DetermineTop50DEGs("head")

# Fat body.

fatBodyData <- DetermineTop50DEGs("fatbody")

# Ovaries.

ovariesData <- DetermineTop50DEGs("ovaries")

# Format gene names.
# Print gene names to check for duplicate names.

# Head.

headData$gene_name

# Add gene name to gene where it is missing.

headData[9, "gene_name"] <- "CG4500"

# Many genes are "uncharacterized proteins", therefore add
# gene symbol to these genes using a loop.

for (ROW in 1:length(headData$symbol)) {
  if (headData[ROW, "gene_name"] == "uncharacterized protein") {
    headData[ROW, "gene_name"] <- 
      paste0(headData[ROW, "gene_name"], " ", headData[ROW, "symbol"])
  }
}

# Fat body.

# Print gene names to check for duplicate names.

fatBodyData$gene_name

# Many genes are "uncharacterized proteins", therefore add
# gene symbol to these genes using a loop.

for (ROW in 1:length(fatBodyData$symbol)) {
  if (fatBodyData[ROW, "gene_name"] == "uncharacterized protein") {
    fatBodyData[ROW, "gene_name"] <- 
      paste0(fatBodyData[ROW, "gene_name"], " ", fatBodyData[ROW, "symbol"])
  }
}

# Ovaries.

# Print gene names to check for duplicate names.

ovariesData$gene_name

# Add gene symbols to genes with very similar names.

# Many genes are "uncharacterized proteins", therefore add
# gene symbol to these genes using a loop.

for (ROW in 1:length(ovariesData$symbol)) {
  if (ovariesData[ROW, "gene_name"] == "uncharacterized protein") {
    ovariesData[ROW, "gene_name"] <- 
      paste0(ovariesData[ROW, "gene_name"], " ", ovariesData[ROW, "symbol"])
  }
}

# Plot heatmaps.

# Head.

headHeatmapPlot <- PlotHeatmap("head")

# Fat body.

fatBodyHeatmapPlot <- PlotHeatmap("fatbody")

# Ovaries.

ovariesHeatmapPlot <- PlotHeatmap("ovaries")

# Save Heatmaps ----

# Head.

# Create directory.

dir.create("../02_outputs/02_head/32_DEG_top50_heatmaps")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../02_outputs/02_head/32_DEG_top50_heatmaps")

# Save plot.

ggsave("00_NER0008751_obj3_exp3_fig_S10_head_heatmap.svg", headHeatmapPlot,
       width = 29, height = 21, units = "cm")

# Fat body.

# Create directory.

dir.create("../../03_fatbody/32_DEG_top50_heatmaps")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../03_fatbody/32_DEG_top50_heatmaps")

# Save plot.

ggsave("00_NER0008751_obj3_exp3_fig_S11_fat_body_heatmap.svg", fatBodyHeatmapPlot,
       width = 29, height = 21, units = "cm")

# Ovaries.

# Create directory.

dir.create("../../01_ovaries/32_DEG_top50_heatmaps")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../01_ovaries/32_DEG_top50_heatmaps")

# Save plot.

ggsave("00_NER0008751_obj3_exp3_fig_S12_ovaries_heatmap.svg", ovariesHeatmapPlot,
       width = 29, height = 21, units = "cm")
