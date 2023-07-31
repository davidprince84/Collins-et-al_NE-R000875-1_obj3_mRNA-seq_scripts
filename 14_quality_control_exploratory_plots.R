#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Produce principal component analysis from Kallisto data.
#-------------------------------------------------------------------------------
# Inputs:
# Dmel_6_22_97_transcripts2genes.txt file and Kallisto abundances.

# Outputs:
# .svg figure of the plots.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ggpubr)  # ggarrange(), annotate_figure()
library(ggrepel)  # geom_text_repel() [also loads ggplot2]
library(reshape2)  # melt()
# DESeq2 and tximport loaded via the script in the "Custom functions" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Transcript to Gene File ----

setwd("../00_data/10_transcripts2genes")

t2g <- read.table(file = "Dmel_6_22_97_transcripts2genes.txt",
                  header = FALSE,
                  col.names = c("TXNAME",
                                "GENEID"))

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_BuildDDSForDEAnalysis.R")

# FUNCTION DEFINITIONS ----

GenerateExploratoryPlots <- function(x, intgroup = "condition", ntop = 2000) {
  # Acknowledgement:
  # The PCA part of this function was adapted from the pcaPlot function in the 
  # DESeq2 R package v1.28.1.
  #
  # Generates a figure showing the normalisation of the gene expression data by 
  # DESeq2 as a boxplot and the clustering of this data as a principal component
  # analysis (PCA) plot.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("fatbody", "head"  
  #      or "ovaries").
  #   intgroup: interesting groups: a character vector of names in colData(x) to 
  #             use for grouping (default of "condition").
  #   ntop: number of top genes to use for principal components, selected by 
  #         highest row variance (default of 2000).
  #
  # Returns:
  #   A ggplot annotated figure with 2 subplots arranged by the ggpubr package. 
  
  # Set variables based on argument.
  
  if (x == "head") {
    dds <- ddsHead
    tissue <- "Head"
  } else if (x == "fatbody") {
    dds <- ddsFatBody
    tissue <- "Fat body"
  } else if (x == "ovaries") {
    dds <- ddsOvaries
    tissue <- "Ovaries"
  } else {
    stop('Argument x must be fatbody", "head" or "ovaries".')
  }
  
  # Transform the data to stabilize the variance across the mean.
  # Select  rlog transformation chosen if there are relatively few 
  # samples (n < 30), or vst transformation if there are more than 30
  # samples.
  
  if (length(dds$sample) < 30) {
    ddsTransform <- rlog(dds)
    yLabelContents <- "rlog_transformed_counts"
  } else if (length(dds$sample) >= 30) {
    ddsTransform <- vst(dds)
    yLabelContents <- "vst_transformed_counts"
  }
  
  # Box plot of normalized, transformed counts.
  
  # Format data for plot.
  
  boxPlotData <- melt(assay(ddsTransform))
  
  colnames(boxPlotData) <- c("gene", "sample", "transformed_counts")
  
  # Calculate median value of counts.
  
  samplesMedian <- median(boxPlotData$transformed_counts)
  
  # Generate box plot.
  
  boxPlot <- ggplot(boxPlotData, aes(x = sample, y = transformed_counts)) + 
    geom_boxplot() +
    geom_hline(yintercept = samplesMedian) +  # Add median line
    ylab(yLabelContents)        
  
  # Format box plot.
  
  formattedBoxPlot <- boxPlot + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5))
  
  # PCA plot.
  
  # Calculate the variance for each gene.
  
  rv <- rowVars(assay(ddsTransform))
  
  # Select the ntop genes by variance.
  
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform a PCA on the data in assay(x) for the selected genes.
  
  pca <- prcomp(t(assay(ddsTransform)[select, ]))
  
  # The contribution to the total variance for each component.
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(ddsTransform)))) {
    stop("The argument 'intgroup' should specify columns of colData(dds).")
  }
  
  intgroup.df <- as.data.frame(colData(ddsTransform)[, intgroup, drop=FALSE])
  
  # Add the intgroup factors together to create a new grouping factor.
  
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(ddsTransform)[[intgroup]]
  }
  
  # Assemble the data for the plot.
  
  d <- data.frame(PC1 = pca$x[,1], 
                  PC2 = pca$x[,2], 
                  group = group, 
                  intgroup.df, 
                  name = ddsTransform$sample)
  
  # Generate plot.
  
  pcaPlot <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group"))
  
  # Add labels to plot.
  
  pcaPlotLabels <- pcaPlot +
    geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() +
    geom_text_repel(aes(label = name), size = 2)
  
  # Format plot.
  
  formattedPCAPlot <- pcaPlotLabels + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Group plots into a single plot and annotate with tissue.
  
  figure <- ggarrange(formattedBoxPlot, formattedPCAPlot, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)
  
  figure <- annotate_figure(figure,
                            top = tissue)
  
  # Return figure.
  
  return(figure)
  
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsHead <- BuildDDSForDEAnalysis("head")

ddsFatBody <- BuildDDSForDEAnalysis("fatbody")

ddsOvaries <- BuildDDSForDEAnalysis("ovaries")

# Generate Exploratory Plots ----

headPlots <- GenerateExploratoryPlots("head")

fatBodyPlots <- GenerateExploratoryPlots("fatbody")

ovariesPlots <- GenerateExploratoryPlots("ovaries")

# Save Plots ----

# Head.

setwd("../02_outputs/02_head/02_quality_control_reports")

ggsave("20_NER0008751_obj3_exp3_fig_S1_head_exploratory_plots.svg",
       headPlots, width = 29.7, height = 21, units = "cm")

# Fat body.

setwd("../../03_fatbody/02_quality_control_reports")

ggsave("20_NER0008751_obj3_exp3_fig_S2_fat_body_exploratory_plots.svg",
       fatBodyPlots, width = 29.7, height = 21, units = "cm")

# Ovaries.

setwd("../../01_ovaries/02_quality_control_reports")

ggsave("20_NER0008751_obj3_exp3_fig_S3_ovaries_exploratory_plots.svg",
       ovariesPlots, width = 29.7, height = 21, units = "cm")
