#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce Euler diagrams of the treatment differentially expressed gene
# (DEG) list comparisons.
#-------------------------------------------------------------------------------
# Inputs:
# Results of the treatment DEG list comparisons.

# Outputs:
# .svg file of the Euler diagrams.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(eulerr)  # euler()
library(ggplotify)  # as.ggplot()
library(ggplot2)
library(ggpubr)  # ggarrange()
library(RColorBrewer)  # source of Euler diagram fill colours.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Change Working Directory ----

setwd("../02_outputs/00_all_tissues/")

treatmentComparisonResults <- 
  read.csv("11_NER0008751_obj3_exp3_table_S9_comparison_stats.csv")

# FUNCTION DEFINITIONS ----

MakeEulerPlot <- function (x, y) {
  # Make a Euler plot of the results of comparing the DEGs between M and 
  # H treatments for a given tissue and direction of differential 
  # expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("fatbody", 
  #      "head" or "ovaries").
  #   y: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in RNA_2 compared to RNA_1) or down-
  #      regulated genes (genes more highly expressed in RNA_1 than RNA_2) 
  #      ("up" or "down").
  #
  # Returns:
  #   ggplot2 object containing a Euler diagram of the results.
  
  # Assign variable based on arguments.
  
  if (x == "head" && y == "up") {
    z <- 1
  } else if (x == "head" && y == "down") {
    z <- 2
  } else if (x == "fatbody" && y == "up") {
    z <- 3
  } else if (x == "fatbody" && y == "down") {
    z <- 4
  } else if (x == "ovaries" && y == "up") {
    z <- 5
  } else if (x == "ovaries" && y == "down") {
    z <- 6
  } else {
    stop('Arguments x and/or y are incorrect')
  }
  
  # Fit Euler Diagram to the data.
  
  eulerData <- euler(c(M = treatmentComparisonResults[z, 6],
                       H = treatmentComparisonResults[z, 7],
                       "M&H" = treatmentComparisonResults[z, 5]))
  
  # Plot Euler diagram.
  
  ePlot <- as.ggplot(plot(eulerData, 
                          fills = c("#0072B2", "#56B4E9", "white"),
                          quantities = TRUE))
  
  # Return Euler diagram.
  
  return(ePlot)
  
} 

# EXECUTED STATEMENTS ----

# Generate Euler Diagram of Each Comparison ----

# Head.

headUpPlot <- MakeEulerPlot("head", "up")

headDownPlot <- MakeEulerPlot("head", "down")

# Fat body.

fatBodyUpPlot <- MakeEulerPlot("fatbody", "up")

fatBodyDownPlot <- MakeEulerPlot("fatbody", "down")

# Ovaries.

ovariesUpPlot <- MakeEulerPlot("ovaries", "up")

ovariesDownPlot <- MakeEulerPlot("ovaries", "down")

# Combine Diagrams into Single Plot ----

plotsCombined <- ggarrange(headUpPlot, headDownPlot,
                           fatBodyUpPlot, fatBodyDownPlot,
                           ovariesUpPlot, ovariesDownPlot,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           ncol = 2, nrow = 3)

# Saving Combined Plot ----

# Save figures as SVG file so that the positioning of the labels
# can be adjusted manually for clarity.

ggsave("13_NER0008751_obj3_exp3_fig_4_euler_diagram.svg",
       plot = plotsCombined, height = 29.7, width = 21, units = "cm")
