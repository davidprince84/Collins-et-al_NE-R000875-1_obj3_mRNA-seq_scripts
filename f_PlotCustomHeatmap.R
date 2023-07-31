#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 06.05.2021
# File last updated: 08.06.2023
# Project: NER0008751 (Obj3) 
# Analysis: Comparative analysis with other species.
# Tasks: Produce a custom heatmap.
#-------------------------------------------------------------------------------
# Inputs:
# A matrix of gene names and expression values.

# Outputs:
# A heatmap with annotations.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(pheatmap)  # pheatmap().
library(ggplot2)  # ggsave(). 

# FUNCTION DEFINITION ----

PlotCustomHeatmap <- function (x) {
  # Plots a custom heat map, depending on the script calling the function.
  #
  # Args:
  #   x: number denoting which script is calling the function (41 or 45).
  #
  # Returns:
  #   A heatmap.
  
  # Set column gaps based on argument x.
  
  if (x == 41) {
    columnGaps <- c(3, 6)
  } else if (x == 45) {
    columnGaps <- c(2, 4)
  }
  
  # Set column and row annotations based on argument x.
  
  if (x == 41) {
    colAnnotations <- 
      data.frame(Treatment = c(NA, "M", "H", NA, "M", "H", "M", "H"),
                 Study = c("Pacifico et al. 2018", "Current study", "Current study",
                           "Chen et al. 2014", "Current study", "Current study",
                           "Current study", "Current study"),
                 Tissue = c("brain", "head", "head", "fat body", "fat body", 
                            "fat body", "ovaries", "ovaries"))

    row.names(colAnnotations) <- colnames(plotMatrix)
    
    rowAnnotations <- 
      data.frame(Longevity_influence = noZeroResults$longevity.influence)
    
    rownames(rowAnnotations) <- noZeroResults$name
    
    annotationColours <- list(Treatment = c("M" = "#0072B2", "H" = "#56B4E9"),
                              Study = c("Pacifico et al. 2018" = "darkorchid1", 
                                        "Current study" = "gray76",
                                        "Chen et al. 2014" = "orange"),
                              Tissue = c("brain" = "pink1", "head" = "bisque", "fat body" = "tan3", "ovaries" = "white"),
                              Longevity_influence =c("Pro-Longevity" = "white", "Anti-Longevity" = "black"))

  } else if (x == 45) {
    colAnnotations <- 
      data.frame(Treatment = c("M", "H", "M", "H", "M", "H"),
                 Tissue = c("head", "head", "fat body", 
                            "fat body", "ovaries", "ovaries"))
    
    row.names(colAnnotations) <- colnames(plotMatrix)
    
    annotationColours <- list(Treatment = c("M" = "#0072B2", "H" = "#56B4E9"),
                              Tissue = c("head" = "bisque", "fat body" = "tan3",
                                         "ovaries" = "white"))
  } else {
    stop('Argument x must equal 41 or 45.')
  }
  
  # Plot heatmap based on argument x.
  
  if (x == 41) {
    heatmapPlot <- pheatmap(plotMatrix,
                            scale = "none",
                            cluster_rows = TRUE,
                            cluster_cols = FALSE,
                            annotation_col = colAnnotations,
                            annotation_row = rowAnnotations,
                            annotation_colors = annotationColours,
                            legend = TRUE,
                            cellwidth=15,
                            cellheight=10,
                            fontsize = 7, 
                            gaps_col = columnGaps)
  } else if (x == 45) {
    heatmapPlot <- pheatmap(plotMatrix,
                            scale = "none",
                            cluster_rows = TRUE,
                            cluster_cols = FALSE,
                            annotation_col = colAnnotations,
                            annotation_colors = annotationColours,
                            legend = TRUE,
                            cellwidth=15,
                            cellheight=10,
                            fontsize = 7, 
                            gaps_col = columnGaps)
  }
  
  # Return heatmap.
  
  return(heatmapPlot)
  
}
