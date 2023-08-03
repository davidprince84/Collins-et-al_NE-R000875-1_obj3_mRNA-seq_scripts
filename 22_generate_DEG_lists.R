#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce lists of differentially expressed genes (DEGs) for each tissue
# and treatment.
#-------------------------------------------------------------------------------
# Inputs:
# Dmel_6_22_97_transcripts2genes.txt file, dmelAnnotations object and Kallisto 
# abundances.

# Outputs:
# .csv files containing lists of DEGs for a given tissue and treatment.
# .csv file containing all genes expressed in the analysis, to be used as a
# background list in downstream comparisons.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
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

# LOAD CUSTOM FUNCTIONS ----

source("f_BuildDDSForDEAnalysis.R")

# FUNCTION DEFINITIONS ----

MakeGeneList <- function (x, LFC = 0) {
  # Produces and saves a list of all expressed genes for a tissue with 
  # differential expression results added.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("fatbody", "head"  
  #      or "ovaries").
  #   LFC: numeric denoting the LFC to use when filter DESeq2 results
  #      (default of 0 (i.e. all significant genes)).
  #
  # Returns:
  #   A .csv file of Drosophila melanogaster genes with associated normalized 
  #   counts for each sample and differential expression results.
  
  # Note: false discovery rate (FDR) of 0.05 used, but this can be changed by
  # changing the alpha level in the results() call below.
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Set variables based on arguments.
  
  if (x == "head") {
    dds <- ddsHead
    tableNumber <- "D1"
  } else if (x == "fatbody") {
    dds <- ddsFatBody
    tableNumber <- "D2"
  } else if (x == "ovaries") {
    dds <- ddsOvaries
    tableNumber <- "D3"
  } else {
    stop('Argument x must be "fatbody", "head" or "ovaries".')
  }
  
  # Generate data.frame of normalized counts for all genes and add gene names
  # as a column (for merger with differential expression results).
  
  normalizedCounts <- as.data.frame(counts(dds, normalized = TRUE))
  
  normalizedCounts$gene_symbol <- row.names(normalizedCounts)
  
  # Extract differential expression results for M (medium) and add gene names
  # as a column (for merger with normalized count data).
  
  allMExtractedResults <- 
    as.data.frame(results(dds, contrast = c("condition", "M_RNA2", "M_RNA1"),
                          alpha = 0.05, lfcThreshold = LFC))
  
  allMExtractedResults$gene_symbol <- row.names(allMExtractedResults)
  
  # Add columns explicitly stating where a gene is differentially expressed and
  # whether it is up-regulated (with age, i.e. more expressed in old flies)
  # or down-regulated (with age, i.e. more expressed in younger flies).
  
  allMExtractedResults[, "DEG"] <- NA
  allMExtractedResults[, "expression_with_age"] <- NA
  
  for (i in (1:length(allMExtractedResults[, 1]))) {
    if (is.na(allMExtractedResults[i, "padj"])) {
      allMExtractedResults[i, "DEG"] <- "NO"
    } else if (allMExtractedResults[i, "padj"] < 0.05) {
      allMExtractedResults[i, "DEG"] <- "YES"
    } else {
      allMExtractedResults[i, "DEG"] <- "NO"
    }
    if (allMExtractedResults[i, "DEG"] == "YES" && allMExtractedResults[i, "log2FoldChange"] > 0) {
      allMExtractedResults[i, "expression_with_age"] <- "up-regulated"
    } else if (allMExtractedResults[i, "DEG"] == "YES" && allMExtractedResults[i, "log2FoldChange"] < 0) {
      allMExtractedResults[i, "expression_with_age"] <- "down-regulated"
    } else {
      allMExtractedResults[i, "expression_with_age"] <- "NA"
    }
  }
  
  # Rename columns.
  
  colnames(allMExtractedResults) <- 
    c("baseMean_M", "log2FoldChange_M", "lfcSE_M", "stat_M", "pvalue_M",
      "padj_M", "gene_symbol", "DEG_in_M", "expression_with_age_in_M")
  
  # Extract differential expression results for H (high) and add gene names
  # as a column (for merger with normalized count data).
  
  allHExtractedResults <- 
    as.data.frame(results(dds, contrast = c("condition", "H_RNA2", "H_RNA1"),
                          alpha = 0.05, lfcThreshold = LFC))
  
  allHExtractedResults$gene_symbol <- row.names(allHExtractedResults)
  
  # Add columns explicitly stating where a gene is differentially expressed and
  # whether it is up-regulated (with age, i.e. more expressed in old bees)
  # or down-regulated (with age, i.e. more expressed in younger bees).
  
  allHExtractedResults[, "DEG"] <- NA
  allHExtractedResults[, "expression_with_age"] <- NA
  
  for (i in (1:length(allHExtractedResults[, 1]))) {
    if (is.na(allHExtractedResults[i, "padj"])) {
      allHExtractedResults[i, "DEG"] <- "NO"
    } else if (allHExtractedResults[i, "padj"] < 0.05) {
      allHExtractedResults[i, "DEG"] <- "YES"
    } else {
      allHExtractedResults[i, "DEG"] <- "NO"
    }
    if (allHExtractedResults[i, "DEG"] == "YES" && allHExtractedResults[i, "log2FoldChange"] > 0) {
      allHExtractedResults[i, "expression_with_age"] <- "up-regulated"
    } else if (allHExtractedResults[i, "DEG"] == "YES" && allHExtractedResults[i, "log2FoldChange"] < 0) {
      allHExtractedResults[i, "expression_with_age"] <- "down-regulated"
    } else {
      allHExtractedResults[i, "expression_with_age"] <- "NA"
    }
  }
  
  # Rename columns.
  
  colnames(allHExtractedResults) <- 
    c("baseMean_H", "log2FoldChange_H", "lfcSE_H", "stat_H", "pvalue_H",
      "padj_H", "gene_symbol", "DEG_in_H", "expression_with_age_in_H")
  
  # Merge normalizedCounts, allMExtractedResults and allHExtractedResults to 
  # attach count data to differential expression results.
  
  countsAndResults <- merge(normalizedCounts, allMExtractedResults, by = "gene_symbol")
  
  countsAndResults <- merge(countsAndResults, allHExtractedResults, by = "gene_symbol")
  
  # Merge countsAndResults with dmelAnnotations to add gene names.
  
  colnames(dmelAnnotations) <- c("gene_symbol", "gene_name")
  
  namedResults <- merge(countsAndResults, dmelAnnotations, by.x = "gene_symbol",
                        all.x = TRUE)
  
  # Write namedResults to file for manuscript supplementary data and 
  # submission to Gene Expression Omnibus (GEO).
  
  write.csv(namedResults, 
            paste0(x, "_all_expressed_genes_and_DE_results.csv"),
            row.names = FALSE)
  
  write.csv(namedResults, 
            paste0("00_NER0008751_obj1_exp1_table_", tableNumber, ".csv"),
            row.names = FALSE)
}

MakeDEGLists <- function (x) {
  # Make DEG lists for a tissue by dividing up the list of all genes into
  # lists of DEGs.
  #
  # Args:
  #   x: string denoting which tissue to generate a DEG list for ("fatbody",
  #      "head" or "ovaries").
  #
  # Returns:
  #   Four .csv files of Drosophila melanogaster DEGs and associated stats.
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.

  # Load overall gene list based on argument.
  
  overallGeneList <- 
    read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Split overallGeneList into individual lists.
  
  # First, split by treatment, rename columns and remove NAs in expression_with_age.
  
  columnNames <- c("symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", 
                   "pvalue", "padj", "DEG", "expression_with_age", "name")
  
  mResults <- overallGeneList[, c(1, 14:21, 30)]
  
  colnames(mResults) <- columnNames
  
  mResults <- mResults[!(is.na(mResults[, "expression_with_age"])), ]
  
  hResults <- overallGeneList[, c(1, 22:30)]
  
  colnames(hResults) <- columnNames
  
  hResults <- hResults[!(is.na(hResults[, "expression_with_age"])), ]
  
  # Split each treatment into up- and down-regulated DEGs.
  # Up-regulated = more expressed in older flies.
  # Down-regulated = more expressed in younger flies.
  
  mUpDEGs <- 
    mResults[mResults[, "expression_with_age"] == "up-regulated", ]
  
  mDownDEGs <- 
    mResults[mResults[, "expression_with_age"] == "down-regulated", ]
  
  hUpDEGs <- 
    hResults[hResults[, "expression_with_age"] == "up-regulated", ]
  
  hDownDEGs <- 
    hResults[hResults[, "expression_with_age"] == "down-regulated", ]
  
  # Remove DEG and expression_with_age columns.
  
  mUpDEGs <- mUpDEGs[, c(1:7, 10)]
  
  mDownDEGs <- mDownDEGs[, c(1:7, 10)]
  
  hUpDEGs <- hUpDEGs[, c(1:7, 10)]
  
  hDownDEGs <- hDownDEGs[, c(1:7, 10)]
  
  # Write results to .csv files.
  
  write.csv(mUpDEGs, paste0(x, "_results_M_treatment_up_regulated_LFC0.csv"),
            row.names = FALSE)
  
  write.csv(mDownDEGs, paste0(x, "_results_M_treatment_down_regulated_LFC0.csv"),
            row.names = FALSE)
  
  write.csv(hUpDEGs, paste0(x, "_results_H_treatment_up_regulated_LFC0.csv"),
            row.names = FALSE)
  
  write.csv(hDownDEGs, paste0(x, "_results_H_treatment_down_regulated_LFC0.csv"),
            row.names = FALSE)
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsHead <- BuildDDSForDEAnalysis("head")

ddsFatBody <- BuildDDSForDEAnalysis("fatbody")

ddsOvaries <- BuildDDSForDEAnalysis("ovaries")

# Make and Save Gene Lists ----

# Head.

dir.create("../02_outputs/02_head/31_DESeq2_DEG_lists")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../02_outputs/02_head/31_DESeq2_DEG_lists")

# Make gene lists.

MakeGeneList("head")

MakeDEGLists("head")

# Fat body.

# Create directory.

dir.create("../../03_fatbody/31_DESeq2_DEG_lists")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists")

# Make gene lists.

MakeGeneList("fatbody")

MakeDEGLists("fatbody")

# Ovaries.

# Create directory.

dir.create("../../01_ovaries/31_DESeq2_DEG_lists")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists")

# Make gene lists.

MakeGeneList("ovaries")

MakeDEGLists("ovaries")
