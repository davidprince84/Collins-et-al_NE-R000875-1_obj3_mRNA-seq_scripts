#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj3) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Generate a data.frame which links Drosophila melanogaster (Dmel)
# gene IDs to gene names.
#-------------------------------------------------------------------------------
# Inputs:
# org.Dm.eg.db package.

# Outputs:
# Data.frame linking Drosophila melanogaster (Dmel) gene IDs to gene names.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(org.Dm.eg.db)  # select(), keys()

# LOAD DATA ----

# Isolating Dmel Genome Information from Annotation Package ----

# Place Entrezids keys in an object.

dmelKeys <- keys(org.Dm.eg.db, keytype = "ENTREZID")

# Select annotations based on the keys.

dmelAnnotations <-
  select(org.Dm.eg.db, columns = c("FLYBASE", "GENENAME"), keys = dmelKeys)

# Remove the ENTREZID column.

dmelAnnotations <- dmelAnnotations[, c("FLYBASE", "GENENAME")]

# Rename the columns.

colnames(dmelAnnotations) <- 
  c("Dmel_gene_ID", "Dmel_gene_description")

# Remove object no longer required.

rm(dmelKeys)
