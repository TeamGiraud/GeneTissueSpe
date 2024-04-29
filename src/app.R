# Load packages ----
library(shiny)
library(shinydashboard)
library(dplyr)

# Load datasets of AIRE dependant genes
AIREdep = read.csv2("data/TRA_AIRE_dependency.csv")
# Load datasets of gene expression in mouse and human
gene_keys = read.delim2("data/Mouse_Human_merged_expression_data.txt", sep = " ")
colnames(gene_keys)[2] = "Gene.name"
gene_keys$nTPM = as.numeric(gene_keys$nTPM)
# load brain tissue list for the pool_brain function
brain_tissue = read.csv2("data/brain_tissues.csv") %>% pull(Tissue)

# Source helper functions -----
source("Genes_to_tissue_enrich.R", local = TRUE)


runApp()




