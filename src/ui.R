# Load packages ----
library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)
library(shinythemes)
library(googlesheets4)

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


# User interface ----
### Experimental section

header = dashboardHeader(title = "GeneTissueSpec")

sidebar = dashboardSidebar(sidebarMenu(
  menuItem("Graphs", tabName = "dashboard", icon = icon("chart-scatter")),
  menuItem("Results table", icon = icon("table-cells"), tabName = "widgets",
           badgeLabel = "data", badgeColor = "green")
),
helpText("Welcome to GTS, a shiny-based application to identify gene tissue specificity and enrichment.
                            version 0.7"),
title = "Parameters :",

selectInput("specie_choice",
            label = "specie :",
            choices = list("Human (hsapiens)" = "Human","Mouse (mmusculus" = "Mouse"),
            selected = "Human",
            width = "200px"),

textInput("gene_list",
          label = "gene list :",
          value = "INS KRT5 ACT2 CHRNA1 TNNT2 TSPAN18 ALB"),

actionButton("action", "Submit genes", width = "200px", icon = icon("refresh")), 
submitButton(text = "Refresh page", width = "200px", icon = icon("refresh")),

helpText("max value for CTM calculation (see Yanai et al, 2004) :"),
numericInput("k",
             label = "k :",
             value = 6,
             width = "150px"),

helpText("SPM threshold over which a gene is classified as tissue restricted (see Pan et al, 2013) :"),
numericInput("spec_thres",
             label = "SPM restricted :",
             value = 0.9,
             width = "150px"),

helpText("SPM threshold over which a gene is classified as tissue enriched (see Pan et al, 2013) :"),
numericInput("select_thres_spm",
             label = "SPM enriched :",
             value = 0.3,
             width = "150px"),

helpText("CTM threshold over which a gene is classified as tissue enriched (see Pan et al, 2013) :"),
numericInput("select_thres_ctm",
             label = "CTM thres :",
             value = 0.9,
             width = "150px"),

helpText("Tau threshold over which a gene is classified as tissue restricted (see Yanai et al, 2004) :"),
numericInput("tau_thres",
             label = "tau thres :",
             value = 0.8,
             width = "150px"),

helpText("number of tissues in which a gene can be expressed over which it is considered non selected (see Yanai et al, 2004) :"),
numericInput("Ib_thres",
             label = "binary index thres :",
             value = 8,
             width = "150px"),

helpText("Because of the multiplicity of brain sub-tissues in the gene expression database, do you want to
                                   simplify using a common 'brain' label?"),
radioButtons("pool_brain", label = "pool all brain tissues :",
             choices = c(T,F))
)




body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dashboard",
            fluidRow(
              box(
                title = "Tau index distribution in gene list", width = 6, solidHeader = TRUE, status = "primary",
                plotOutput("plottau"),
              ),
              box(
                title = "Binary index distribution in gene list", width = 6, solidHeader = TRUE, status = "primary",
                plotOutput("plotib"),
              )
            ),
            fluidRow(
              box(
                title = "Tissue restricted genes repartition", width = 6, solidHeader = TRUE, status = "primary",
                plotOutput("plottissue1"), radioButtons("method1", label = "method :", choices = c("Tau", "SPM"))
              ),
              box(
                title = "Tissue enriched genes repartition", width = 6, solidHeader = TRUE, status = "primary",
                plotOutput("plottissue2"), radioButtons("method2", label = "method :", choices = c("Tau", "SPM"))
              ),
            ),
            fluidRow(
              column(width = 6,
                     box(
                       title = "Tissue expression", width = NULL, solidHeader = TRUE, status = "primary",
                       plotOutput("plotgene"),textInput("single_gene",
                                                        label = "gene to plot :",
                                                        value = "INS")
                     )
              ),
              column(width = 3,
                     box(
                       title = "Number of genes :", width = NULL, background = "blue",
                       textOutput("box_len_genlist")
                     ),
                     box(
                       title = "Number of genes with Tau above chosen threshold :", width = NULL, background = "blue",
                       textOutput("n_tau_genes")
                     ),
                     box(
                       title = "Proportion of TRA in list (Tau) :", width = NULL, background = "maroon",
                       textOutput("stats_restr")
                     ),
                     box(
                       title = "Proportion of TSA in list (Tau) :", width = NULL, background = "orange",
                       textOutput("stats_enrich")
                     ),
                     box(
                       title = "Proportion of Aire dependent genes :", width = NULL, background = "red",
                       textOutput("Aire_dep")
                     )
              )
            )
            
    ),
    tabItem(tabName = "widgets",
            box(
              title = "results dataset", width = 12, solidHeader = TRUE, status = "primary",
              dataTableOutput("table"),
            ) ,
            box(downloadButton("downloaddata", label = "Download results"), width = 2)
    )
  )
)



ui = dashboardPage(
  header,
  sidebar,
  body
)