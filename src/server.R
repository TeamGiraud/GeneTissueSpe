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


# Server logic ----
server <- function(input, output) {
  
  #Store Initial genes from User Inputs
  values <- reactiveValues(gene_list = NULL)
  
  observeEvent(input$action, {
    values$gene_list <- unlist(strsplit(c(input$gene_list), split = " "))
  })
  
  ds1 = reactive({
    Tissue_spec(
      gene_keys = gene_keys,
      gene_list = values$gene_list,
      k=input$k, 
      spec_thres = input$spec_thres,  
      select_thres_spm = input$select_thres_spm, 
      select_thres_ctm = input$select_thres_ctm,
      tau_thres = input$tau_thres, 
      Ib_thres = input$Ib_thres) %>% 
      mutate(SPM = round(SPM,5)) %>%
      mutate(CTM = round(CTM,5)) %>%
      mutate(Tau = round(Tau,5)) %>%
      mutate(Enrich_SPM = ifelse(SPM >= as.numeric(input$spec_thres), "Tissue_restricted",
                                 ifelse( (SPM < as.numeric(input$spec_thres)) & (SPM>as.numeric(input$select_thres_spm)) & (CTM>as.numeric(input$select_thres_ctm)), "Tissue_enriched",
                                         "Low_spe"))) %>%
      mutate(Enrich_Tau = ifelse( Tau >= as.numeric(input$tau_thres),
                                 ifelse( Ib == 1 & rank == 1, "Tissue_restricted",
                                        ifelse( (Ib > 1) & (Ib <= as.numeric(input$Ib_thres)) & (bin == 1), "Tissue_enriched", 
                                               "Low_spe")),
                                 ifelse( (Tau<as.numeric(input$tau_thres)) & (Tau>(as.numeric(input$tau_thres)-0.15) ),
                                         ifelse( (Ib <= as.numeric(input$Ib_thres)) & (bin == 1),
                                                 ifelse( Ib == 1 & rank == 1, "Tissue_restricted","Low_Spe"),
                                                 "Low_spe"),
                                         "Low_spe"))
      ) %>%
      mutate(Tissue = ifelse(input$pool_brain == T & Tissue %in% brain_tissue,"brain", Tissue )) %>%
      group_by(Gene.name)
  })
  
  ds = reactive({change_brain(data = ds1(), pool_brain = input$pool_brain)})
  
  output$table <- renderDataTable(
    ds()
  )
  
  output$TEST <- renderPrint({
    values$gene_list
  })
  
  output$plottau <- renderPlot({
    plot_tau(ds())
  }, res = 96)
  
  output$plotib <- renderPlot({
    plot_Ib(ds())
  }, res = 96)
  
  output$box_len_genlist <- renderText({
    paste0("number of genes in list :   ", length(values$gene_list)) 
  })
  
  stats_restr = reactive({ round((ds() %>% group_by(Gene.name) %>% count(Enrich_Tau) %>% filter(Enrich_Tau == "Tissue_restricted") %>% nrow())/length(values$gene_list),3)})
  stats_enrich = reactive({ round((ds() %>% group_by(Gene.name) %>% count(Enrich_Tau) %>% filter(Enrich_Tau == "Tissue_enriched") %>% nrow())/length(values$gene_list),3)})
  Aire_dep = reactive({ (ds() %>% group_by(Gene.name) %>% count(AIREdep) %>% ungroup() %>% pull(AIREdep) %>% as.factor() %>% table() %>% prop.table())[["yes"]] %>% round(3)})
  n_tau_genes = reactive({ds() %>% dplyr::select(Tau) %>% unique() %>% filter(Tau > input$tau_thres) %>% count(Tau) %>% pull(Tau) %>% length() %>% as.character()  })
  
  output$stats_enrich <- renderText({
    stats_enrich()
  })
  
  output$stats_restr <- renderText({
    stats_restr()
  })
  
  output$n_tau_genes <- renderText({
    n_tau_genes()
  })
  
  output$Aire_dep <- renderText({
    Aire_dep()
  })
  
  values <- reactiveValues(single_gene = NULL)
  observeEvent(input$action, {
    values$single_gene <- input$single_gene
  })
  
  output$plotgene <- renderPlot({
    plot_tissue_expr(data =  ds(), gene = values$single_gene)
  }, res = 96)
  
  output$plottissue1 <- renderPlot({
    plot_by_tissue(data = ds(), TRA_or_TSA = "TRA", Tau_or_SPM = input$method1)
  }, res = 96)
  
  output$plottissue2 <- renderPlot({
    plot_by_tissue(data = ds(), TRA_or_TSA = "TSA", Tau_or_SPM = input$method2)
  }, res = 96)
  
  output$downloaddata <- downloadHandler(
    filename = "TissueSpecResults.csv",
    content = function(file) {
      write.csv(ds(), file)
    }
  )
  
  
}