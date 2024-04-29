library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)

# Tissue_spec : 
# input : a vector of gene names and the gene tissue expression dataset from HPA
# output : a dataframe with SPM and CTM values from Pan et al (2013) and Tau value from
# Yanai et al (2004) and a classification of gene tissue specificity

Tissue_spec = function(gene_list,  # chr vector with the gene to test tissue specificity
                       k = 6, #threhold of tissues for gap housekeeping gene
                       gene_keys, # dataset off tissue expression from HPA
                       spec_thres = 0.9,  # SPM threshold to classify gene as tissue restriced
                       select_thres_spm = 0.3, # SPM threshold to classify gene as tissue enriched
                       select_thres_ctm = 0.9,# CTM threshold to classify gene as tissue enriched
                       tau_thres = 0.7, # tau threshold to classify gene as tissue restricted
                       Ib_thres = 5  #threhold of tissues to consider for tissue enrichment
                       ){ 
# intersect between th 2
exprdat = gene_keys[gene_keys$Gene.name %in% gene_list,]

# log transform of the HPA dataset
#exprdat = exprdat %>% mutate(nTPM = nTPM * 10) # add a 10x to capture expression of low (<1) expressed genes
exprdat = exprdat %>% mutate(logTPM = ifelse(nTPM < 1, 1, log2(nTPM))) # pass to the log2 
### SPM calculation
# for each gene, calculate the sum of all squared gene expression in TPM in all samples
sum_squrd_expr =  exprdat %>%
  group_by(Gene.name) %>%
  summarise(sum_squrd_expr = sum(nTPM^2))
exprdat = left_join(exprdat, sum_squrd_expr, by = "Gene.name")
# calculate squared gene expression in TPM by sample
exprdat = exprdat %>% mutate(squrd_nTPM = nTPM^2)
# calculate SPM
exprdat = exprdat %>% mutate(SPM = squrd_nTPM/sum_squrd_expr)
# drop intermediate variables
exprdat = exprdat %>% dplyr::select(-squrd_nTPM,-sum_squrd_expr)
### CTM calculation
# calculate rank of expression in sammples by gene
exprdat <- exprdat %>% 
  arrange(Gene.name, desc(nTPM)) %>%
  group_by(Gene.name) %>% 
  mutate(rank = rank(desc(nTPM), ties.method = "first"))
# calculate CTM
varnames = c() # initiate the names for the CTM for each value of k 
for (i in (2:k)){
  varnames[i] = paste0("CTM_",i)
}
varnames = varnames[-c(1)] # remove first NA
for (i in (2:k)){   # calculer le CTM dans une nouvelle variable pour chaque valeur de k
  CTM = exprdat %>% 
    filter(rank <= i) %>%
    group_by(Gene.name) %>%
    summarise(CTM = sqrt(sum(SPM^2)))
  exprdat = left_join(exprdat,CTM, by = "Gene.name") 
}
colnames(exprdat)[(ncol(exprdat)-k+2):ncol(exprdat)] = varnames # renommer de maniere iterative
for (i in (2:(k))){   # mettre CTM = 0 si au dessus de k
  vec = data.frame(exprdat[,"rank"]<=i, exprdat[,varnames[i-1]])
  vec[,2] = ifelse(vec$rank, vec[,2], 0)
  exprdat[,paste0("CTM_",i)] = vec[,2]
}
# calculate max CTM
exprdat = exprdat %>% rowwise() %>%
  mutate(CTM=max(c_across(varnames[1]:last(varnames))))
# remove intermediary variables 
exprdat = exprdat %>% select(-varnames)

### Tau calculation
# calculate number of samples
exprdat = exprdat %>% count(Gene.name) %>% left_join(exprdat, by = "Gene.name")
# calculate the maximum value of gene TPM
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  summarise(max_logTPM = max(logTPM)) %>%
  left_join(exprdat, by = "Gene.name")
# calculate the Xi chapeau
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  mutate(Xi_chapeau = logTPM/max_logTPM) %>%
  mutate(un_moin_Xi_chapeau = 1 - Xi_chapeau)
# calculate the nominator
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  summarise(nomi = sum(un_moin_Xi_chapeau)) %>%
  left_join(exprdat, by = "Gene.name")
# calculate the Tau
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  mutate(Tau = nomi/(n-1))
# drop intermediate variables
exprdat = exprdat %>% dplyr::select(-nomi,-un_moin_Xi_chapeau, -max_logTPM, -Xi_chapeau )
#calculate gap
exprdat = exprdat %>% mutate(gap= -(nTPM - lag(nTPM, default = nTPM[1])))
# calculate the maximum value of gap by group
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  summarise(max_gap = max(gap)) %>%
  left_join(exprdat, by = "Gene.name")
# calculate the rank at which there is the gap
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  filter(max_gap == gap) %>%
  mutate(rank_gap = rank) %>%
  dplyr::select(Gene.name, rank_gap) %>%
  left_join(exprdat, by = "Gene.name")
# calculate the binary classifier if rank < rank_gap
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  mutate(bin = ifelse(rank_gap > Ib_thres, 0, ifelse(rank<rank_gap,1,0)))
# calculate binary index Ib
exprdat = exprdat %>%
  group_by(Gene.name) %>%
  summarise(Ib = sum(bin)) %>%
  left_join(exprdat, by = "Gene.name")
# drop useless variables
exprdat = exprdat %>% dplyr::select(-max_gap)
# add AIRE dep
exprdat = left_join(exprdat, filter(AIREdep, !is.na(AIREdep)), by = "Gene.name")
exprdat = exprdat %>% mutate(AIREdep = ifelse(is.na(AIREdep), "unknown", AIREdep))

### Class genes:

#exprdat = exprdat %>% mutate(Enrich_SPM = ifelse(SPM > spec_thres, "Tissue_restricted", ifelse( (SPM>select_thres_spm)&(CTM>select_thres_ctm), "Tissue_enriched", "Low_spe")))
# classify as tissue restricted is tau is above tau_thr
# classify as tissue enriched if 
#exprdat = exprdat %>% mutate(Enrich_Tau = ifelse((Tau>tau_thres)&(rank == 1),"Tissue_restricted", ifelse((Tau>(tau_thres-0.1))&(rank <= Ib_thres),"Tissue_enriched", "Low_spe")))

#return exprdat
return(exprdat)
}




### Plot_tissue_expr
#input : data: dataframe returned by the above function with gene tissue expression
#        gene : chr with the gene name to plot
# output : ggplot ordered barchart of the tissue expression
plot_tissue_expr = function(data, gene){
  df = data %>% filter(Gene.name == gene, !is.na(Tissue))  %>%
  arrange(desc(nTPM)) %>%
  ungroup()
  # plot en barchart
  ggplot(data=df, aes(reorder(Tissue, nTPM), nTPM, fill = Enrich_Tau)) +
    geom_bar(stat="identity", width=0.5) +
    xlab("Tissue") +
    ggtitle(gene) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=10))
}


# brouillon
plot_Ib = function(data){
  temp = data %>% distinct(Gene.name, Ib)
  p = ggplot(temp, aes(Ib))  + 
    geom_histogram() +
    ggtitle("Ib distribution") +
    theme_prism(base_size = 16) 
  return(p)
}

# brouillon
plot_tau = function(data){
temp = data %>% distinct(Gene.name, Tau)
p = ggplot(temp, aes(Tau))  + 
  geom_histogram() +
  ggtitle("Tau distribution") +
  theme_prism(base_size = 16) 
return(p)
}


plot_by_tissue = function(data, TRA_or_TSA = "TRA", Tau_or_SPM = "Tau"){
  if (Tau_or_SPM == "Tau"){
    data = data %>% filter(Enrich_Tau == ifelse(TRA_or_TSA == "TRA", "Tissue_restricted", "Tissue_enriched")) %>%
      ungroup() %>% count(Tissue) %>% filter(Tissue != "thymus")
  }else{data = data %>% filter(Enrich_SPM == ifelse(TRA_or_TSA == "TRA", "Tissue_restricted", "Tissue_enriched")) %>%
    ungroup() %>% count(Tissue)%>% filter(Tissue != "thymus") 
  }
  p = ggplot(data, aes(y=reorder(Tissue, n),x = n)) +
    geom_bar(stat="identity", width=0.5) +
    ylab("") +
    ggtitle(paste0(TRA_or_TSA," by Tissue with ",Tau_or_SPM," method")) +
    scale_x_continuous(expand=c(0,0)) +
    theme_prism(base_size = 16) + 
    theme(legend.position = "none", axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size=10)) 
  return(p)
}


import_hpa = function(gene_list){
  gene_keys <- data.frame()
  for (i in 1:length(gene_list)){
    temp = gene_list[i]
    gene_keys <- rbind(gene_keys, fromJSON(getURL(paste0("www.proteinatlas.org/api/search_download.php?search=",temp,"&format=json&columns=g,t_RNA_adipose_tissue,t_RNA_adrenal_gland,t_RNA_amygdala,t_RNA_appendix,t_RNA_basal_ganglia,t_RNA_bone_marrow,t_RNA_breast,t_RNA_cerebellum,t_RNA_cerebral_cortex,t_RNA_cervix,t_RNA_choroid_plexus,t_RNA_colon,t_RNA_duodenum,t_RNA_endometrium_1,t_RNA_epididymis,t_RNA_esophagus,t_RNA_fallopian_tube,t_RNA_gallbladder,t_RNA_heart_muscle,t_RNA_hippocampal_formation,t_RNA_hypothalamus,t_RNA_kidney,t_RNA_liver,t_RNA_lung,t_RNA_lymph_node,t_RNA_medulla_oblongata,t_RNA_midbrain,t_RNA_ovary,t_RNA_pancreas,t_RNA_parathyroid_gland,t_RNA_pituitary_gland,t_RNA_placenta,t_RNA_pons,t_RNA_prostate,t_RNA_rectum,t_RNA_retina,t_RNA_salivary_gland,t_RNA_seminal_vesicle,t_RNA_skeletal_muscle,t_RNA_skin_1,t_RNA_small_intestine,t_RNA_smooth_muscle,t_RNA_spinal_cord,t_RNA_spleen,t_RNA_stomach_1,t_RNA_testis,t_RNA_thalamus,t_RNA_thymus,t_RNA_thyroid_gland,t_RNA_tongue,t_RNA_tonsil,t_RNA_urinary_bladder,t_RNA_vagina,t_RNA_white_matter&compress=no"))))
  }
  gene_keys = filter(gene_keys, Gene %in% gene_list)
  for (i in 3:dim(gene_keys)[2]){
    colnames(gene_keys)[i] <-  gsub(" \\[.*", "", gsub(".*- ","",colnames(gene_keys[i])))
  }
  gene_keys <- pivot_longer(gene_keys, names_to = "Tissue", values_to = "nTPM", cols = c(3:dim(gene_keys)[2]))
  colnames(gene_keys)[2] <- "Gene.name"
  gene_keys$nTPM = as.numeric(gene_keys$nTPM)
  return(gene_keys)
}


# change the brain samples
change_brain = function(data, pool_brain = T){
  if (pool_brain){
    data_b = data %>% filter(Tissue == "brain") %>%
    group_by(Gene.name) %>%
    mutate(rank_brain = min(rank)) %>%
    mutate(Tissue2 = ifelse(rank_brain == rank, "brain", NA))
  data3  = left_join(data, data_b)
  data3 = data3 %>% mutate(Tissue3 = ifelse(is.na(Tissue2)&is.na(rank_brain), Tissue, Tissue2))
 data3 = data3 %>% mutate(Tissue = Tissue3)
  data3 = data3 %>% select(-Tissue3, -rank_brain, -Tissue2)
  }
  else{data3 = data}
  return(data3)
}


