## This part should follow the previous part which contains essential function and data

##------- loading the package -------
## gene name and ID annotation
library(org.Hs.eg.db)
library(biomaRt)  ## for gene name/id annotation
## Tidyverse and supporting packages for Tidyverse
library(dplyr); library(tidyr); library(tibble); library(stringr); library(readr)
library(ggplot2); library(ggpubr); library(ggrepel)
library(patchwork)  ## plot multiple plot together by direct adding
library(gridExtra)  ## multiple plots
library(scales)  ## axis scales
library(GGally); library(extrafont)
## Machine Learning
library(caret) ## ML assistant tool
## survival analysis
library(survival); library(survminer)

########## Machine Learning ################


##------- TCGA data pretreatment: fixing the gene name and ID in TCGA expression data ----------
## TCGA.PCA18.cbp_exp.rsem_RAW is a version whose gene name and ID issue have been resolved by the following code:
## Hugo symbol of 642846 will be LOC642846 according to NCBI
## Hugo symbol of 9403 will be SELENOF according to NCBI
for(i in 1:length(TCGA.PCA18.cbp_exp.rsem_RAW)) {
  if(sum(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id == 642846) >1) {stop("Duplicated target gene ID!")}
  
  if (TCGA.PCA18.cbp_exp.rsem_RAW[[i]][TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id == 642846,"Hugo_Symbol"] == "DDX12B") {
    TCGA.PCA18.cbp_exp.rsem_RAW[[i]][TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id == 642846,"Hugo_Symbol"] <- "LOC642846"
  } else {warning("Name unchanged on ID 642846")}
  if (TCGA.PCA18.cbp_exp.rsem_RAW[[i]][TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id == 9403,"Hugo_Symbol"] == "SEP15") {
    TCGA.PCA18.cbp_exp.rsem_RAW[[i]][TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id == 9403,"Hugo_Symbol"] <- "SELENOF"
  } else {warning("Name unchanged on ID 9403")}
}

## ID of KIAA1107 is 284697 (23285 is out of date)
for(i in 1:length(TCGA.PCA18.cbp_exp.rsem_RAW)) {
  if(sum(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol == "KIAA1107") >1) {stop("Duplicated target gene ID!")}
  if (TCGA.PCA18.cbp_exp.rsem_RAW[[i]][TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol == "KIAA1107","Entrez_Gene_Id"] == 23285) {
    TCGA.PCA18.cbp_exp.rsem_RAW[[i]][TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol == "KIAA1107","Entrez_Gene_Id"] <- 284697
  } else {warning("Name unchanged")}
}

## ID of QSOX1 is 5768 (200058 is out of date)
for(i in 1:length(TCGA.PCA18.cbp_exp.rsem_RAW)) {
  ## location of QSOX1
  ind <- which(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol == "QSOX1")
  
  if(length(ind) == 1) {  ## duplication does not exist
    warning("single QSOX1 gene.")
    ## check which version of QSOX1 in this case
    if(paste(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol[ind],TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id[ind],sep = "_") == "QSOX1_200058"){
      warning(paste("contain only QSOX1_200058 at", names(TCGA.PCA18.cbp_exp.rsem_RAW)[[i]],sep = " ")); next
    } else if (paste(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol[ind],TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id[ind],sep = "_") == "QSOX1_5768"){
      warning(paste("unchanged: contain only QSOX1_5768 at", names(TCGA.PCA18.cbp_exp.rsem_RAW)[[i]],sep = " ")); next
    } else {stop("Condition unknown. Please check your code.")}
    
  } else if (length(ind) == 2) {  ## in the case that we expect: two duplications
    ## the one we want to discard: the one that is not QSOX1_5768
    if (!("QSOX1_5768" %in% paste(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol[ind],TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id[ind],sep = "_"))) {stop("QSOX1_5768 is missing!")}
    ind2 <- paste(TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Hugo_Symbol[ind],TCGA.PCA18.cbp_exp.rsem_RAW[[i]]$Entrez_Gene_Id[ind],sep = "_") != "QSOX1_5768"
    TCGA.PCA18.cbp_exp.rsem_RAW[[i]] <- TCGA.PCA18.cbp_exp.rsem_RAW[[i]][-ind[ind2],]
  } else {warning(paste("more than two QSOX1 at", names(TCGA.PCA18.cbp_exp.rsem_RAW)[[i]],sep = " "))} ## the case of more than 2 duplications
}

## ready-to-be-used-formate TCGA RNA expression
## make a working version
TCGA.PCA18.cbp_exp.rsem <- TCGA.PCA18.cbp_exp.rsem_RAW

for (i in 1:length(TCGA.PCA18.cbp_exp.rsem)){
  ## preprocess on working_table_TCGA
  TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Entrez <- paste(TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol, TCGA.PCA18.cbp_exp.rsem[[i]]$Entrez_Gene_Id, sep = "_")  ## combine Hugo and id for later name matching
  TCGA.PCA18.cbp_exp.rsem[[i]] <- relocate(TCGA.PCA18.cbp_exp.rsem[[i]], Hugo_Entrez, .after = "Entrez_Gene_Id")  ## move it to the front
  ## check duplication. Name and ID are identical at the same time
  if (any(duplicated(TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Entrez))) {
    duplicated_TCGA <- TCGA.PCA18.cbp_exp.rsem[[i]][TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol %in% TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol[which(duplicated(TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Entrez))],]
  }
  ## after preview on the duplication, remove them from TCGA
  TCGA.PCA18.cbp_exp.rsem[[i]] <- TCGA.PCA18.cbp_exp.rsem[[i]][-which(TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol %in% TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol[which(duplicated(TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Entrez))]),]
}


##------- TCGA data pretreatment: gene identifier duplication ---------
## duplication on Hugo Symbol
duplicated_TCGA <- TCGA.PCA18.cbp_exp.rsem[[i]][TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol %in% TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol[duplicated(TCGA.PCA18.cbp_exp.rsem[[i]]$Hugo_Symbol)],]
duplicated_TCGA <- duplicated_TCGA[order(duplicated_TCGA$Hugo_Symbol),]

## duplication on Entrez ID
duplicated_TCGA <- TCGA.PCA18.cbp_exp.rsem$LUAD[TCGA.PCA18.cbp_exp.rsem$LUAD$Entrez_Gene_Id %in% TCGA.PCA18.cbp_exp.rsem$LUAD$Entrez_Gene_Id[duplicated(TCGA.PCA18.cbp_exp.rsem$LUAD$Entrez_Gene_Id)],]
duplicated_TCGA <- duplicated_TCGA[order(duplicated_TCGA$Entrez_Gene_Id),]

##------- CCLE data pretreatment: data cleaning -----------
## this aims to deal with the problematic gene items in CCLE
table(duplicated(meta_gene_GExp$Entrez_Gene_Id))
## pick the genes that are problematic
genelist_CCLE_Exp_problematic <- meta_gene_GExp[meta_gene_GExp$Entrez_Gene_Id %in% meta_gene_GExp$Entrez_Gene_Id[duplicated(meta_gene_GExp$Entrez_Gene_Id)],]
genelist_CCLE_Exp_problematic <- genelist_CCLE_Exp_problematic[order(genelist_CCLE_Exp_problematic$Entrez_Gene_Id),]
genelist_CCLE_Exp_problematic$Hugo_Entrez <- paste(genelist_CCLE_Exp_problematic$Hugo_Symbol,genelist_CCLE_Exp_problematic$Entrez_Gene_Id,sep = "_")
## location of problematic gene iterms in CCLE data, according to the order of genelist
ind <- match(genelist_CCLE_Exp_problematic$Hugo_Symbol, meta_gene_GExp$Hugo_Symbol)
## review on prolematics gene and their data
RNA_expression_problematic <- RNA_expression[,(ind+1)]

## a view from correlation significance
summary_RNA_expression_problematic <- summary_RNA_expression[summary_RNA_expression$Gene_Transcript %in% genelist_CCLE_Exp_problematic$Hugo_Symbol,]
## correct the ID of RLN2 in order to have a proper selection from later
genelist_CCLE_Exp_problematic$Entrez_Gene_Id[genelist_CCLE_Exp_problematic$Hugo_Symbol == "RLN2"] <- 6019
## from the duplicates items, we pick the genes with a proper official Hugo symbol
library(org.Hs.eg.db)
genelist_CCLE_Exp_problematic$Hugo_Symbol_official <- mapIds(org.Hs.eg.db,
                                                             keys= as.character(genelist_CCLE_Exp_problematic$Entrez_Gene_Id),
                                                             column="SYMBOL",
                                                             keytype="ENTREZID",
                                                             multiVals="first") %>% as.character()
## a preview before removing
genelist_CCLE_Exp_problematic[genelist_CCLE_Exp_problematic$Hugo_Symbol == genelist_CCLE_Exp_problematic$Hugo_Symbol_official,]
genelist_CCLE_Exp_problematic[genelist_CCLE_Exp_problematic$Hugo_Symbol != genelist_CCLE_Exp_problematic$Hugo_Symbol_official,]

genelist_CCLE_Exp_remove <- genelist_CCLE_Exp_problematic[genelist_CCLE_Exp_problematic$Hugo_Symbol != genelist_CCLE_Exp_problematic$Hugo_Symbol_official,]
## make sure the number of items is correct
dim(genelist_CCLE_Exp_problematic); dim(genelist_CCLE_Exp_remove)

## a preview on the genes that are going to be remove from CCLE
colnames(RNA_expression)[which(colnames(RNA_expression) %in% genelist_CCLE_Exp_remove$Hugo_Entrez)]; dim(RNA_expression)
## remove
RNA_expression <- RNA_expression[,-which(colnames(RNA_expression) %in% genelist_CCLE_Exp_remove$Hugo_Entrez)]; dim(RNA_expression)

## resolve the mistake on gene ID of RLN2, which is 6019 instead of 6013
colnames(RNA_expression)[colnames(RNA_expression) == "RLN2_6013"] <- "RLN2_6019"

## generate the new gene list
meta_gene_GExp <- data.frame(Hugo_Symbol = sub("_.*","",colnames(RNA_expression)[-1]), 
                             Entrez_Gene_Id = as.numeric(sub(".*_","",colnames(RNA_expression)[-1])),
                             Hugo_Entrez = colnames(RNA_expression)[-1])

meta_gene_GExp$Hugo_Symbol_official <- mapIds(org.Hs.eg.db,
                                              keys= as.character(meta_gene_GExp$Entrez_Gene_Id),
                                              column="SYMBOL",
                                              keytype="ENTREZID",
                                              multiVals="first") %>% as.character()

## re-check on gene name and ID uniqueness
if (any(duplicated(meta_gene_GExp$Hugo_Symbol)) | any(duplicated(meta_gene_GExp$Entrez_Gene_Id))) {warning("CCLE expression data are not ready: maybe duplication in Entrez ID.")}

summary_RNA_expression[summary_RNA_expression$Hugo_Entrez == "RLN2_6013",]
summary_RNA_expression$Hugo_Entrez[summary_RNA_expression$Hugo_Entrez == "RLN2_6013"] <- "RLN2_6019"
summary_RNA_expression$Entrez_Gene_id[summary_RNA_expression$Hugo_Entrez == "RLN2_6013"] <- 6019

table(summary_RNA_expression$Hugo_Entrez %in% meta_gene_GExp$Hugo_Entrez); table(meta_gene_GExp$Hugo_Entrez %in% summary_RNA_expression$Hugo_Entrez)
summary_RNA_expression <- summary_RNA_expression[summary_RNA_expression$Hugo_Entrez %in% meta_gene_GExp$Hugo_Entrez,]



##------- CCLE --------
## RNA_expression should be pre-processed in mother script (MachineLearning.R)
if (any(duplicated(sub("_.*","",colnames(RNA_expression)[-1]))) | any(duplicated(sub(".*_","",colnames(RNA_expression)[-1])))) {warning("CCLE expression data are not ready: maybe duplication in Entrez ID.")}
DataTable <- DT.generator(RNA_expression, DEPMAPID_DT, c("DEPMAPID","PRIMARY_SITE","DOUBLING_TIME"))
## further processing
rownames(DataTable) <- DataTable$DEPMAPID
DataTable <- filter(DataTable, DOUBLING_TIME <400) ## excluding unusual samples

## create a working version of meta_gene_GExp
table(colnames(RNA_expression)[-1] %in% meta_gene_GExp$Hugo_Entrez)
genelist_CCLE_Exp <- meta_gene_GExp; dim(genelist_CCLE_Exp)
# genelist_CCLE_Exp$Hugo_Entrez <- paste(genelist_CCLE_Exp$Hugo_Symbol,genelist_CCLE_Exp$Entrez_Gene_Id, sep = "_") ## for old method
if(any(duplicated(genelist_CCLE_Exp$Hugo_Entrez)) |
   any(duplicated(genelist_CCLE_Exp$Hugo_Symbol)) |
   any(duplicated(genelist_CCLE_Exp$Entrez_Gene_Id))
) {warning("Check the uniqueness of the gene name or ID")}

##------- compatible between CCLE and TCGA -------
## select the missing genes from genelist_CCLE_Exp that is mssing in TCGA (can't be found by name or ID)
## use processed TCGA data
working_table_TCGA <- TCGA.PCA18.cbp_exp.rsem
missing <- missing.gene.compatible(genelist_CCLE_Exp, working_table_TCGA); dim(missing) ## genes can't be found in table2 in either name or ID
genelist_CCLE_Exp <- genelist_CCLE_Exp[!(genelist_CCLE_Exp$Hugo_Entrez %in% missing$Hugo_Entrez),] ;View(genelist_CCLE_Exp)  ## excluding all genes from the missing list
table((duplicated(genelist_CCLE_Exp$Hugo_Symbol)))
table((duplicated(genelist_CCLE_Exp$Entrez_Gene_Id)))
## at this point genelist_CCLE_Exp is a list of genes based on CCLE, and it can be fully mapped into TCGA
## gene names and IDs mapping list between CCLE and TCGA
genelist_CCLE_X_TCGA <- gene.identif.mapping(genelist_CCLE_Exp, working_table_TCGA, ".CCLE",".TCGA") #;View(genelist_CCLE_X_TCGA)
genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA <- paste(genelist_CCLE_X_TCGA$Hugo_Symbol.TCGA, genelist_CCLE_X_TCGA$Entrez_Gene_Id.TCGA, sep = "_")

## CHECK POINT
if (all(genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA %in% working_table_TCGA$Hugo_Entrez)) {cat("All genes from genelist can be found in TCGA.")} else {warning("Genes from genelist are missing in TCGA")}

## both Hugo symbol and Entrez ID are matched
table(genelist_CCLE_X_TCGA$Hugo_Entrez == genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA)
## only Entrez ID are matched
table(genelist_CCLE_X_TCGA$Entrez_Gene_Id.CCLE == genelist_CCLE_X_TCGA$Entrez_Gene_Id.TCGA)

table((genelist_CCLE_X_TCGA$Entrez_Gene_Id.CCLE == genelist_CCLE_X_TCGA$Entrez_Gene_Id.TCGA)+
        (genelist_CCLE_X_TCGA$Hugo_Symbol.CCLE == genelist_CCLE_X_TCGA$Hugo_Symbol.TCGA))

genelist_CCLE_X_TCGA_symbolORid <- genelist_CCLE_X_TCGA
genelist_CCLE_X_TCGA_symbolANDid <- genelist_CCLE_X_TCGA[genelist_CCLE_X_TCGA$Hugo_Entrez == genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA,]

## Pick your CCLE-TCGA mapping strategy to work with
cat("Working with genes mapped into TCGA via name or ID")
genelist_CCLE_X_TCGA <- genelist_CCLE_X_TCGA_symbolORid
cat("Working with genes mapped into TCGA via name and ID simutaneously.")
genelist_CCLE_X_TCGA <- genelist_CCLE_X_TCGA_symbolANDid

## genelist storage
# genelist_CCLE_Exp_collection <- list()
# genelist_CCLE_X_TCGA_collection <- list()
current_working_object; current_working_object_abb
genelist_CCLE_Exp_collection$LUAD <- genelist_CCLE_Exp
genelist_CCLE_X_TCGA_collection$LUAD <- genelist_CCLE_X_TCGA

current_working_object; current_working_object_abb
genelist_CCLE_Exp_collection$LUSC <- genelist_CCLE_Exp
genelist_CCLE_X_TCGA_collection$LUSC <- genelist_CCLE_X_TCGA

current_working_object; current_working_object_abb
genelist_CCLE_Exp_collection$BRCA <- genelist_CCLE_Exp
genelist_CCLE_X_TCGA_collection$BRCA <- genelist_CCLE_X_TCGA




##------- Therapy response data (Freeman, 4 cohorts) pretreatment -----------
## gene annotation
## build connection
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_list_exp <- CPB4_exp_pretr[,1] ## extract gene list
gene_list_exp$ENSEMBL_id <- sub("\\..*","",gene_list_exp$Name) ## remove version number

## via biomaRt package method
## build connection to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl ## current status
## all the filters
ensembl_filters <- listFilters(ensembl)
searchFilters(mart = ensembl, pattern = "gene")
## all the attributes
ensembl_attributes <- listAttributes(ensembl)
searchAttributes(mart = ensembl, pattern = "hgnc")
## search by Ensembl id only
output_id <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_gene_id_version","entrezgene_id"),
                   filters = 'ensembl_gene_id',
                   values = sub("\\..*","",CPB4_exp_pretr$Name),
                   mart = ensembl)
## join the output into gene list
gene_list_exp <- gene_list_exp %>% left_join(output_id, by = c(ENSEMBL_id = "ensembl_gene_id"))

## overview on missing gene
ind <- which(!(summary_RNA_expression$Entrez_Gene_id[1:100] %in% gene_list_exp$entrezgene_id))
ind <- which(!(summary_RNA_expression$Entrez_Gene_id[1:300] %in% gene_list_exp$entrezgene_id))
summary_RNA_expression[ind,]



##------- Therapy response data (Lee, Immunotherapy, targeted drug therapy) pretreatment -----------
## extract list of genes that are overlapping between CCLE and each cohorts, and list of samples
## for immuno, list of overlapping genes
list_gene_immuno_13 <- c()
for(i in 1:length(Immuno_13_opt)){
  ## extract table
  current_table <- Immuno_13_opt[[i]]$mRNA
  ## filter missing values
  ind1 <- apply(Immuno_13_opt[[i]]$mRNA,1,function(X) !any(is.na(X)))  ## check NA
  ind2 <- apply(Immuno_13_opt[[i]]$mRNA,1,function(X) !any(X %in% c("NaN",-Inf)))  ## check character "NaN"
  ind <- (ind1+ind2)==2
  current_list_gene <- rownames(Immuno_13_opt[[i]]$mRNA)[ind]  ## pick non-NA for all samples
  
  ## tailored correlation gene list
  list_gene_immuno_13[[i]] <- summary_RNA_expression[which(summary_RNA_expression$Gene_Transcript %in% current_list_gene),c("Gene_Transcript","Entrez_Gene_id","Hugo_Entrez","order_SPq")]
  names(list_gene_immuno_13)[i] <- names(Immuno_13_opt)[i]
  rm(ind1,ind2,ind,current_table,current_list_gene)
} ## data with available gene # less than 10K should be excluded
## output
list_gene_immuno_13_DT <- list_gene_immuno_13; rm(list_gene_immuno_13)


## for immuno, list of patients
working_table_collection <- Immuno_13
list_patient_collection <- list()
for(i in 1:length(working_table_collection)){
  current_dataset <- names(working_table_collection)[i]
  list_patient <- data.frame(Patient_id = working_table_collection[[current_dataset]]$samples,
                             Response = NA)
  if("surv.dt" %in% names(working_table_collection[[current_dataset]])){
    list_patient$surv_status <- working_table_collection[[current_dataset]]$surv.dt$status
    list_patient$surv_time <- working_table_collection[[current_dataset]]$surv.dt$time
  }
  if("treatment" %in% names(working_table_collection[[current_dataset]])){
    list_patient$Treatment <- working_table_collection[[current_dataset]]$treatment}
  if(current_dataset %in% c("Liu","Riaz")){
    list_patient$Response <- working_table_collection[[current_dataset]]$response
    list_patient_collection[[current_dataset]] <- list_patient; next}
  
  list_patient$Response[working_table_collection[[current_dataset]]$indR] <- "Responder"
  list_patient$Response[working_table_collection[[current_dataset]]$indNR] <- "Non-responder"
  list_patient_collection[[current_dataset]] <- list_patient
}
list_patient_immuno_collection <- list_patient_collection
list_patient_immuno_collection$Miao$Patient_id2 <- Immuno_13_opt$Miao$samples
list_patient_immuno_collection$Miao$Patient_id <- list_patient_immuno_collection$Miao$Patient_id2
list_patient_immuno_collection$Miao <- list_patient_immuno_collection$Miao[,1:2]
list_patient_immuno_collection$Liu <- list_patient_immuno_collection$Liu %>% left_join(Immuno_13_opt$Liu$clin[,c("samples","Primary_Type")], by=c(Patient_id="samples"))


## for targeted drug, list of overlapping genes
list_gene_targeteddrug_10 <- c()
for(i in 1:length(Targeted_10)){
  ## extract table
  current_table <- Targeted_10[[i]]$mRNA
  ## filter missing values
  ind1 <- apply(Targeted_10[[i]]$mRNA,1,function(X) !any(is.na(X)))  ## check NA
  ind2 <- apply(Targeted_10[[i]]$mRNA,1,function(X) !any(X == "NaN"))  ## check character "NaN"
  ind <- (ind1+ind2)==2
  current_list_gene <- rownames(Targeted_10[[i]]$mRNA)[ind]  ## pick non-NA for all samples
  
  ## tailored correlation gene list
  list_gene_targeteddrug_10[[i]] <- summary_RNA_expression[which(summary_RNA_expression$Gene_Transcript %in% current_list_gene),c("Gene_Transcript","Entrez_Gene_id","Hugo_Entrez","order_SPq")]
  names(list_gene_targeteddrug_10)[i] <- names(Targeted_10)[i]
  rm(ind1,ind2,ind,current_table,current_list_gene)
} ## data with available gene # less than 10K should be excluded
## output
list_gene_targeteddrug_10_DT <- list_gene_targeteddrug_10; rm(list_gene_targeteddrug_10)


## for targeted drug, list of patients
working_table_collection <- Targeted_10
list_patient_collection <- list()
for(i in 1:length(working_table_collection)){
  current_dataset <- names(working_table_collection)[i]
  list_patient <- data.frame(Patient_id = working_table_collection[[current_dataset]]$samples,
                             Response = NA)
  if("surv.dt" %in% names(working_table_collection[[current_dataset]])){
    list_patient$surv_status <- working_table_collection[[current_dataset]]$surv.dt$status
    list_patient$surv_time <- working_table_collection[[current_dataset]]$surv.dt$time
  }
  if("treatment" %in% names(working_table_collection[[current_dataset]])){
    list_patient$Treatment <- working_table_collection[[current_dataset]]$treatment}
  if(current_dataset=="GSE32603"){list_patient$Response <- t(Targeted_10[[current_dataset]]$response)
  ind <- list_patient$Response==1
  indd <- list_patient$Response==0
  list_patient$Response[ind] <- "Responder"
  list_patient$Response[indd] <- "Non-responder"
  list_patient_collection[[current_dataset]] <- list_patient; rm(ind,indd); next}
  
  list_patient$Response[working_table_collection[[current_dataset]]$indR] <- "Responder"
  list_patient$Response[working_table_collection[[current_dataset]]$indNR] <- "Non-responder"
  list_patient_collection[[current_dataset]] <- list_patient
}
## output
list_patient_targeteddrug_collection <- list_patient_collection
list_patient_targeteddrug_collection$GSE109211 <- list_patient_targeteddrug_collection$GSE109211 %>% 
  left_join(Targeted_10$GSE109211$clin[,c("geo_accession","treatment:ch1")], by=c(Patient_id="geo_accession"))

##----------- Therapy response data (CTR database) pretreatment --------------
## no tailored model for CTR RNA-seq data, so no need to create a gene list for RNA-seq

## for Microarray:
## gene list
gene_list_MC20161 <- data.frame(X=MC_data$CTR_Microarray_1$X)
gene_list_MC12399 <- data.frame(X=MC_data$CTR_Microarray_10$X)
for(i in 1:length(MC_data)){
  current_dataset <- names(MC_data)[i]
  current_table <- MC_data[[i]]
  if(nrow(current_table)==20161){  ## save it according to number of genes
    gene_list_MC20161$Y <- gene_list_MC20161$X %in% current_table$X
    colnames(gene_list_MC20161)[ncol(gene_list_MC20161)] <- current_dataset
  } else if(nrow(current_table)==12399){
    gene_list_MC12399$Y <- gene_list_MC12399$X %in% current_table$X
    colnames(gene_list_MC12399)[ncol(gene_list_MC12399)] <- current_dataset
  }
}


## for RNA-seq:
## log2(TPM+1) normalization
output_data_tpm <- list()
for(i in 1:length(output_data)){
  working_table <- output_data[[i]]
  if(i==1){
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    gene_info <- getBM(attributes = c("hgnc_symbol", "start_position","end_position"),
                       filters = 'hgnc_symbol',
                       values = working_table$X,
                       mart = ensembl)
    gene_info$gene_length <- gene_info$end_position - gene_info$start_position
  } else if (any(!(output_data[[i]]$X %in% output_data[[1]]$X))){
    gene_info <- getBM(attributes = c("hgnc_symbol", "start_position","end_position"),
                       filters = 'hgnc_symbol',
                       values = working_table$X,
                       mart = ensembl)
    gene_info$gene_length <- gene_info$end_position - gene_info$start_position
  }
  ## delete everything unwanted
  # rm(ensembl,ensembl_ds,ensembl_attributes,ensembl_filters)
  rownames(working_table) <- working_table$X
  working_table$gene_length <- gene_info$gene_length[match(rownames(working_table),gene_info$hgnc_symbol)]
  
  if(ncol(working_table)==3){
    colname <- colnames(working_table)[2]
    rowname <- rownames(working_table)
    # working_table <- data.frame(working_table[,2]/(working_table$gene_length/1000))
    working_table <- working_table[,2]/(working_table$gene_length/1000)
    # working_table <- working_table[apply(working_table, 1, function(X) !any(is.na(X))),]
    working_table <- working_table/sum(as.numeric(unlist(na.omit(working_table))))*10^6 #%>% as.data.frame()
    working_table <- data.frame(working_table)
    colnames(working_table) <- colname; rm(colname)
    rownames(working_table) <- rowname; rm(rowname)
  } else {
    working_table <- apply(working_table[,c(-1,-ncol(working_table))], 2,function(X) X/(working_table$gene_length/1000))
    # working_table <- working_table[apply(working_table, 1, function(X) !any(is.na(X))),]
    working_table <- apply(working_table, 2, function(X) X/sum(as.numeric(na.omit(X)))*10^6) %>% as.data.frame()
  }
  
  working_table <- log2(working_table+1)
  
  output_data_tpm[[i]] <- working_table
  names(output_data_tpm)[i] <- names(output_data)[i]
}; rm(i)
RNAseq_data_tpm <- output_data_tpm


## merge seq data into single data frame
## use the same code for all cases - pick one of the followings:
## for RNA-seq:
data_group <- output_data_tpm
output_data_tpm_df <- output_data_tpm[[1]]
## for microarray:
## MC20161
data_group <- MC_data[colnames(gene_list_MC20161)[-1]]
output_data_tpm_df <- data_group[[1]]
## MC12399
data_group <- MC_data[colnames(gene_list_MC12399)[-1]]
output_data_tpm_df <- data_group[[1]]

## then run this:
for(i in 2:length(data_group)){
  if(any(!(rownames(data_group[[i]])==rownames(output_data_tpm_df)))){stop()}  ## make sure the genes are consistent
  # output_data_tpm_df <- cbind(output_data_tpm_df,data.frame(data_group[[i]][,-1]))
  if(ncol(data_group[[i]])==2){ ## for single sample dataset: has trouble to use cbind()
    output_data_tpm_df <- cbind(output_data_tpm_df,data.frame(Z=data_group[[i]][,-1])); colnames(output_data_tpm_df)[colnames(output_data_tpm_df)=="Z"]<-colnames(data_group[[i]])[2]
  } else {
    output_data_tpm_df <- cbind(output_data_tpm_df,data.frame(data_group[[i]][,-1]))
  }
}
## output - pick one of the following:
## RNA-seq
RNAseq_data_tpm_df <- output_data_tpm_df; rm(output_data_tpm_df)
working_table = RNAseq_data_tpm_df
## microarray
MC20161_data_tpm_df <- output_data_tpm_df; rm(output_data_tpm_df)
working_table = MC20161_data_tpm_df
## microarray
MC12399_data_tpm_df <- output_data_tpm_df; rm(output_data_tpm_df)
colnames(MC12399_data_tpm_df)[grep("CyR",colnames(MC12399_data_tpm_df))] <- gsub("\\.","-",colnames(MC12399_data_tpm_df)[grep("CyR",colnames(MC12399_data_tpm_df))]) ## reverse "check.name"
working_table = MC12399_data_tpm_df



##------- data pre-processing and splitting --------
#### input main working table, DataTable
DataTable <- DT.generator(RNA_expression, DEPMAPID_DT, c("DEPMAPID","PRIMARY_SITE","DOUBLING_TIME")); current_DT_GR <- "DT"
MyCorrList <- summary_RNA_expression

## further processing
rownames(DataTable) <- DataTable$DEPMAPID
if (current_DT_GR == "DT") {
  DataTable <- filter(DataTable, DOUBLING_TIME <400)} ## excluding unusual samples

##### splitting index 
## (splitting on samples)
set.seed(2022)
ind_split <- createDataPartition(y = DataTable$DOUBLING_TIME, p= 0.7)[[1]]  ## splitting index
## train_mother and test_mother have all the predictors/genes, but the samples have been splitted
train_mother <- DataTable[ind_split,-1:-2]
test_mother <- DataTable[-ind_split,-1:-2]

## check sample with NA
ind <- which(!(apply(train_mother, 1,function(X) all(!is.na(X)))))
ind <- which(!(apply(test_mother, 1,function(X) all(!is.na(X)))))
## discard NA sample
test_mother <- test_mother[-ind,]; rm(ind)

## remove genes that was excluded from common gene list
table(MyCorrList$Hugo_Entrez %in% genelist_CCLE_X_TCGA$Hugo_Entrez)
MyCorrList <- MyCorrList[MyCorrList$Hugo_Entrez %in% genelist_CCLE_X_TCGA$Hugo_Entrez,]

## (splitting on genes/predictors via column names)
## data preprocessing is also done in the meanwhile
massive <- c(2,5,10,15,20,30,40,50,70,100,200,300)
j <- 1; train_dt <- list(); test_dt <- list()
gene_mut <- c() ## when mutation is not applied
for (i in massive) {
  ind_col <- col.selector(MyCorrList, "pvalue_spearman", 1:i, c("DOUBLING_TIME"), colnum_id = "Hugo_Entrez")
  
  ## the train_dt and test_dt have the predictors selected, while the samples has been splitted before
  train_dt[[j]] <- dplyr::select(train_mother, all_of(ind_col))
  test_dt[[j]] <- dplyr::select(test_mother, all_of(ind_col))
  
  ## scale (make sure label variables are before numeric variables)
  processCenter <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = c("scale", "center"))
  train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(processCenter,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  ## for regular test_dt splitted from the same DataTable as train_dt
  processCenter <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = c("scale", "center"))
  test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(processCenter,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  #### near-zero variation
  # remove near zero variation for the columns at least
  # 85% of the values are the same
  # this function creates the filter but doesn't apply it yet
  nzv <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))],method="nzv",uniqueCut = 15)
  train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(nzv,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  nzv <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))],method="nzv",uniqueCut = 15)
  test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(nzv,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  ##### correlated predictors
  # create a filter for removing higly correlated variables
  # if two variables are highly correlated only one of them is removed
  corrFilt <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = "corr", cutoff = 0.9)
  train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(corrFilt,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  corrFilt <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = "corr", cutoff = 0.9)
  test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(corrFilt,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  names(train_dt)[j] <- i; names(test_dt)[j] <- i
  j <- j+1
}; rm(i,j,ind_col, processCenter)


# trctrl <- trainControl(method = "none")
trctrl <- trainControl(method='repeatedcv', 
                       number=10, 
                       repeats=6,
                       search = "grid")



##------- Model Training, Random Forest ---------
## massive series training
model_ExpDT <- list(); predicitons_in_ExpDT <- list(); predictions_ExpDT <- list(); j <- 1
for (i in massive) {
  model_ExpDT[[j]] <- train(DOUBLING_TIME ~.,
                            data = train_dt[[j]],
                            method = "ranger",  ## random forest
                            trControl = trctrl
  )
  names(model_ExpDT)[j] <- paste("model_ExpDT",i,sep = "_")
  
  ## in-training validation
  predicitons_in_ExpDT[[j]] <- predict(model_ExpDT[[j]], train_dt[[j]][,-1])
  ## test-data validation
  predictions_ExpDT[[j]] <- predict(model_ExpDT[[j]], test_dt[[j]][,-1])
  
  j <- j+1
}; rm(i,j)
names(predictions_ExpDT) <- massive
print(model_ExpDT[[5]])

## put more information into predictions_ExpDT
for(i in 1:length(predictions_ExpDT)) {
  names(predictions_ExpDT[[i]]) <- rownames(test_dt[[1]])  ## marked with patient id
  predictions_ExpDT[[i]] <- data.frame(predicted_DT = predictions_ExpDT[[i]])  ## convert it into data frame
  predictions_ExpDT[[i]]$DEPMAPID <- rownames(predictions_ExpDT[[i]])
  ## further annotation
  predictions_ExpDT[[i]] <- left_join(predictions_ExpDT[[i]], DEPMAPID_DT[,c("DEPMAPID","DOUBLING_TIME","PRIMARY_SITE")], by = "DEPMAPID")
}

## output
## This is the regular DTLearner model trained by 70/30 splitting on CCLE; it is also used for RNA-seq datasets from CTR database
model_ExpDT_reg <- model_ExpDT; rm(model_ExpDT)
predicitons_ExpDT_reg <- predictions_ExpDT; rm(predictions_ExpDT)
saveRDS(model_ExpDT_reg,"model_ExpDT_reg.rds")

## model testing
val_ExpDT <- data.frame(Number_of_predictor = as.numeric(sub(".*_","",names(model_ExpDT))), 
                        model.name = names(model_ExpDT))
for (i in 1:length(model_ExpDT)) {
  ## validation on testing data
  val_ExpDT$t.test[i] <- t.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME)$p.value
  val_ExpDT$Pearson_p[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "pearson")$p.value
  val_ExpDT$Spearman_p[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "spearman")$p.value
  val_ExpDT$Pearson_est[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "pearson")$estimate
  val_ExpDT$Spearman_est[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "spearman")$estimate
  ## validation on training data
  val_ExpDT$t.test_inTrain[i] <- t.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME)$p.value
  val_ExpDT$Pearson_p_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "pearson")$p.value
  val_ExpDT$Spearman_p_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "spearman")$p.value
  val_ExpDT$Pearson_est_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "pearson")$estimate
  val_ExpDT$Spearman_est_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "spearman")$estimate
}

## output
val_ExpDT_reg <- val_ExpDT; rm(val_ExpDT)

### plots between prediction and actual DT on test set cell lines
plot_data <- predicitons_ExpDT_reg$'50'
ptext_loca <- axis_location(plot_data,1,3,"rightup") %>% log10()
input_bquote <- bquote(paste('Spr. p-value ='~ .(sub("e.*","",signif(val_ExpDT_reg$Spearman_p[val_ExpDT_reg$Number_of_predictor==50],digits=4))%>%as.numeric())
                             %*% 10^.(sub(".*e","",signif(val_ExpDT_reg$Spearman_p[val_ExpDT_reg$Number_of_predictor==50],digits=4))%>%as.numeric())))

plot <-
  ggplot(plot_data, aes(log10(predicted_DT),log10(DOUBLING_TIME))) + geom_point(color = "dodgerblue4",size=2,alpha=0.7) + #scale_fill_distiller(palette = "Blues")+
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(hjust = 0.5), legend.position = "none") +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('Predicted Doubling Time log' ['10'] ('hr')), y = bquote('Doubling Time log' ['10'] ('hr')), title = "pDT by Classical Model on Test Set") +
  geom_smooth(method=lm, se=FALSE, col='grey', linewidth = 1, linetype = 2) +
  scale_x_continuous(limits = c(min(log10(plot_data$predicted_DT)),max(log10(plot_data$predicted_DT)*1.05))) +
  stat_regline_equation(label.x = ptext_loca[1]*0.95, label.y = ptext_loca[2]*1.07, size = 6,family = "Times New Roman")+ ## trendline equation
  stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
           label.x = ptext_loca[1]*0.97, label.y = ptext_loca[2]*1.04, size = 6,family = "Times New Roman")+ ## R^2
  annotate("text",family = "Times New Roman", x=ptext_loca[1], y=ptext_loca[2], size = 6,
           label = input_bquote)+  ## pvalue
  annotate("text",family = "Times New Roman", x=ptext_loca[1], y=ptext_loca[2]*0.96, size = 6,
           label = paste("Rho =",signif(val_ExpDT_reg$Spearman_est[val_ExpDT_reg$Number_of_predictor==50],digits=4)))  ## est.

## This is Figure 3A ##
ggsave(filename = "plot/classical_model_50_DT_test.pdf",
       plot,
       width = 8.1, height = 6
)


### MKI67 for DT
plot_data <- RNA_expression[,c(1,which(sub("_.*","",colnames(RNA_expression)) %in% "MKI67"))] %>% 
  left_join(DEPMAPID_DT[,c("DEPMAPID","DOUBLING_TIME","PRIMARY_SITE")]) %>%
  filter(!is.na(DOUBLING_TIME) & DEPMAPID %in% predicitons_ExpDT_reg$'50'$DEPMAPID)
ptext_loca <- axis_location(plot_data,2,3,"leftup"); ptext_loca[2] <- log10(ptext_loca[2])

plot_2 <-
  ggplot(plot_data, aes(MKI67_4288,log10(DOUBLING_TIME)))  + geom_point(color = "dodgerblue4",size=2,alpha=0.7) + #scale_colour_gradient(trans="reverse")+
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(hjust = 0.5), legend.position="none") +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x="MKI67 Expression Level",y = bquote('Doubling Time log' ['10'] ('hr')), title = "MKI67 and Cell Line DT on Test Set") +
  geom_smooth(method=lm, se=FALSE, col='grey', linewidth = 1, linetype = 2) +
  stat_regline_equation(label.x = ptext_loca[1]*0.95, label.y = ptext_loca[2]*1.06, size = 6,family = "Times New Roman")+ ## trendline equation
  stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
           label.x = ptext_loca[1]*0.97, label.y = ptext_loca[2]*1.035, size = 6,family = "Times New Roman")+ ## R^2
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*1.1, y=ptext_loca[2], size = 6,
           label = paste("Spr. p-value =",signif(cor.test(plot_data$MKI67_4288,plot_data$DOUBLING_TIME,method="spearman")$p.value, digits = 4)))+  ## pvalue
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*1.1, y=ptext_loca[2]*0.965, size = 6,
           label = paste("Rho =",signif(cor.test(plot_data$MKI67_4288,plot_data$DOUBLING_TIME,method="spearman")$est,digits=4)))  ## est.

## This is Figure 3B ##
ggsave(filename = "plot/MKI67_50_DT_test.pdf",
       plot_2,
       width = 7, height = 6
)


##----- gene list from the dataset for tailored model -------
## in order to have the tailored model, a list of expression data genes from targeted dataset is needed ahead. 
## The list is used to exclude the CCLE expression gene for training that is not in the list 

## for GENT2
list_gene_P2

## for CTR database
gene_list_MC20161 ## CTR-MC20161
gene_list_MC12399 ## CTR-MC12399

## for targeted drug datasets from Lee
list_gene_targeteddrug_10_DT
## for immunotherapy datasets from Lee
list_gene_immuno_13_DT

##------- tailored model for GENT2, and CTR database for therapy response --------
## data pre-processing and splitting
#### input main working table, DataTable
DataTable <- DT.generator(RNA_expression, DEPMAPID_DT, c("DEPMAPID","PRIMARY_SITE","DOUBLING_TIME")); current_DT_GR <- "DT"
MyCorrList <- summary_RNA_expression

## further processing
rownames(DataTable) <- DataTable$DEPMAPID
if (current_DT_GR == "DT") {
  DataTable <- filter(DataTable, DOUBLING_TIME <400)} ## excluding unusual samples

##### splitting index 
## (splitting on samples)
set.seed(2022)
ind_split <- createDataPartition(y = DataTable$DOUBLING_TIME, p= 0.7)[[1]]  ## splitting index
## train_mother and test_mother have all the predictors/genes, but the samples have been splitted
train_mother <- DataTable[ind_split,-1:-2]
test_mother <- DataTable[-ind_split,-1:-2]

## check sample with NA
ind <- which(!(apply(train_mother, 1,function(X) all(!is.na(X)))))
ind <- which(!(apply(test_mother, 1,function(X) all(!is.na(X)))))
## discard NA sample
test_mother <- test_mother[-ind,]; rm(ind)


## remove genes that was excluded from common gene list
table(MyCorrList$Hugo_Entrez %in% genelist_CCLE_X_TCGA$Hugo_Entrez)
MyCorrList <- MyCorrList[MyCorrList$Hugo_Entrez %in% genelist_CCLE_X_TCGA$Hugo_Entrez,]

## if any filter on genes is applied: (e.g. tailored model): pick one of the cases
## for GENT2
MyCorrList <- MyCorrList[MyCorrList$Gene_Transcript %in% list_gene_P2$list_gene_P2,] ## GENT2
## for CTR database
## pick one of the following:
working_list_gene2 <- gene_list_MC20161 ## CTR-MC20161
working_list_gene2 <- gene_list_MC12399 ## CTR-MC12399
## then run this:
colnames(working_list_gene2)[1] <-"Gene_Transcript"
MyCorrList <- MyCorrList[MyCorrList$Gene_Transcript %in% working_list_gene2$Gene_Transcript,]  ## single table (not list)


## (splitting on genes/predictors via column names)
## data preprocessing is also done in the meanwhile
massive <- c(50)
j <- 1; train_dt <- list(); test_dt <- list()
gene_mut <- c() ## when mutation is not applied
for (i in massive) {
  ind_col <- col.selector(MyCorrList, "pvalue_spearman", 1:i, c("DOUBLING_TIME"), colnum_id = "Hugo_Entrez")
  
  ## the train_dt and test_dt have the predictors selected, while the samples has been splitted before
  train_dt[[j]] <- dplyr::select(train_mother, all_of(ind_col))
  test_dt[[j]] <- dplyr::select(test_mother, all_of(ind_col))
  
  ## scale (make sure label variables are before numeric variables)
  processCenter <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = c("scale", "center"))
  train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(processCenter,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  ## for regular test_dt splitted from the same DataTable as train_dt
  processCenter <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = c("scale", "center"))
  test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(processCenter,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  #### near-zero variation
  ## remove near zero variation for the columns at least
  ## 85% of the values are the same
  ## this function creates the filter but doesn't apply it yet
  nzv <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))],method="nzv",uniqueCut = 15)
  train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(nzv,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  nzv <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))],method="nzv",uniqueCut = 15)
  test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(nzv,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  ##### correlated predictors
  ## create a filter for removing higly correlated variables
  ## if two variables are highly correlated only one of them is removed
  corrFilt <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = "corr", cutoff = 0.9)
  train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(corrFilt,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  corrFilt <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = "corr", cutoff = 0.9)
  test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(corrFilt,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
  
  names(train_dt)[j] <- i; names(test_dt)[j] <- i
  j <- j+1
}; rm(i,j,ind_col, processCenter)


trctrl <- trainControl(method='repeatedcv', 
                       number=10, 
                       repeats=6,
                       search = "grid")

## Model Training, Random Forest
## massive series training
model_ExpDT <- list(); predicitons_in_ExpDT <- list(); predictions_ExpDT <- list(); j <- 1
for (i in massive) {
  model_ExpDT[[j]] <- train(DOUBLING_TIME ~.,
                            data = train_dt[[j]],
                            method = "ranger",  ## random forest
                            trControl = trctrl
  )
  names(model_ExpDT)[j] <- paste("model_ExpDT",i,sep = "_")
  
  ## in-training validation
  predicitons_in_ExpDT[[j]] <- predict(model_ExpDT[[j]], train_dt[[j]][,-1])
  ## test-data validation
  predictions_ExpDT[[j]] <- predict(model_ExpDT[[j]], test_dt[[j]][,-1])
  
  j <- j+1
}; rm(i,j)
names(predictions_ExpDT) <- massive

## put more information into predictions_ExpDT
for(i in 1:length(predictions_ExpDT)) {
  names(predictions_ExpDT[[i]]) <- rownames(test_dt[[1]])  ## marked with patient id
  predictions_ExpDT[[i]] <- data.frame(predicted_DT = predictions_ExpDT[[i]])  ## convert it into data frame
  predictions_ExpDT[[i]]$DEPMAPID <- rownames(predictions_ExpDT[[i]])
  ## further annotation
  predictions_ExpDT[[i]] <- left_join(predictions_ExpDT[[i]], DEPMAPID_DT[,c("DEPMAPID","DOUBLING_TIME","PRIMARY_SITE")], by = "DEPMAPID")
}

## model testing
val_ExpDT <- data.frame(Number_of_predictor = as.numeric(sub(".*_","",names(model_ExpDT))), 
                        model.name = names(model_ExpDT))
for (i in 1:length(model_ExpDT)) {
  ## validation on testing data
  val_ExpDT$t.test[i] <- t.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME)$p.value
  val_ExpDT$Pearson_p[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "pearson")$p.value
  val_ExpDT$Spearman_p[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "spearman")$p.value
  val_ExpDT$Pearson_est[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "pearson")$estimate
  val_ExpDT$Spearman_est[i] <- cor.test(predictions_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "spearman")$estimate
  ## validation on training data
  val_ExpDT$t.test_inTrain[i] <- t.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME)$p.value
  val_ExpDT$Pearson_p_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "pearson")$p.value
  val_ExpDT$Spearman_p_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "spearman")$p.value
  val_ExpDT$Pearson_est_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "pearson")$estimate
  val_ExpDT$Spearman_est_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "spearman")$estimate
}


## output
## GENT2 tailored model
model_ExpDT_gent2_reg <- model_ExpDT; rm(model_ExpDT)
predicitons_ExpDT_gent2_reg <- predictions_ExpDT; rm(predictions_ExpDT)
val_ExpDT_gent2_reg <- val_ExpDT; rm(val_ExpDT)

## CTR-DB
## MC20161
model_ExpDT_ctr_MC20161 <- model_ExpDT; rm(model_ExpDT)
predicitons_ctr_ExpDT_MC20161 <- predictions_ExpDT; rm(predictions_ExpDT)
val_ExpDT_ctr_MC20161 <- val_ExpDT; rm(val_ExpDT)

## MC12399
model_ExpDT_ctr_MC12399 <- model_ExpDT; rm(model_ExpDT)
predicitons_ctr_ExpDT_MC12399 <- predictions_ExpDT; rm(predictions_ExpDT)
val_ExpDT_ctr_MC12399 <- val_ExpDT; rm(val_ExpDT)


##------- tailored model for targeted drug and immunotherapy for therapy response --------
## this code is for training the tailored models for therapy response: 
## similar to before, except working_list_gene is used to select the CCLE genes that also appear in targeted dataset
#### input main working table, DataTable
DataTable <- DT.generator(RNA_expression, DEPMAPID_DT, c("DEPMAPID","PRIMARY_SITE","DOUBLING_TIME")); current_DT_GR <- "DT"
MyCorrList <- summary_RNA_expression

## further processing
rownames(DataTable) <- DataTable$DEPMAPID
if (current_DT_GR == "DT") {
  DataTable <- filter(DataTable, DOUBLING_TIME <400)} ## excluding unusual samples

##### splitting index 
## (splitting on samples)
set.seed(2022)
ind_split <- createDataPartition(y = DataTable$DOUBLING_TIME, p= 0.7)[[1]]  ## splitting index
## train_mother and test_mother have all the predictors/genes, but the samples have been splitted
train_mother <- DataTable[ind_split,-1:-2]
test_mother <- DataTable[-ind_split,-1:-2]

## check sample with NA
ind <- which(!(apply(train_mother, 1,function(X) all(!is.na(X)))))
ind <- which(!(apply(test_mother, 1,function(X) all(!is.na(X)))))
## discard NA sample
test_mother <- test_mother[-ind,]; rm(ind)

## pick one of the following based on which tailored model is going to be trained
## for targeted drug datasets
working_list_gene <- list_gene_targeteddrug_10_DT
## for immunotherapy datasets
working_list_gene <- list_gene_immuno_13_DT

model_ExpDT_collection <- list(); predicitons_ExpDT_collection <- list(); val_ExpDT_collection <- list()
for(k in 1:length(working_list_gene)){
  name_working_list_gene <- names(working_list_gene)[k]
  
  MyCorrList <- summary_RNA_expression
  
  MyCorrList <- MyCorrList[MyCorrList$Hugo_Entrez %in% genelist_CCLE_X_TCGA$Hugo_Entrez,]
  
  ## filtering missing genes in targeted datasets
  MyCorrList <- MyCorrList[MyCorrList$Gene_Transcript %in% working_list_gene[[i]]$Gene_Transcript,]
  
  ## (splitting on genes/predictors via column names)
  ## data preprocessing is also done in the meanwhile
  massive <- c(50)  ## the pipeline has been settled down with 50 genes (of course, it is okay to do a full trials)
  j <- 1; train_dt <- list(); test_dt <- list()
  gene_mut <- c() ## when mutation is not applied
  for (i in massive) {
    ind_col <- col.selector(MyCorrList, "pvalue_spearman", 1:i, c("DOUBLING_TIME"), colnum_id = "Hugo_Entrez")
    
    ## the train_dt and test_dt have the predictors selected, while the samples has been splitted before
    train_dt[[j]] <- dplyr::select(train_mother, all_of(ind_col))
    test_dt[[j]] <- dplyr::select(test_mother, all_of(ind_col))
    
    ## scale (make sure label variables are before numeric variables)
    processCenter <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = c("scale", "center"))
    train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(processCenter,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
    
    ## for regular test_dt splitted from the same DataTable as train_dt
    processCenter <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = c("scale", "center"))
    test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(processCenter,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
    
    #### near-zero variation
    # remove near zero variation for the columns at least
    # 85% of the values are the same
    # this function creates the filter but doesn't apply it yet
    nzv <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))],method="nzv",uniqueCut = 15)
    train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(nzv,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
    nzv <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))],method="nzv",uniqueCut = 15)
    test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(nzv,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
    
    ##### correlated predictors
    # create a filter for removing higly correlated variables
    # if two variables are highly correlated only one of them is removed
    corrFilt <- preProcess(train_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = "corr", cutoff = 0.9)
    train_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(corrFilt,train_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
    corrFilt <- preProcess(test_dt[[j]][,-1:(-1*(length(gene_mut)+1))], method = "corr", cutoff = 0.9)
    test_dt[[j]][,-1:(-1*(length(gene_mut)+1))] <- predict(corrFilt,test_dt[[j]][,-1:(-1*(length(gene_mut)+1))])
    
    names(train_dt)[j] <- i; names(test_dt)[j] <- i
    j <- j+1
  }; rm(i,j,ind_col, processCenter)
  
  trctrl <- trainControl(method='repeatedcv',
                         number=10,
                         repeats=6,
                         search = "grid")
  
  ## Model Training
  ## massive series training
  model_ExpDT <- list(); predicitons_in_ExpDT <- list(); predicitons_ExpDT <- list(); j <- 1
  for (i in massive) {
    model_ExpDT[[j]] <- train(DOUBLING_TIME ~.,
                              data = train_dt[[j]],
                              method = "ranger",  ## random forest
                              trControl = trctrl
    )
    names(model_ExpDT)[j] <- paste("model_ExpDT",i,sep = "_")
    
    ## in-training validation
    predicitons_in_ExpDT[[j]] <- predict(model_ExpDT[[j]], train_dt[[j]][,-1])
    ## test-data validation
    ## for regular test_dt (train and test are splitting from same DataTable)
    predicitons_ExpDT[[j]] <- predict(model_ExpDT[[j]], test_dt[[j]][,-1])
    
    j <- j+1
  }; rm(i,j)
  names(predicitons_ExpDT) <- massive
  
  ## put more information into predicitons_ExpDT
  for(i in 1:length(predicitons_ExpDT)) {
    names(predicitons_ExpDT[[i]]) <- rownames(test_dt[[1]])  ## marked with patient id
    predicitons_ExpDT[[i]] <- data.frame(predicted_DT = predicitons_ExpDT[[i]])  ## convert it into data frame
    predicitons_ExpDT[[i]]$DEPMAPID <- rownames(predicitons_ExpDT[[i]])
    ## further annotation
    predicitons_ExpDT[[i]] <- left_join(predicitons_ExpDT[[i]], DEPMAPID_DT[,c("DEPMAPID","DOUBLING_TIME","PRIMARY_SITE")], by = "DEPMAPID")
  }; rm(i)
  
  val_ExpDT <- data.frame(Number_of_predictor = as.numeric(sub(".*_","",names(model_ExpDT))), 
                          model.name = names(model_ExpDT))
  for (i in 1:length(model_ExpDT)) {
    ## validation on testing data
    val_ExpDT$t.test[i] <- t.test(predicitons_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME)$p.value
    val_ExpDT$Pearson_p[i] <- cor.test(predicitons_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "pearson")$p.value
    val_ExpDT$Spearman_p[i] <- cor.test(predicitons_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "spearman")$p.value
    val_ExpDT$Pearson_est[i] <- cor.test(predicitons_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "pearson")$estimate
    val_ExpDT$Spearman_est[i] <- cor.test(predicitons_ExpDT[[i]]$predicted_DT, test_dt[[i]]$DOUBLING_TIME, method = "spearman")$estimate
    ## validation on training data
    val_ExpDT$t.test_inTrain[i] <- t.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME)$p.value
    val_ExpDT$Pearson_p_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "pearson")$p.value
    val_ExpDT$Spearman_p_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "spearman")$p.value
    val_ExpDT$Pearson_est_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "pearson")$estimate
    val_ExpDT$Spearman_est_inTrain[i] <- cor.test(predicitons_in_ExpDT[[i]], train_dt[[i]]$DOUBLING_TIME, method = "spearman")$estimate
  }; rm(i)
  
  ## output
  ## targeted drug tailored model
  model_ExpDT_collection[[name_working_list_gene]] <- model_ExpDT; rm(model_ExpDT)
  predicitons_ExpDT_collection[[name_working_list_gene]] <- predicitons_ExpDT; rm(predicitons_ExpDT)
  val_ExpDT_collection[[name_working_list_gene]] <- val_ExpDT; rm(val_ExpDT)
}

## output
## targeted drug by DT: one model for one cohorts (missing genes vary)
model_ExpDT_targeteddrug_collection <- model_ExpDT_collection; rm(model_ExpDT_collection)
predicitons_targeteddrug_ExpDT_collection <- predicitons_ExpDT_collection; rm(predicitons_ExpDT_collection)
val_ExpDT_targeteddrug_collection <- val_ExpDT_collection; rm(val_ExpDT_collection)

## immuno: one model for one cohorts (missing genes vary)
model_ExpDT_immuno_collection <- model_ExpDT_collection; rm(model_ExpDT_collection)
predicitons_immuno_ExpDT_collection <- predicitons_ExpDT_collection; rm(predicitons_ExpDT_collection)
val_ExpDT_immuno_collection <- val_ExpDT_collection; rm(val_ExpDT_collection)



##------- prediction on TCGA-test data (by loop) ---------
TCGAdata_list <- names(TCGA.PCA18.cbp_exp.rsem)
TCGAdata_list <- paste(TCGAdata_list,"TCGA",sep = "_")
# which(TCGAdata_list == "LAML_TCGA")
# TCGAdata_list <- TCGAdata_list[-6]
MyCorrList <- summary_RNA_expression
MyCorrList <- summary_DataTable_fullCCLE

## DT
model_ExpDT <- model_ExpDT_reg; current_DT_GR <- "DT"
MyCorrList <- summary_RNA_expression

## if lineage is involved
lineage_list_CCLE <- unique(DEPMAPID_DT$PRIMARY_SITE)
unique(CCLE_sample$PRIMARY_SITE)
unique(CCLE_sample$SITE_SUBTYPE1)
lineage_list_CCLE <- lineage_list_CCLE[!(lineage_list_CCLE %in% c("Na",NA))]  ## excluding missing value
lineage_list_TCGA <- data.frame(Abbreviation = c("LUAD_TCGA","LUSC_TCGA","BRCA_TCGA","DLBCell_TCGA","GBM_TCGA",
                                                 "LAML_TCGA","SKCM_TCGA","ACC_TCGA","CHOL_TCGA","BLCA_TCGA","COAD_TCGA",
                                                 "LGG_TCGA","CESC_TCGA","ESCA_TCGA","STAD_TCGA","UVM_TCGA","HNSC_TCGA",
                                                 "KIRP_TCGA","KICH_TCGA","KIRC_TCGA","LIHC_TCGA","OV_TCGA","PAAD_TCGA",
                                                 "MESO_TCGA","PRAD_TCGA","PCPG_TCGA","SARC_TCGA","TGCT_TCGA","THYM_TCGA",
                                                 "THCA_TCGA","UCEC_TCGA","UCS_TCGA"),
                                Name = c("Lung adenocarcinoma","Lung squamous cell carcinoma","Breast invasive carcinoma",
                                         "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma","Glioblastoma multiforme",
                                         "Acute Myeloid Leukemia","Skin Cutaneous Melanoma","Adrenocortical carcinoma","Cholangiocarcinoma",
                                         "Bladder Urothelial Carcinoma","Colon adenocarcinoma","Brain Lower Grade Glioma",
                                         "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
                                         "Esophageal carcinoma","Stomach adenocarcinoma","Uveal Melanoma",
                                         "Head and Neck squamous cell carcinoma","Kidney renal papillary cell carcinoma",
                                         "Kidney Chromophobe","Kidney renal clear cell carcinoma","Liver hepatocellular carcinoma",
                                         "Ovarian serous cystadenocarcinoma","Pancreatic adenocarcinoma","Mesothelioma",
                                         "Prostate adenocarcinoma","Pheochromocytoma and Paraganglioma","Sarcoma","Testicular Germ Cell Tumors",
                                         "Thymoma","Thyroid carcinoma","Uterine Corpus Endometrial Carcinoma","Uterine Carcinosarcoma"),
                                PRIMARY_SITE.CCLE = c("Lung","Lung","Breast","Haematopoietic_And_Lymphoid_Tissue","Central_Nervous_System",
                                                      "Haematopoietic_And_Lymphoid_Tissue", "Skin","Kidney","Biliary_Tract","Urinary_Tract",
                                                      "Large_Intestine","Central_Nervous_System","Endometrium","Oesophagus","Stomach",
                                                      "Eye", "Upper_Aerodigestive_Tract","Kidney","Kidney","Kidney","Liver","Ovary","Pancreas",
                                                      "Pleura","Prostate","Soft_Tissue","Soft_Tissue","Testis","Thymus","Thyroid",
                                                      "Endometrium","Endometrium"))
## those that are missing
lineage_list_TCGA$PRIMARY_SITE.CCLE[!(lineage_list_TCGA$PRIMARY_SITE.CCLE %in% lineage_list_CCLE)]
lineage_list_CCLE[!(lineage_list_CCLE %in% lineage_list_TCGA$PRIMARY_SITE.CCLE)]


working_table_TCGA_X_ExpSpTop_collection_temp <- list()
working_table_TCGA_X_ExpSpTop_sub_collection_temp <- list()
working_table_TCGA_X_ExpSpTop_predict_collection_temp <- list()
TCGA.PCA18_predict_score_collection_temp <- list()
plot_working_table_TCGA_X_ExpSpTop_predict_temp <- list()

surv_test_collection_temp <- list()

## survival test score container
template <- data.frame(matrix(NA, length(TCGAdata_list), 6))
colnames(template) <- c("cancer_type","KM_pvalue","KM_pvalue_neglog10","cox_wald_pvalue", "cox_wald_pvalue_neglog10", "cox_hazard_ratio")
template$cancer_type <- TCGAdata_list
summary_surv_score_temp <- list(template,template,template,template,template,
                                template,template,template,template,template,
                                template,template); rm(template)
names(summary_surv_score_temp) <- names(model_ExpDT_reg)
for (k in 1:length(TCGAdata_list)) {
  working_table_TCGA <- TCGA.PCA18.cbp_exp.rsem[[k]]; working_table_TCGA_patient <- TCGA.PCA18.cbp_patient_RAW[[k]]
  current_working_object_abb <- TCGAdata_list[k]
  
  ## in case that lineage model is used but TCGA tissue type can't be found in CCLE: skip
  if ("PRIMARY_SITE" %in% colnames(model_ExpDT[[1]]$trainingData) & 
      !(lineage_list_TCGA$PRIMARY_SITE.CCLE[lineage_list_TCGA$Abbreviation == current_working_object_abb] %in% lineage_list_CCLE)) {
    predicitons_ExpDT_fullCCLE_lineage_collection[[current_working_object_abb]] <- NA; val_ExpDT_fullCCLE_lineage_collection[[current_working_object_abb]] <- NA
    next
  }
  
  ## object to hold the result temporary - before final output to regular/GR-classified object
  working_table_TCGA_X_ExpSpTop_temp <- list()
  working_table_TCGA_X_ExpSpTop_sub_temp <- list()
  working_table_TCGA_X_ExpSpTop_predict_temp <- list()
  TCGA.PCA18_predict_score_temp <- list()
  plot_temp <- list()
  surv_test_temp <- list()
  ## hold the statistical test result/score
  TCGA.PCA18_predict_score_temp <- list(pDT_OS_spearman_p = as.data.frame(matrix(NA, 1, length(model_ExpDT))), 
                                        pDT_DSS_spearman_p = as.data.frame(matrix(NA, 1, length(model_ExpDT))),
                                        pDT_PFS_spearman_p = as.data.frame(matrix(NA, 1, length(model_ExpDT))),
                                        pDT_OS_spearman_est = as.data.frame(matrix(NA, 1, length(model_ExpDT))),
                                        pDT_DSS_spearman_est = as.data.frame(matrix(NA, 1, length(model_ExpDT))),
                                        pDT_PFS_spearman_est = as.data.frame(matrix(NA, 1, length(model_ExpDT))))
  for(i in 1:length(TCGA.PCA18_predict_score_temp)){
    colnames(TCGA.PCA18_predict_score_temp[[i]]) <- names(model_ExpDT)
    # rownames(TCGA.PCA18_predict_score_temp[[i]])
  }; rm(i)
  
  j<-1 ## index for tracking
  ## from now on, it depends on model/predictor used
  for (model in names(model_ExpDT)) {
    predictorNumber <- as.numeric(sub(".*_","",model))
    
    MyCorrList <- summary_RNA_expression
    
    ## now exclude the genes from the CCLE list that is missing in TCGA (according to gene list)
    MyCorrList <- MyCorrList[MyCorrList$Hugo_Entrez %in% genelist_CCLE_X_TCGA$Hugo_Entrez,]
    
    ## a version of TCGA that compatible to CCLE
    working_table_TCGA_X_ExpSpTop <- working_table_TCGA; #current_working_object
    working_table_TCGA_X_ExpSpTop$Hugo_Symbol <- genelist_CCLE_X_TCGA[match(working_table_TCGA_X_ExpSpTop$Hugo_Entrez ,genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA), "Hugo_Symbol.CCLE"]
    working_table_TCGA_X_ExpSpTop$Entrez_Gene_Id <- genelist_CCLE_X_TCGA[match(working_table_TCGA_X_ExpSpTop$Hugo_Entrez ,genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA), "Entrez_Gene_Id.CCLE"]
    working_table_TCGA_X_ExpSpTop$Hugo_Entrez <- genelist_CCLE_X_TCGA[match(working_table_TCGA_X_ExpSpTop$Hugo_Entrez ,genelist_CCLE_X_TCGA$Hugo_Entrez.TCGA), "Hugo_Entrez"]
    
    working_table_TCGA_X_ExpSpTop <- working_table_TCGA_X_ExpSpTop[working_table_TCGA_X_ExpSpTop$Hugo_Entrez %in% MyCorrList$Hugo_Entrez[order(MyCorrList$pvalue_spearman)[1:predictorNumber]],]
    ## now exclude the genes from the CCLE list that is missing in TCGA
    MyCorrList <- MyCorrList[MyCorrList$Hugo_Entrez %in% working_table_TCGA_X_ExpSpTop$Hugo_Entrez,]
    
    ## double check: the subset contains every genes we want according to symbol and id
    if (!all(working_table_TCGA_X_ExpSpTop$Hugo_Entrez %in% MyCorrList$Hugo_Entrez[order(MyCorrList$pvalue_spearman)[1:predictorNumber]])) {
      stop("Key genes are missing in TCGA data")
    }
    
    ## mark with the correlation result
    working_table_TCGA_X_ExpSpTop <- left_join(working_table_TCGA_X_ExpSpTop,  
                                               MyCorrList[,c("Hugo_Entrez", "pvalue_pearson","pvalue_spearman","order_SPq")], 
                                               by = "Hugo_Entrez") %>%
      relocate(pvalue_pearson, pvalue_spearman, order_SPq, .after = Hugo_Entrez) ## move the annotation columns forward, but after third column
    
    ## remove the cases with NA data on predictor
    ind_NAcase <- apply(working_table_TCGA_X_ExpSpTop, 2, function(X) any(is.na(X)))  ## which cases have at least one NA data
    ## final input data
    working_table_TCGA_X_ExpSpTop_sub <- working_table_TCGA_X_ExpSpTop[,!(ind_NAcase)]  ## keep those without any NA data
    rm(ind_NAcase)
    
    ## choose one of the following. depends on gene name or Hugo_Entrez is involved (remember to also change the one ~20 lines after this one)
    # rownames(working_table_TCGA_X_ExpSpTop_sub) <- working_table_TCGA_X_ExpSpTop$Hugo_Symbol
    rownames(working_table_TCGA_X_ExpSpTop_sub) <- working_table_TCGA_X_ExpSpTop$Hugo_Entrez
    if (length(base::grep("TCGA", colnames(working_table_TCGA_X_ExpSpTop_sub))) >= 1) {  ## at least 1 available cases: continue
      working_table_TCGA_X_ExpSpTop_sub <- dplyr::select(working_table_TCGA_X_ExpSpTop_sub, 
                                                         base::grep("TCGA", colnames(working_table_TCGA_X_ExpSpTop_sub))) %>%  ## give up the none-value part
        t()
    } else { ## no available cases: ends ahead and start next TCGA data
      # working_table_TCGA_X_ExpSpTop_sub
      ## testing output
      for(i in 1:length(TCGA.PCA18_predict_score_temp)) {TCGA.PCA18_predict_score_temp[[i]][1,model] <- NA}
      ## data output
      working_table_TCGA_X_ExpSpTop_temp[[model]] <- working_table_TCGA_X_ExpSpTop
      working_table_TCGA_X_ExpSpTop_sub_temp[[model]] <- working_table_TCGA_X_ExpSpTop_sub
      working_table_TCGA_X_ExpSpTop_predict_temp[[model]] <- NA
      plot_temp[[model]] <- NA
      j <- j+1
      next
    }
    
    ## choose one of the following. depends on gene name or Hugo_Entrez is involved
    # working_table_TCGA_X_ExpSpTop_sub <- working_table_TCGA_X_ExpSpTop_sub[,MyCorrList$Gene_Transcript[order(MyCorrList$pvalue_spearman)[1:predictorNumber]]]  ## reorder the columns/predictors
    working_table_TCGA_X_ExpSpTop_sub <- working_table_TCGA_X_ExpSpTop_sub[,MyCorrList$Hugo_Entrez[order(MyCorrList$pvalue_spearman)[1:predictorNumber]]]  ## reorder the columns/predictors
    
    ## CHECK POINT
    # if (!all(MyCorrList$Gene_Transcript[order(MyCorrList$pvalue_spearman)[1:predictorNumber]] == colnames(working_table_TCGA_X_ExpSpTop_sub))) {stop("TCGA gene names and order are not compatible to correlation list")}
    # if (!all(MyCorrList$Hugo_Entrez[order(MyCorrList$pvalue_spearman)[1:predictorNumber]] == colnames(working_table_TCGA_X_ExpSpTop_sub))) {stop("TCGA gene names and order are not compatible to correlation list")}
    
    ## preprocessing
    working_table_TCGA_X_ExpSpTop_sub <- predict(preProcess(working_table_TCGA_X_ExpSpTop_sub, method = c("scale", "center")), 
                                                 working_table_TCGA_X_ExpSpTop_sub)
    
    ## in case lineage is included in model 
    if ("PRIMARY_SITE" %in% colnames(model_ExpDT[[1]]$trainingData)) {
      working_table_TCGA_X_ExpSpTop_sub$PRIMARY_SITE <- lineage_list_TCGA$PRIMARY_SITE.CCLE[lineage_list_TCGA$Abbreviation == current_working_object_abb]
      working_table_TCGA_X_ExpSpTop_sub <- relocate(working_table_TCGA_X_ExpSpTop_sub, PRIMARY_SITE)
    }
    
    
    ## apply model
    #----- CHECK POINT -----#
    # if(!all(model_ExpDT[[model]]$coefnames == colnames(working_table_TCGA_X_ExpSpTop_sub))){stop("ML model gene names or order are not compatible with input data (TCGA)")}
    
    working_table_TCGA_X_ExpSpTop_predict <- data.frame(pateint = rownames(working_table_TCGA_X_ExpSpTop_sub),
                                                        predicted_DT = predict(model_ExpDT[[model]], working_table_TCGA_X_ExpSpTop_sub))
    ## annotation on patients
    # table(working_table_TCGA_patient$DSS_STATUS)  ## preview
    working_table_TCGA_X_ExpSpTop_predict$pateint <- substr(working_table_TCGA_X_ExpSpTop_predict$pateint,1,12)
    
    ## check any missing
    if(!all(working_table_TCGA_X_ExpSpTop_predict$pateint %in% working_table_TCGA_patient$PATIENT_ID)){warning("Mismatch between prediction list and meta patient list")}  
    
    working_table_TCGA_X_ExpSpTop_predict <- working_table_TCGA_X_ExpSpTop_predict %>%
      left_join(working_table_TCGA_patient, by = c("pateint" = "PATIENT_ID")) %>%
      # left_join(filter(working_table_TCGA_patient, DSS_STATUS == "1:DEAD WITH TUMOR"), by = c("pateint" = "PATIENT_ID")) %>%
      # filter(!is.na(DSS_STATUS))  %>%
      dplyr::select(c("pateint", "predicted_DT", "AGE", "SEX", 
                      "OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS",
                      "SUBTYPE","PATH_M_STAGE","PATH_N_STAGE","PATH_T_STAGE"))
    
    # table(duplicated(working_table_TCGA_X_ExpSpTop_predict$predicted_DT))
    
    # KM and Cox
    # srvl <- SurvivalData.generator(working_table_TCGA_X_ExpSpTop_predict, 2, 7, col_checkNA = 1:2, F) ## DSS as censored
    srvl <- SurvivalData.generator(working_table_TCGA_X_ExpSpTop_predict, 2, 11, col_checkNA = 3:4, F) ## PFS as censored
    
    # srvl <- SurvivalData.generator(working_table_TCGA_X_ExpSpTop_predict, 2, 7, col_checkNA = 1:2, T) ## DSS as censored, taking quantile
    # srvl <- SurvivalData.generator(working_table_TCGA_X_ExpSpTop_predict, 2, 11, col_checkNA = 3:4, T) ## PFS as censored, taking quantile
    if (nrow(srvl)>0) {
      ## KM
      # surv_object <- Surv(time = srvl$DSS_MONTHS, event = srvl$fustat) ## DSS
      surv_object <- Surv(time = srvl$PFS_MONTHS, event = srvl$fustat) ## PFS
      surv_fit_PI <- survfit(surv_object ~ pDT_group, data = srvl)
      surv_fit_sum <- summary(surv_fit_PI)
      ## p-value
      surv_diff <- survdiff(surv_object ~ pDT_group, data = srvl)
      ## plot KM curves
      surv_plot <- ggsurvplot(surv_fit_PI, data = srvl,
                              pval = TRUE, pval.method = T, pval.size = 4.5, pval.method.size = 4.5,
                              title = paste("Survival Curves (Kaplan-Meie) on", sub("_.*","",current_working_object_abb)),
                              xlab = "Time (month)",
                              ncensor.plot = TRUE,
                              ggtheme = theme_light() + theme(text = element_text(size = 13, family = "Times New Roman"),
                                                              plot.title = element_text(size = 13, hjust = 0.5)))
      
      ## cox
      surv_cox <- coxph(surv_object ~ pDT_group, data = srvl)
      surv_cox_sum <- summary(surv_cox)
      summary_surv_score_temp[[model]]$KM_pvalue[k] <- surv_diff$pvalue  ## p-value from KM curves
      summary_surv_score_temp[[model]]$KM_pvalue_neglog10[k] <- -log10(summary_surv_score_temp[[model]]$KM_pvalue[k])
      summary_surv_score_temp[[model]]$cox_wald_pvalue[k] <- surv_cox_sum$waldtest["pvalue"]  ## p-value from wald test
      summary_surv_score_temp[[model]]$cox_wald_pvalue_neglog10[k] <- -log10(summary_surv_score_temp[[model]]$cox_wald_pvalue[k])
      summary_surv_score_temp[[model]]$cox_hazard_ratio[k] <- surv_cox_sum$coefficients[,"exp(coef)"]  ## exp(coef) is the hazard ratio
      
    } else {surv_object <- NA; surv_fit_PI <- NA; surv_fit_sum <- NA; surv_diff <- NA; surv_plot <- NA;
    surv_cox <- NA; surv_cox_sum <- NA}
    
    ## data output
    surv_test_temp[[model]] <- list(srvl, surv_object, surv_fit_PI, surv_fit_sum, surv_diff, surv_plot,
                                    surv_cox, surv_cox_sum)
    names(surv_test_temp[[model]]) <- c("srvl", "surv_object", "surv_fit_PI", "surv_fit_sum", "surv_diff", "surv_plot",
                                        "surv_cox", "surv_cox_sum")
    rm(srvl, surv_object, surv_fit_PI, surv_fit_sum, surv_diff, surv_plot, surv_cox, surv_cox_sum)
    
    j <- j+1
    
    # ## data output
    # working_table_TCGA_X_ExpSpTop_temp[[model]] <- working_table_TCGA_X_ExpSpTop; rm(working_table_TCGA_X_ExpSpTop)
    # working_table_TCGA_X_ExpSpTop_sub_temp[[model]] <- working_table_TCGA_X_ExpSpTop_sub; rm(working_table_TCGA_X_ExpSpTop_sub)
    # working_table_TCGA_X_ExpSpTop_predict_temp[[model]] <- working_table_TCGA_X_ExpSpTop_predict; rm(working_table_TCGA_X_ExpSpTop_predict)
  }
  
  # working_table_TCGA_X_ExpSpTop_collection_temp[[current_working_object_abb]] <- working_table_TCGA_X_ExpSpTop_temp; rm(working_table_TCGA_X_ExpSpTop_temp)
  # working_table_TCGA_X_ExpSpTop_sub_collection_temp[[current_working_object_abb]] <- working_table_TCGA_X_ExpSpTop_sub_temp; rm(working_table_TCGA_X_ExpSpTop_sub_temp)
  # working_table_TCGA_X_ExpSpTop_predict_collection_temp[[current_working_object_abb]] <- working_table_TCGA_X_ExpSpTop_predict_temp; rm(working_table_TCGA_X_ExpSpTop_predict_temp)
  # TCGA.PCA18_predict_score_collection_temp[[current_working_object_abb]] <- TCGA.PCA18_predict_score_temp; rm(TCGA.PCA18_predict_score_temp)
  # plot_working_table_TCGA_X_ExpSpTop_predict_temp[[current_working_object_abb]] <- plot_temp; rm(plot_temp)
  surv_test_collection_temp[[current_working_object_abb]] <- surv_test_temp; rm(surv_test_temp)
}

## a data frame version of survival test score
template <- data.frame(matrix(NA,32,13)); rownames(template) <- summary_surv_score_temp[[1]]$cancer_type; colnames(template) <- c(names(summary_surv_score_temp),"n")
summary_surv_score_df_temp <- list(template,template,template,template,template); names(summary_surv_score_df_temp) <- colnames(summary_surv_score_temp[[1]])[-1]; rm(template)
for (i in 1:length(summary_surv_score_temp)) {
  summary_surv_score_df_temp$KM_pvalue[,i] <- c(summary_surv_score_temp[[i]]$KM_pvalue)
  summary_surv_score_df_temp$KM_pvalue_neglog10[,i] <- c(summary_surv_score_temp[[i]]$KM_pvalue_neglog10)
  summary_surv_score_df_temp$cox_wald_pvalue[,i] <- c(summary_surv_score_temp[[i]]$cox_wald_pvalue)
  summary_surv_score_df_temp$cox_wald_pvalue_neglog10[,i] <- c(summary_surv_score_temp[[i]]$cox_wald_pvalue_neglog10)
  summary_surv_score_df_temp$cox_hazard_ratio[,i] <- c(summary_surv_score_temp[[i]]$cox_hazard_ratio)
}
summary_surv_score_df_temp$KM_pvalue$n <- sapply(1:length(surv_test_collection_temp),function(X) length(surv_test_collection_temp[[X]]$model_ExpDT_50$surv_object))
summary_surv_score_df_temp$KM_pvalue_neglog10$n <- sapply(1:length(surv_test_collection_temp),function(X) length(surv_test_collection_temp[[X]]$model_ExpDT_50$surv_object))
summary_surv_score_df_temp$cox_wald_pvalue$n <- sapply(1:length(surv_test_collection_temp),function(X) length(surv_test_collection_temp[[X]]$model_ExpDT_50$surv_object))
summary_surv_score_df_temp$cox_wald_pvalue_neglog10$n <- sapply(1:length(surv_test_collection_temp),function(X) length(surv_test_collection_temp[[X]]$model_ExpDT_50$surv_object))
summary_surv_score_df_temp$cox_hazard_ratio$n <- sapply(1:length(surv_test_collection_temp),function(X) length(surv_test_collection_temp[[X]]$model_ExpDT_50$surv_object))

## final output

## regular method - DT gene
working_table_TCGA_X_ExpSpTop_reg_collection <- working_table_TCGA_X_ExpSpTop_collection_temp
working_table_TCGA_X_ExpSpTop_sub_reg_collection <- working_table_TCGA_X_ExpSpTop_sub_collection_temp
working_table_TCGA_X_ExpSpTop_predict_reg_collection <- working_table_TCGA_X_ExpSpTop_predict_collection_temp; saveRDS(working_table_TCGA_X_ExpSpTop_predict_reg_collection, "working_table_TCGA_X_ExpSpTop_predict_reg_collection.rds")
TCGA.PCA18_predict_score_reg_collection <- TCGA.PCA18_predict_score_collection_temp
plot_working_table_TCGA_X_ExpSpTop_predict_reg_collection <- plot_working_table_TCGA_X_ExpSpTop_predict_temp

## for DSS
surv_test_reg_DSS_collection <- surv_test_collection_temp; rm(surv_test_collection_temp)
summary_surv_score_DSS_reg <- summary_surv_score_temp; rm(summary_surv_score_temp)
summary_surv_score_DSS_df_reg <- summary_surv_score_df_temp; rm(summary_surv_score_df_temp)

## for PFS
surv_test_reg_PFS_collection <- surv_test_collection_temp; rm(surv_test_collection_temp)
summary_surv_score_PFS_reg <- summary_surv_score_temp; rm(summary_surv_score_temp)
summary_surv_score_PFS_df_reg <- summary_surv_score_df_temp; rm(summary_surv_score_df_temp)

## summary with 50-model
temp_KMlog10 <- cbind(summary_surv_score_DSS_df_reg$KM_pvalue_neglog10$model_ExpDT_50,
                      summary_surv_score_DSS_Q_df_reg$KM_pvalue_neglog10$model_ExpDT_50,
                      summary_surv_score_PFS_df_reg$KM_pvalue_neglog10$model_ExpDT_50,
                      summary_surv_score_PFS_Q_df_reg$KM_pvalue_neglog10$model_ExpDT_50) %>% as.data.frame()
rownames(temp_KMlog10) <- rownames(summary_surv_score_DSS_df_reg$KM_pvalue_neglog10)
colnames(temp_KMlog10) <- c("DSS","DSS_Q","PFS","PFS_Q")

temp_coxlog10 <- cbind(summary_surv_score_DSS_df_reg$cox_wald_pvalue_neglog10$model_ExpDT_50,
                       summary_surv_score_DSS_Q_df_reg$cox_wald_pvalue_neglog10$model_ExpDT_50,
                       summary_surv_score_PFS_df_reg$cox_wald_pvalue_neglog10$model_ExpDT_50,
                       summary_surv_score_PFS_Q_df_reg$cox_wald_pvalue_neglog10$model_ExpDT_50) %>% as.data.frame()
rownames(temp_coxlog10) <- rownames(summary_surv_score_DSS_df_reg$cox_wald_pvalue_neglog10)
colnames(temp_coxlog10) <- c("DSS","DSS_Q","PFS","PFS_Q")


### MKI67 as biomarker (for GENT2) - survival 
summary_surv_score_temp <- data.frame(matrix(NA,length(survival_type_set),6)); colnames(summary_surv_score_temp) <- c("KM_pvalue", "KM_pvalue_neglog10", "cox_wald_pvalue", "cox_wald_pvalue_neglog10", "cox_hazard_ratio","n")
surv_test_temp <- list()
for(j in 1:length(survival_type_set)){
  survival_type <- survival_type_set[j]
  
  srvl <- GeneSymbolP2_df_MKI67_df
  colnames(srvl)[colnames(srvl) == "Primary_site"] <- "Tissue"
  if(survival_type == "OS") {
    srvl <- SurvivalData.generator(srvl, 2, 8, col_checkNA = 5:6, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
    surv_object <- Surv(time = srvl$OS_MONTHS, event = srvl$fustat)
  } else if(survival_type == "DFS") {
    srvl <- SurvivalData.generator(srvl, 2, 10, col_checkNA = 7:8, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
    surv_object <- Surv(time = srvl$DFS_MONTHS, event = srvl$fustat)
  } else if(survival_type == "DSS") {
    srvl <- SurvivalData.generator(srvl, 2, 12, col_checkNA = 1:2, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
    surv_object <- Surv(time = srvl$DSS_MONTHS, event = srvl$fustat)
  } else if(survival_type == "PFS") {
    srvl <- SurvivalData.generator(srvl, 2, 14, col_checkNA = 3:4, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
    surv_object <- Surv(time = srvl$PFS_MONTHS, event = srvl$fustat)
  }
  
  if (nrow(srvl)>0) {
    ### MKI67
    # surv_fit_PI <- survfit(surv_object ~ MKI67_Group, data = srvl)
    surv_fit_PI <- survfit(surv_object ~ MKI67_Group + Tissue, data = srvl)
    surv_fit_sum <- summary(surv_fit_PI)
    # p-value
    # surv_diff <- survdiff(surv_object ~ MKI67_Group, data = srvl)
    surv_diff <- survdiff(surv_object ~ MKI67_Group + Tissue, data = srvl)
    # cox
    # surv_cox <- coxph(surv_object ~ MKI67_Group, data = srvl)  ## MKI67
    surv_cox <- coxph(surv_object ~ MKI67_Group + Tissue, data = srvl)  ## MKI67
    surv_cox_sum <- summary(surv_cox)
    
    ## plot KM curves
    names(surv_fit_PI$strata) <- gsub("_"," ",names(surv_fit_PI$strata)) ## changing "MKI67_Group=Low KI67" into "MKI67 group=Low KI67". can also apply for other cases
    names(surv_fit_PI$strata) <- gsub("="," = ",names(surv_fit_PI$strata)) ## changing "MKI67 group=Low KI67" into "MKI67 group = Low KI67". can also apply for other cases
    surv_plot <- ggsurvplot(surv_fit_PI, data = srvl,
                            # pval = TRUE, pval.method = T, pval.size = 4.5, pval.method.size = 4.5,
                            # title = paste("Survival Curves (Kaplan-Meie) (",survival_type,") on Cancers from GENT2", sep = ""),  ## whole data
                            title = paste("Survival Curves (Kaplan-Meie) (",survival_type,") Based on MKI67 Level", sep = ""),  ## MKI67
                            subtitle = paste("n =", length(surv_object)),
                            xlab = "Time (Month)", ylab = "Survival Probability",
                            ncensor.plot = TRUE,
                            ggtheme = theme_light() + theme(text = element_text(size = 18, family = "Times New Roman"),
                                                            plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 17, hjust = 0.5),
                                                            axis.title = element_text(size = 18), axis.text = element_text(size = 17), legend.text = element_text(size = 17),
                                                            axis.text.y = element_text(angle = 45)))
    surv_plot$ncensor.plot$labels$title <- "Number of Censoring"  ## normalize the censor plot title
    surv_plot$plot <- surv_plot$plot + annotate("text",label= paste("Log-rank\n",signif(surv_diff$pvalue,digits=4),sep=""),x=max(surv_fit_PI$time)*0.75,y=0.9, size = 6, family = "Times New Roman")
    ## for multi-variable
    surv_plot$plot <- surv_plot$plot + guides(col = guide_legend(ncol = 2))
    
    summary_surv_score_temp$KM_pvalue[j] <- surv_diff$pvalue  ## p-value from KM curves
    summary_surv_score_temp$KM_pvalue_neglog10[j] <- -log10(summary_surv_score_temp$KM_pvalue[j])
    summary_surv_score_temp$cox_wald_pvalue[j] <- surv_cox_sum$waldtest["pvalue"]  ## p-value from wald test
    summary_surv_score_temp$cox_wald_pvalue_neglog10[j] <- -log10(summary_surv_score_temp$cox_wald_pvalue[j])
    # summary_surv_score_temp$cox_hazard_ratio[j] <- surv_cox_sum$coefficients[,"exp(coef)"]  ## exp(coef) is the hazard ratio
    
  }else {cat(survival_type, "doesn't have infor.\n");
    surv_object <- NA; surv_fit_PI <- NA; surv_fit_sum <- NA; surv_diff <- NA; surv_plot <- NA;
    surv_cox <- NA; surv_cox_sum <- NA}
  
  summary_surv_score_temp$n[j] <- nrow(srvl)
  rownames(summary_surv_score_temp)[j] <- survival_type
  
  ## data output
  surv_test_temp[[survival_type]] <- list(srvl, surv_object, surv_fit_PI, surv_fit_sum, surv_diff, surv_plot,
                                          surv_cox, surv_cox_sum)
  names(surv_test_temp[[survival_type]]) <- c("srvl", "surv_object", "surv_fit_PI", "surv_fit_sum", "surv_diff", "surv_plot",
                                              "surv_cox", "surv_cox_sum")
  rm(srvl, surv_object, surv_fit_PI, surv_fit_sum, surv_diff, surv_plot, surv_cox, surv_cox_sum)
}
## output
## one variable
## MKI67
summary_surv_score_GENT2_MKI67 <- summary_surv_score_temp; surv_test_GENT2_MKI67 <- surv_test_temp; rm(summary_surv_score_temp, surv_test_temp)
surv_test_GENT2_MKI67$OS$surv_plot$ncensor.plot <- surv_test_GENT2_MKI67$OS$surv_plot$ncensor.plot + scale_y_continuous(breaks = c(seq(0,11,by = 4),11)) ## changing the y axis manually to avoid over-crowded
saveRDS(summary_surv_score_GENT2_MKI67, surv_test_GENT2_MKI67, "surv_GENT2_MKI67.rds")
## multi-variable
## MKI67
summary_surv_score_GENT2_MKI67_2var <- summary_surv_score_temp; surv_test_GENT2_MKI67_2var <- surv_test_temp; rm(summary_surv_score_temp, surv_test_temp)
surv_test_GENT2_MKI67_2var$OS$surv_plot$ncensor.plot <- surv_test_GENT2_MKI67_2var$OS$surv_plot$ncensor.plot + scale_y_continuous(breaks = c(seq(0,11,by = 4),11)) ## changing the y axis manually to avoid over-crowded



### KM - regular model, PFS, PAAD
plot_data <- surv_test_reg_PFS_collection$PAAD_TCGA$model_ExpDT_50; survival_type<- "PFS"; Tissue <- "PAAD"
surv_object<-plot_data$surv_object
names(plot_data$surv_fit_PI$strata) <- gsub("_"," ",names(plot_data$surv_fit_PI$strata)) ## changing "MKI67_Group=Low KI67" into "MKI67 group=Low KI67". can also apply for other cases
names(plot_data$surv_fit_PI$strata) <- gsub("="," = ",names(plot_data$surv_fit_PI$strata)) ## changing "MKI67 group=Low KI67" into "MKI67 group = Low KI67". can also apply for other cases
if(plot_data$surv_diff$pvalue<0.01){
  input_bquote <- bquote(paste(.(sub("e.*","",signif(plot_data$surv_diff$pvalue,digits=4)%>%scientific())%>%as.numeric())
                               %*% 10^.(sub(".*e","",signif(plot_data$surv_diff$pvalue,digits=4)%>%scientific())%>%as.numeric())))
} else {input_bquote <- signif(plot_data$surv_diff$pvalue,digits=4)}


surv_plot <- ggsurvplot(plot_data$surv_fit_PI, data = plot_data$srvl,
                        title = paste("Kaplan-Meier Curves (",survival_type,") on ",Tissue," Cancer", sep = ""),  ## Model
                        subtitle = paste("n =", length(plot_data$surv_object)),
                        xlab = "Time (Month)", ylab = "Survival Probability",
                        ncensor.plot = TRUE,
                        ggtheme = theme_classic() + theme(text = element_text(size = 19, family = "Times New Roman"),
                                                          plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 18, hjust = 0.5),
                                                          axis.title = element_text(size = 19), axis.text = element_text(size = 18), legend.text = element_text(size = 18),
                                                          axis.text.y = element_text(angle = 45))+
                          theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
                                axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm")) ## axis title 10 units from two sides; right plot margin=0.7cm
)
surv_plot$ncensor.plot$labels$title <- "Number of Censoring"  ## normalize the censor plot title
surv_plot$ncensor.plot$labels$y <- "n"  ## normalize the censor plot title
surv_plot$plot <- surv_plot$plot + 
  annotate("text",label= paste("Log-rank",sep=""),x=max(plot_data$surv_fit_PI$time)*0.75,y=0.9, size = 6, family = "Times New Roman")+
  annotate("text",label= input_bquote,x=max(plot_data$surv_fit_PI$time)*0.75,y=0.75, size = 6, family = "Times New Roman")

hold_plot <- surv_plot[1:2]
for(i in 1:length(hold_plot)){ ## remove legend from censor plot
  if(i%%2==0){hold_plot[[i]]<-hold_plot[[i]]+theme(legend.position="none")} }

## This integrated KM plot is Figure 6A ##
ggexport(
  marrangeGrob(hold_plot, 
               layout_matrix = layout_srvl_withcensor_2in1[,1:4]),
  filename = "TCGA_PAAD_50_KM.pdf",
  width = 8, height = 8.5)#; rm(surv_object)


## survival analysis barplot
### focus on 50 only
## cancer type vs vs -log10Cox-pvalue
## regular
plot_table <- summary_surv_score_PFS_reg$model_ExpDT_50 ## PFS
plot_table$cox_wald_qvalue <- p.adjust(plot_table$cox_wald_pvalue,method = "BH")
plot_table$cox_wald_qvalue_neglog10 <- -log10(plot_table$cox_wald_qvalue)
plot_table$cancer_type <- sub("_.*","",plot_table$cancer_type)
## construct the data for plotting: med pDT vs -log10Cox-pvalue
median_temp <- lapply(working_table_TCGA_X_ExpSpTop_predict_reg_collection,"[[",8); median_temp <- median_temp[!is.na(median_temp)]

median_temp <- lapply(median_temp, "[[", 2)
plot_table2 <- data.frame(pDT_median = sapply(median_temp, median)) %>% rownames_to_column("cancer_type")
plot_table2$cancer_type <- sub("_.*","",plot_table2$cancer_type)
plot_table2 <- plot_table2 %>% left_join(plot_table[,c("cancer_type","cox_wald_qvalue_neglog10")]) %>% drop_na(cox_wald_qvalue_neglog10)

plot_table <- plot_table %>% filter(!(cancer_type %in% c("LAML"))) ## for reg 
plot_table$cox_hazard_ratio[plot_table$cancer_type%in%c("DLBCell")] <- NA  ## for reg ## for PFS
plot_table$cancer_type[plot_table$cancer_type == "DLBCell"] <- "DLBC"

## This is Figure 6B ##
plot_hold <-
  ggplot(plot_table, aes(cancer_type,cox_wald_qvalue_neglog10)) + geom_bar(stat="identity", aes(fill=cox_hazard_ratio))+coord_flip()+
  scale_fill_gradient(low="purple", high="red") +
  scale_x_discrete(limits = plot_table$cancer_type[order(plot_table$cox_wald_qvalue_neglog10, decreasing = F)[(nrow(plot_table)-16):nrow(plot_table)]]) +
  geom_hline(yintercept = -log10(0.1), linetype="dashed", color = "red") +
  theme_classic() + theme(text = element_text(size = 18, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 19, hjust = 0.5),
                          axis.text = element_text(size = 18),legend.position = c(0.99,0.01), legend.justification = c("right", "bottom"), legend.text = element_text(size = 19)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x="Cancer Type", y=bquote('-log' ['10'] ('Cox q-value')), title = "Cox Regression on PFS among Diff. Cancers", fill = "Hazard Ratio")+
  scale_y_continuous(limits = c(0,2.5), breaks = seq(0,3,0.5))+  ## for reg. PFS
  labs(subtitle = "50's Classical Model Prediction")

### DTLearner for SNV #####
##------- SNV for all TCGA
TCGA_SNV_count_pDT_reg_collection <- list()
FIG_SMutate <- list()
summary_TCGA_SNV_count_pDT <- data.frame(matrix(NA,length(TCGA.PCA18.cbp_mutation_RAW),12))
rownames(summary_TCGA_SNV_count_pDT)<- names(TCGA.PCA18.cbp_mutation_RAW)
colnames(summary_TCGA_SNV_count_pDT) <- c("spearman_p_reg","spearman_est_reg","n_reg","spearman_p_fcl","spearman_est_fcl","n_fcl","spearman_p_sfl","spearman_est_sfl","n_sfl","spearman_p_blnc","spearman_est_blnc","n_blnc")
for (i in 1:length(TCGA.PCA18.cbp_mutation_RAW)) {
  current_working_object_abb <- names(TCGA.PCA18.cbp_mutation_RAW)[i]
  working_table_TCGA <- TCGA.PCA18.cbp_exp.rsem[[current_working_object_abb]]
  working_table_TCGA_patient <- TCGA.PCA18.cbp_patient_RAW[[current_working_object_abb]]
  
  ## SNV on TCGA
  TCGA_SNV <- TCGA.PCA18.cbp_mutation_RAW[[current_working_object_abb]] %>% 
    filter(Variant_Type == "SNP") %>% 
    dplyr::select(c("Hugo_Symbol","Entrez_Gene_Id","Variant_Classification","Variant_Type","Tumor_Sample_Barcode","Transcript_ID"))
  
  TCGA_SNV$Patient_ID <- substr(TCGA_SNV$Tumor_Sample_Barcode,1,12)
  
  ## counting the number of SNV on each genes for each patients
  list_gene <- unique(TCGA_SNV$Hugo_Symbol)
  list_patient <- unique(TCGA_SNV$Patient_ID)
  ## everytime the SNV is counted via table(), the data from each patient will be added into df via left_join()
  TCGA_SNV_count <- data.frame(matrix(NA, length(list_gene), 1))
  colnames(TCGA_SNV_count) <- "Hugo_Symbol"
  TCGA_SNV_count$Hugo_Symbol <- list_gene
  for (patient in list_patient) {  ## treat on patient level
    df_patient <- TCGA_SNV %>% filter(Patient_ID == patient)  ## extract the data from current patient
    TCGA_SNV_count <- TCGA_SNV_count %>% left_join(data.frame(table((df_patient$Hugo_Symbol))), by = c("Hugo_Symbol" = "Var1"))
    colnames(TCGA_SNV_count)[ncol(TCGA_SNV_count)] <- patient  ## named with patient ID
  }; rm(patient,df_patient)
  
  ## table of sum of SNV and pDT
  #### regular model
  TCGA_SNV_count_pDT_reg <- data.frame(sum_SNV = apply(TCGA_SNV_count[,-1],2,function(X) sum(X, na.rm = T))) %>%
    rownames_to_column(var = "Patient_ID") %>%
    left_join(predicitons_ExpDT_reg_collection_50[[paste(current_working_object_abb,"_TCGA",sep="")]][,c("pateint", "predicted_DT", "SUBTYPE")], by = c("Patient_ID" = "pateint")) %>%
    drop_na(predicted_DT)
  TCGA_SNV_count_pDT_reg$sum_SNV_log10 <- log10(TCGA_SNV_count_pDT_reg$sum_SNV)  ## log10 on SNV
  TCGA_SNV_count_pDT_reg <- TCGA_SNV_count_pDT_reg[,c(1,2,5,3,4)]  ## reorder columns
  TCGA_SNV_count_pDT_reg$SUBTYPE[TCGA_SNV_count_pDT_reg$SUBTYPE == ""] <- "N.A."  ## for those subtype marked with nothing - change into N.A.
  
  label_location <- c(x = max(TCGA_SNV_count_pDT_reg$sum_SNV_log10)*0.8, y = max(TCGA_SNV_count_pDT_reg$predicted_DT)*0.98)
  if ((filter(TCGA_SNV_count_pDT_reg,SUBTYPE!="N.A.") %>% count(SUBTYPE) %>% nrow()) > 1) { ## in case there are subtypes ( =1 means subtype is the cancer type itself)
    ind <- which(TCGA_SNV_count_pDT_reg$SUBTYPE != "N.A.")
    # ggplot_recipe <- ggplot(dplyr::filter(TCGA_SNV_count_pDT_reg, SUBTYPE != "N.A."), aes(sum_SNV_log10,predicted_DT)) +  ## if N.A. is not included
    ggplot_recipe <- ggplot(TCGA_SNV_count_pDT_reg, aes(sum_SNV_log10,predicted_DT)) +  ## include everything
      geom_point(aes(color=SUBTYPE), size = 1.5, alpha = 0.7) + scale_color_brewer(palette = 'Set1')
  } else { ## in case there is no subtypes
    ggplot_recipe <- ggplot(TCGA_SNV_count_pDT_reg, aes(sum_SNV_log10,predicted_DT)) +
      geom_point(size = 1.5, alpha = 0.7)
  }
  
  ## corr. test
  summary_TCGA_SNV_count_pDT$spearman_p_reg[i] <- cor.test(TCGA_SNV_count_pDT_reg$sum_SNV_log10, TCGA_SNV_count_pDT_reg$predicted_DT, method = "spearman")$p.value
  summary_TCGA_SNV_count_pDT$spearman_est_reg[i] <- cor.test(TCGA_SNV_count_pDT_reg$sum_SNV_log10, TCGA_SNV_count_pDT_reg$predicted_DT, method = "spearman")$est
  summary_TCGA_SNV_count_pDT$n_reg[i] <- nrow(TCGA_SNV_count_pDT_reg)  ## total cases
  
  FIG_SMutate[[(i*4-3)]] <- ggplot_recipe +
    theme_light() + theme(text = element_text(size = 18, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5)) +
    labs(x="Total SNV # in each patient in log10",y="pDT", title = paste("Total SNV Number vs Predicted DT on", current_working_object_abb, "with Classical Model")) +
    # geom_point(aes(shape=SUBTYPE), alpha = 0.875) +
    geom_smooth(method=lm, se=FALSE, col='grey', linewidth = 0.5) +
    annotate("text",family = "Times New Roman", x=label_location[1], y=label_location[2], 
             label = paste("Spr. p-value =",
                           signif(summary_TCGA_SNV_count_pDT$spearman_p_reg[i], 4),
                           ",\nRho = ", signif(summary_TCGA_SNV_count_pDT$spearman_est_reg[i], 4))) +
    stat_regline_equation(aes(label = paste(after_stat(eq.label), after_stat(rr.label),sep = "~~~~")),
                          label.x = label_location[1]*0.9, label.y = (label_location[2]*0.97),family = "Times New Roman")
  rm(label_location)
  
  names(FIG_SMutate)[[(i*4-3)]] <- paste(current_working_object_abb,"reg",sep = "_")
  
  ## data output
  TCGA_SNV_count_pDT_reg_collection[[i]] <- TCGA_SNV_count_pDT_reg
  names(TCGA_SNV_count_pDT_reg_collection)[i] <- current_working_object_abb
}

summary_TCGA_SNV_count_pDT$spearman_q_reg <- p.adjust(summary_TCGA_SNV_count_pDT$spearman_p_reg, method = "BH")

# plot_data <- FIG_SMutate$STAD_reg$data
plot_data <- TCGA_SNV_count_pDT_reg_collection$STAD; current_working_object_abb="STAD"
ptext_loca <- axis_location(plot_data,3,4,"rightup"); ptext_loca[2] <- log10(ptext_loca[2])
input_bquote <- bquote(paste('Spr. q-value ='~ .(sub("e.*","",signif(summary_TCGA_SNV_count_pDT[current_working_object_abb,"spearman_q_reg"],digits=4))%>%as.numeric())
                             %*% 10^.(sub(".*e","",signif(summary_TCGA_SNV_count_pDT[current_working_object_abb,"spearman_q_reg"],digits=4))%>%as.numeric())))

hold_plot <-
  ggplot(plot_data, aes(sum_SNV_log10,log10(predicted_DT))) + geom_point(color = "dodgerblue4",size=2,alpha=0.7) + scale_fill_distiller(palette = "Blues")+
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), legend.position = "none", plot.subtitle = element_text(size = 20, hjust = 0.5)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm 
  labs(x=bquote('log' ['10'] ('Total SNV Number')),y=bquote('Predicted Doubling Time log' ['10'] ('hr')), title = paste("SNV - Predicted DT on", current_working_object_abb)) +
  labs(subtitle = "Classical Model") +
  geom_smooth(method=lm, se=FALSE, col='grey', linewidth = 1, linetype = 2) +
  scale_x_continuous(limits = c(min(plot_data$sum_SNV_log10),max(plot_data$sum_SNV_log10)*1.06)) +
  
  stat_regline_equation(label.x = ptext_loca[1]*0.7, label.y = ptext_loca[2]*1.07, size = 6,family = "Times New Roman")+ ## trendline equation
  stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
           label.x = ptext_loca[1]*0.83, label.y = ptext_loca[2]*1.055, size = 6,family = "Times New Roman") + ## R^2
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*1.035, size = 6,
           label = input_bquote)+  ## qvalue
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*1.015, size = 6,
           label = paste("Rho =",signif(summary_TCGA_SNV_count_pDT[current_working_object_abb,"spearman_est_reg"],digits=4)))  ## est.

## This is Figure 4A ##
ggsave(filename = "SNV_50_model_STAD.pdf",
       hold_plot,
       width = 6.5, height = 6
)

## histogram of SNV
## classical
plot_table <- summary_TCGA_SNV_count_pDT %>% rownames_to_column("cancer_type") %>% select(cancer_type,spearman_q_reg, n_reg)
plot_table$spearman_q_reg <- -log10(plot_table$spearman_q_reg)
plot_table$n_reg <- log10(plot_table$n_reg)

hold_plot <-
  ggplot(plot_table,aes(spearman_q_reg,cancer_type)) + geom_bar(stat="identity", aes(fill=n_reg)) + scale_fill_gradient() +
  scale_y_discrete(limits = plot_table$cancer_type[order(plot_table$spearman_q_reg, decreasing = F)]) +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color = "red", alpha = 0.8) +
  theme_classic() + theme(text = element_text(size = 18, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), legend.position = c(0.99,0.01), legend.justification = c("right", "bottom")) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(y ="Cancer Type", x =bquote('-log' ['10'] ('Spearman q-value')), title = "Spearman Corr. between pDT - SNV\nClassical Model", fill = bquote('log' ['10'] ('Sample Size')))

## This is Figure 4B ##
ggsave(filename = "plot/SNV_50_model_spearman_barplot.pdf",
       hold_plot,
       width = 6.5, height = 7)


### DTLearner for cancer stages #####
## df for result storage
predicitons_ExpDT_reg_collection_50
predicitons_ExpDT_reg_collection_50_df
pDT_TCGA_exp_stage_regular <- data.frame(matrix(NA, length(predicitons_ExpDT_reg_collection_50), 7))
colnames(pDT_TCGA_exp_stage_regular) <- c("cancer_type","T4_T1_t.test","N1_N0_t.test","M1_M0_t.test","T4_T1_wilcox","N1_N0_wilcox","M1_M0_wilcox")
pDT_TCGA_exp_stage_regular[,1] <- names(predicitons_ExpDT_reg_collection_50)

for(i in 1:nrow(pDT_TCGA_exp_stage_regular)) {
  group_low <- predicitons_ExpDT_reg_collection_50[[i]][grep("T1",predicitons_ExpDT_reg_collection_50[[i]]$PATH_T_STAGE),1:2]
  group_high <- predicitons_ExpDT_reg_collection_50[[i]][grep("T4",predicitons_ExpDT_reg_collection_50[[i]]$PATH_T_STAGE),1:2]
  if(nrow(group_low) >= 2 & nrow(group_high)>=2) {
    pDT_TCGA_exp_stage_regular[i,2] <- t.test(group_low[,2], group_high[,2], alternative = "two.sided")$p.value
    pDT_TCGA_exp_stage_regular[i,5] <- wilcox.test(group_low[,2], group_high[,2], alternative = "two.sided")$p.value
  } else {pDT_TCGA_exp_stage_regular[i,c(2,5)] <- NA}
  
  group_low <- predicitons_ExpDT_reg_collection_50[[i]][grep("N0",predicitons_ExpDT_reg_collection_50[[i]]$PATH_N_STAGE),1:2]
  group_high <- predicitons_ExpDT_reg_collection_50[[i]][grep("N1",predicitons_ExpDT_reg_collection_50[[i]]$PATH_N_STAGE),1:2]
  if(nrow(group_low) >= 2 & nrow(group_high)>=2) {
    pDT_TCGA_exp_stage_regular[i,3] <- t.test(group_low[,2], group_high[,2], alternative = "two.sided")$p.value
    pDT_TCGA_exp_stage_regular[i,6] <- wilcox.test(group_low[,2], group_high[,2], alternative = "two.sided")$p.value
  } else {pDT_TCGA_exp_stage_regular[i,c(3,6)] <- NA}
  
  group_low <- predicitons_ExpDT_reg_collection_50[[i]][grep("M0",predicitons_ExpDT_reg_collection_50[[i]]$PATH_M_STAGE),1:2]
  group_high <- predicitons_ExpDT_reg_collection_50[[i]][grep("M1",predicitons_ExpDT_reg_collection_50[[i]]$PATH_M_STAGE),1:2]
  if(nrow(group_low) >= 2 & nrow(group_high)>=2) {
    pDT_TCGA_exp_stage_regular[i,4] <- t.test(group_low[,2], group_high[,2], alternative = "two.sided")$p.value
    pDT_TCGA_exp_stage_regular[i,7] <- wilcox.test(group_low[,2], group_high[,2], alternative = "two.sided")$p.value
  } else {pDT_TCGA_exp_stage_regular[i,c(4,7)] <- NA}
  
}; rm(group_low, group_high, i)

## pick one of the following:
## for Figure 5B:
plot_data <- pDT_TCGA_exp_stage_regular[,c("cancer_type","T4_T1_t.test")]; colnames(plot_data)[2] ="stages"; type_stage="T4/T1 Stages"; type_test <- "t-test"
## for Figure 5D:
plot_data <- pDT_TCGA_exp_stage_regular[,c("cancer_type","N1_N0_t.test")]; colnames(plot_data)[2] ="stages"; type_stage="N1/N0 Stages"; type_test <- "t-test"
## for Figure 5F:
plot_data <- pDT_TCGA_exp_stage_regular[,c("cancer_type","M1_M0_t.test")]; colnames(plot_data)[2] ="stages"; type_stage="M1/M0 Stages"; type_test <- "t-test"

## Then run the following:
plot_data$cancer_type <- gsub("_TCGA","",plot_data$cancer_type)
plot_data <- plot_data %>% filter(!is.na(stages))
plot_data$stages <- p.adjust(plot_data$stages, "BH")
plot_data$stages <- -log10(plot_data$stages)

## This is Figure 5B/D/F ##
ggplot(plot_data, aes(stages,cancer_type)) + geom_bar(stat="identity", fill="grey50") +
  scale_y_discrete(limits = plot_data$cancer_type[order(plot_data$stages, decreasing = F)]) +
  geom_vline(xintercept = -log10(0.1), linetype="dashed", color = "red") +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('-log' ['10'] ('q-value')), y = "Cancer Type",
       title = paste(type_test, "on pDT btween", type_stage, sep=" "))



## pick one of the following:
## for Figure 5A:
plot_data <- predicitons_ExpDT_reg_collection_50$KIRP_TCGA; type_cancer = "KIRP"
plot_data$PATH_T_STAGE_simple <- substr(plot_data$PATH_T_STAGE, 1,2)
plot_data_sub <- plot_data[,c("PATH_T_STAGE_simple","predicted_DT")]; colnames(plot_data_sub)[1] ="stages"; type_stage="T Stage"
## for Figure 5C:
plot_data <- predicitons_ExpDT_reg_collection_50$LUAD_TCGA; type_cancer = "LUAD"
plot_data$PATH_N_STAGE_simple <- substr(plot_data$PATH_N_STAGE, 1,2)
plot_data_sub <- plot_data[,c("PATH_N_STAGE_simple","predicted_DT")]; colnames(plot_data_sub)[1] ="stages"; type_stage="N Stage"
## for Figure 5E:
plot_data <- predicitons_ExpDT_reg_collection_50$KIRP_TCGA; type_cancer = "KIRP"
plot_data$PATH_M_STAGE_simple <- substr(plot_data$PATH_M_STAGE, 1,2)
plot_data_sub <- plot_data[,c("PATH_M_STAGE_simple","predicted_DT")]; colnames(plot_data_sub)[1] ="stages"; type_stage="M Stage"

## Then run the following: (silence and un-silence the code)
plot_data_sub <- plot_data_sub %>% filter(!is.na(stages),nchar(stages)!=0)
plot_data_sub <- plot_data_sub[-grep("X",plot_data_sub$stages),]

## This is Figure 5A/C/E ##
ggplot(plot_data_sub, aes(stages,predicted_DT)) +
  ## for T stage:
  # geom_violin(data=plot_data_sub%>%filter(stages!="T4"),aes(fill = stages), width=0.8) + geom_boxplot(aes(fill = stages),width=0.1) + 
  ## for N and M stage:
  geom_violin(data=plot_data_sub,aes(fill = stages), width=1) + geom_boxplot(aes(fill = stages),width=0.1) + 
  geom_jitter(shape=1, position=position_jitter(0.15), alpha=0.7)+
  scale_y_continuous(limits = c(min(plot_data_sub$predicted_DT),max(plot_data_sub$predicted_DT)))+
  theme_classic()+
  theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5),
        axis.text =  element_text(size = 20), axis.text.x = element_text(size =20), #axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),
        legend.position = "none")  +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,1,0.5,0.5), "cm"))+ ## axis title 10 units from two sides; right plot margin=1cm
  labs(y= "pDT hr", x = type_stage, 
       title = paste("pDT on ", type_stage,"s", " from ",type_cancer, sep=""))

##------- GENT2 data treatment and prediciton -----------
## the object GeneSymbolP2_df_sub_df has included all the 50 genes for the training from all samples
## apply model
model <- model_ExpDT_gent2_reg$model_ExpDT_50
## rescale
processCenter <- preProcess(GeneSymbolP2_df_sub_df, method = c("scale", "center"))
GeneSymbolP2_df_sub_df <- predict(processCenter, GeneSymbolP2_df_sub_df)

## predict
colnames(GeneSymbolP2_df_sub_df)[colnames(GeneSymbolP2_df_sub_df)=="`KRTAP5-2_440021`"] <- "KRTAP5-2_440021" ## change the format to avoid bug
GeneSymbolP2_df_predict <- data.frame(GSM_samples = unlist(list_sample_P2), 
                                      pDT = predict(model, GeneSymbolP2_df_sub_df))
## annotate with patient data
GeneSymbolP2_df_predict <- GeneSymbolP2_df_predict %>% left_join(sql_survival[,c("Sample_id","Primary_site","Subtype","Stage","TNM_Stage","OS_Month","OS_Status",
                                                                                 "DFS_Month","DFS_Status",
                                                                                 "DSS_Month","DSS_Status",
                                                                                 "PFS_Month","PFS_Status")], by = c(GSM_samples = "Sample_id"))
ind <- apply(GeneSymbolP2_df_predict[,-1:-2], 1, function(X) !all(is.na(X))) %>% which()  ## remove samples without available clinical information
GeneSymbolP2_df_predict <- GeneSymbolP2_df_predict[ind,]

## change "None" into NA, and make the number into numeric
for(i in 3:ncol(GeneSymbolP2_df_predict)){
  ind <- which(GeneSymbolP2_df_predict[,i] == "None")
  GeneSymbolP2_df_predict[ind,i] <- NA
}
GeneSymbolP2_df_predict[,c("OS_Month")] <- as.numeric(GeneSymbolP2_df_predict[,c("OS_Month")])
GeneSymbolP2_df_predict[,c("DFS_Month")] <- as.numeric(GeneSymbolP2_df_predict[,c("DFS_Month")])
GeneSymbolP2_df_predict[,c("DSS_Month")] <- as.numeric(GeneSymbolP2_df_predict[,c("DSS_Month")])
GeneSymbolP2_df_predict[,c("PFS_Month")] <- as.numeric(GeneSymbolP2_df_predict[,c("PFS_Month")])

## make the column name consistent
colnames(GeneSymbolP2_df_predict)[colnames(GeneSymbolP2_df_predict) %in% c("OS_Month","OS_Status","DFS_Month","DFS_Status","DSS_Month","DSS_Status","PFS_Month","PFS_Status")] <- 
  c("OS_MONTHS","OS_STATUS","DFS_MONTHS","DFS_STATUS","DSS_MONTHS","DSS_STATUS","PFS_MONTHS","PFS_STATUS")

## This is the prediction (pDT) for GENT2 samples
GeneSymbolP2_df_predict


## similarly, for MKI67 instead of pDT index
## MKI67 expression level is extracted already, and it is stored it in three data frame
GeneSymbolP2_df_MKI67_df <- rbind(GeneSymbolP2_df_MKI67[[1]],GeneSymbolP2_df_MKI67[[2]],GeneSymbolP2_df_MKI67[[3]])
GeneSymbolP2_df_MKI67_df <- GeneSymbolP2_df_MKI67_df %>% rownames_to_column("GSM_samples")

## annotate with patient data
GeneSymbolP2_df_MKI67_df <- GeneSymbolP2_df_MKI67_df %>% left_join(sql_survival[,c("Sample_id","Primary_site","Subtype","Stage","TNM_Stage","OS_Month","OS_Status",
                                                                                   "DFS_Month","DFS_Status",
                                                                                   "DSS_Month","DSS_Status",
                                                                                   "PFS_Month","PFS_Status")], by = c(GSM_samples = "Sample_id"))

ind <- apply(GeneSymbolP2_df_MKI67_df[,-1:-2], 1, function(X) !all(is.na(X))) %>% which()  ## remove samples without available clinical information
GeneSymbolP2_df_MKI67_df <- GeneSymbolP2_df_MKI67_df[ind,]; rm(ind)

## change "None" into NA, and make the number into numeric
for(i in 3:ncol(GeneSymbolP2_df_MKI67_df)){
  ind <- which(GeneSymbolP2_df_MKI67_df[,i] == "None")
  GeneSymbolP2_df_MKI67_df[ind,i] <- NA
}
GeneSymbolP2_df_MKI67_df[,c("OS_Month")] <- as.numeric(GeneSymbolP2_df_MKI67_df[,c("OS_Month")])
GeneSymbolP2_df_MKI67_df[,c("DFS_Month")] <- as.numeric(GeneSymbolP2_df_MKI67_df[,c("DFS_Month")])
GeneSymbolP2_df_MKI67_df[,c("DSS_Month")] <- as.numeric(GeneSymbolP2_df_MKI67_df[,c("DSS_Month")])
GeneSymbolP2_df_MKI67_df[,c("PFS_Month")] <- as.numeric(GeneSymbolP2_df_MKI67_df[,c("PFS_Month")])

## make the column name consistent
colnames(GeneSymbolP2_df_MKI67_df)[colnames(GeneSymbolP2_df_MKI67_df) %in% c("OS_Month","OS_Status","DFS_Month","DFS_Status","DSS_Month","DSS_Status","PFS_Month","PFS_Status")] <- 
  c("OS_MONTHS","OS_STATUS","DFS_MONTHS","DFS_STATUS","DSS_MONTHS","DSS_STATUS","PFS_MONTHS","PFS_STATUS")

## This is the KI-67 level for GENT2 samples
GeneSymbolP2_df_MKI67_df


## survival test
survival_type_set <- c("OS","DFS","DSS","PFS")
summary_surv_score_temp <- data.frame(matrix(NA,length(survival_type_set),6)); colnames(summary_surv_score_temp) <- c("KM_pvalue", "KM_pvalue_neglog10", "cox_wald_pvalue", "cox_wald_pvalue_neglog10", "cox_hazard_ratio","n")
surv_test_temp <- list()
## pick one of the following data_source:
data_source <- "WholeData"  ## pDT index as biomarker
data_source <- "MKI67"  ## KI-67 as biomarker
for(j in 1:length(survival_type_set)){
  survival_type <- survival_type_set[j]
  
  
  if(data_source == "WholeData"){  ## for whole analysis (without splitting by tissue)
    srvl <- GeneSymbolP2_df_predict
    colnames(srvl)[colnames(srvl) == "Primary_site"] <- "Tissue"
    if(survival_type == "OS") {
      srvl <- SurvivalData.generator(srvl, 2, 8, col_checkNA = 5:6)
      surv_object <- Surv(time = srvl$OS_MONTHS, event = srvl$fustat)
    } else if(survival_type == "DFS") {
      srvl <- SurvivalData.generator(srvl, 2, 10, col_checkNA = 7:8)
      surv_object <- Surv(time = srvl$DFS_MONTHS, event = srvl$fustat)
    } else if(survival_type == "DSS") {
      srvl <- SurvivalData.generator(srvl, 2, 12, col_checkNA = 1:2)
      surv_object <- Surv(time = srvl$DSS_MONTHS, event = srvl$fustat)
    } else if(survival_type == "PFS") {
      srvl <- SurvivalData.generator(srvl, 2, 14, col_checkNA = 3:4)
      surv_object <- Surv(time = srvl$PFS_MONTHS, event = srvl$fustat)
    }
    
  } else if (data_source == "MKI67"){  ## for MKI67 biomarker
    srvl <- GeneSymbolP2_df_MKI67_df
    colnames(srvl)[colnames(srvl) == "Primary_site"] <- "Tissue"
    if(survival_type == "OS") {
      srvl <- SurvivalData.generator(srvl, 2, 8, col_checkNA = 5:6, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
      surv_object <- Surv(time = srvl$OS_MONTHS, event = srvl$fustat)
    } else if(survival_type == "DFS") {
      srvl <- SurvivalData.generator(srvl, 2, 10, col_checkNA = 7:8, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
      surv_object <- Surv(time = srvl$DFS_MONTHS, event = srvl$fustat)
    } else if(survival_type == "DSS") {
      srvl <- SurvivalData.generator(srvl, 2, 12, col_checkNA = 1:2, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
      surv_object <- Surv(time = srvl$DSS_MONTHS, event = srvl$fustat)
    } else if(survival_type == "PFS") {
      srvl <- SurvivalData.generator(srvl, 2, 14, col_checkNA = 3:4, custom_fustat = T, custom_fustat_label = c("MKI67_Group","Low","High"))
      surv_object <- Surv(time = srvl$PFS_MONTHS, event = srvl$fustat)
    }
    
  } else {stop("Specify data_source.")}
  
  
  if (nrow(srvl)>0) {
    if(data_source == "WholeData"){  ## for whole analysis (without splitting by tissue)
      ### whole data
      surv_fit_PI <- survfit(surv_object ~ pDT_Group, data = srvl)
      surv_fit_sum <- summary(surv_fit_PI)
      ## p-value
      surv_diff <- survdiff(surv_object ~ pDT_Group, data = srvl)
      ## cox
      surv_cox <- coxph(surv_object ~ pDT_Group, data = srvl)
      surv_cox_sum <- summary(surv_cox)
      
    } else if (data_source == "MKI67"){  ## for MKI67 biomarker
      ### MKI67
      surv_fit_PI <- survfit(surv_object ~ MKI67_Group, data = srvl)
      surv_fit_sum <- summary(surv_fit_PI)
      # p-value
      surv_diff <- survdiff(surv_object ~ MKI67_Group, data = srvl)
      # cox
      surv_cox <- coxph(surv_object ~ MKI67_Group, data = srvl)  ## MKI67
      surv_cox_sum <- summary(surv_cox)
      
    } else {stop("Specify data_source.")}
    
    
    ## plot KM curves
    names(surv_fit_PI$strata) <- gsub("_"," ",names(surv_fit_PI$strata)) ## changing "MKI67_Group=Low KI67" into "MKI67 group=Low KI67". can also apply for other cases
    names(surv_fit_PI$strata) <- gsub("="," = ",names(surv_fit_PI$strata)) ## changing "MKI67 group=Low KI67" into "MKI67 group = Low KI67". can also apply for other cases
    surv_plot <- ggsurvplot(surv_fit_PI, data = srvl,
                            title = paste("Survival Curves (Kaplan-Meie) (",survival_type,") Based on MKI67 Level", sep = ""),  ## MKI67
                            subtitle = paste("n =", length(surv_object)),
                            xlab = "Time (Month)", ylab = "Survival Probability",
                            ncensor.plot = TRUE,
                            ggtheme = theme_light() + theme(text = element_text(size = 18, family = "Times New Roman"),
                                                            plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 17, hjust = 0.5),
                                                            axis.title = element_text(size = 18), axis.text = element_text(size = 17), legend.text = element_text(size = 17),
                                                            axis.text.y = element_text(angle = 45))+
                              theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
                                    axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm")) ## axis title 10 units from two sides; right plot margin=0.7cm
    )
    surv_plot$ncensor.plot$labels$title <- "Number of Censoring"  ## normalize the censor plot title
    surv_plot$plot <- surv_plot$plot + annotate("text",label= paste("Log-rank\n",signif(surv_diff$pvalue,digits=4),sep=""),x=max(surv_fit_PI$time)*0.75,y=0.9, size = 6, family = "Times New Roman")
    
    summary_surv_score_temp$KM_pvalue[j] <- surv_diff$pvalue  ## p-value from KM curves
    summary_surv_score_temp$KM_pvalue_neglog10[j] <- -log10(summary_surv_score_temp$KM_pvalue[j])
    summary_surv_score_temp$cox_wald_pvalue[j] <- surv_cox_sum$waldtest["pvalue"]  ## p-value from wald test
    summary_surv_score_temp$cox_wald_pvalue_neglog10[j] <- -log10(summary_surv_score_temp$cox_wald_pvalue[j])
    summary_surv_score_temp$cox_hazard_ratio[j] <- surv_cox_sum$coefficients[,"exp(coef)"]  ## exp(coef) is the hazard ratio
    
  }else {cat(survival_type, "doesn't have infor.\n");
    surv_object <- NA; surv_fit_PI <- NA; surv_fit_sum <- NA; surv_diff <- NA; surv_plot <- NA;
    surv_cox <- NA; surv_cox_sum <- NA}
  
  summary_surv_score_temp$n[j] <- nrow(srvl)
  rownames(summary_surv_score_temp)[j] <- survival_type
  
  ## data output
  surv_test_temp[[survival_type]] <- list(srvl, surv_object, surv_fit_PI, surv_fit_sum, surv_diff, surv_plot,
                                          surv_cox, surv_cox_sum)
  names(surv_test_temp[[survival_type]]) <- c("srvl", "surv_object", "surv_fit_PI", "surv_fit_sum", "surv_diff", "surv_plot",
                                              "surv_cox", "surv_cox_sum")
  rm(srvl, surv_object, surv_fit_PI, surv_fit_sum, surv_diff, surv_plot, surv_cox, surv_cox_sum)
}

## output
## whole data
summary_surv_score_all <- summary_surv_score_temp; surv_result_all <- surv_test_temp; rm(summary_surv_score_temp, surv_test_temp)
surv_result_all$OS$surv_plot$ncensor.plot <- surv_result_all$OS$surv_plot$ncensor.plot + scale_y_continuous(breaks = seq(0,12,by = 4)) ## changing the y axis manually to avoid over-crowded
surv_result_all$DFS$surv_plot$ncensor.plot <- surv_result_all$DFS$surv_plot$ncensor.plot + scale_y_continuous(breaks = seq(0,6,by = 2)) ## changing the y axis manually to avoid over-crowded
## MKI67
summary_surv_score_GENT2_MKI67 <- summary_surv_score_temp; surv_test_GENT2_MKI67 <- surv_test_temp; rm(summary_surv_score_temp, surv_test_temp)
surv_test_GENT2_MKI67$OS$surv_plot$ncensor.plot <- surv_test_GENT2_MKI67$OS$surv_plot$ncensor.plot + scale_y_continuous(breaks = c(0,4,7,11)) ## changing the y axis manually to avoid over-crowded


## plot - pick one of the following:
## for GENT2
plot_data <- surv_result_all$DFS; survival_type<- "DFS"; Tissue <- "GENT2"
## for GENT2 with KI67
plot_data <- surv_test_GENT2_MKI67$DFS; survival_type<- "DFS"; Tissue <- "GENT2"

## then run the following: (silence and un-silence part of the code)
surv_object<-plot_data$surv_object
names(plot_data$surv_fit_PI$strata) <- gsub("_"," ",names(plot_data$surv_fit_PI$strata)) ## changing "MKI67_Group=Low KI67" into "MKI67 group=Low KI67". can also apply for other cases
names(plot_data$surv_fit_PI$strata) <- gsub("="," = ",names(plot_data$surv_fit_PI$strata)) ## changing "MKI67 group=Low KI67" into "MKI67 group = Low KI67". can also apply for other cases

if(plot_data$surv_diff$pvalue<0.01){
  input_bquote <- bquote(paste(.(sub("e.*","",signif(plot_data$surv_diff$pvalue,digits=4)%>%scientific())%>%as.numeric())
                               %*% 10^.(sub(".*e","",signif(plot_data$surv_diff$pvalue,digits=4)%>%scientific())%>%as.numeric())))
} else {input_bquote <- signif(plot_data$surv_diff$pvalue,digits=4)}

surv_plot <- ggsurvplot(plot_data$surv_fit_PI, data = plot_data$srvl,
                        title = paste("Kaplan-Meier Curves (",survival_type,") on ",Tissue," Cancer", sep = ""),  ## pDT
                        # title = paste("Kaplan-Meier Curves (",survival_type,") Based on MKI67 Level", sep = ""),  ## MKI67
                        subtitle = paste("n =", length(plot_data$surv_object)),
                        xlab = "Time (Month)", ylab = "Survival Probability",
                        ncensor.plot = TRUE,
                        ggtheme = theme_classic() + theme(text = element_text(size = 19, family = "Times New Roman"),
                                                          plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 18, hjust = 0.5),
                                                          axis.title = element_text(size = 19), axis.text = element_text(size = 18), legend.text = element_text(size = 18),
                                                          axis.text.y = element_text(angle = 45))+
                          theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
                                axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm")) ## axis title 10 units from two sides; right plot margin=0.7cm
)
surv_plot$ncensor.plot$labels$title <- "Number of Censoring"  ## normalize the censor plot title
surv_plot$ncensor.plot$labels$y <- "n"  ## normalize the censor plot title
surv_plot$plot <- surv_plot$plot + 
  annotate("text",label= paste("Log-rank",sep=""),x=max(plot_data$surv_fit_PI$time)*0.75,y=0.9, size = 6, family = "Times New Roman")+
  annotate("text",label= input_bquote,x=max(plot_data$surv_fit_PI$time)*0.75,y=0.75, size = 6, family = "Times New Roman")
hold_plot <- surv_plot[1:2]
for(i in 1:length(hold_plot)){ ## remove legend from censor plot
  if(i%%2==0){hold_plot[[i]]<-hold_plot[[i]]+theme(legend.position="none")} }

## This is Figure 7A or Figure 7B ##
hold_plot
ggexport(
  marrangeGrob(hold_plot, 
               layout_matrix = layout_srvl_withcensor_2in1[,1:4]),
  filename = "plot/GENT2.pdf",
  width = 8, height = 8.5)#; rm(surv_object)


##------- therapy analysis (Freeman, 4 cohorts): data treatment and prediction -----------
## apply model
CPB4_exp_pretr <- CPB4_exp_pretr %>% left_join(gene_list_exp[,c(1,3,5)]) %>%  ## annotate the gene with the gene list from before
  relocate(hgnc_symbol,entrezgene_id, .after = "Name")
## pick genes
number_gene <- 50; model <- model_ExpDT_reg$model_ExpDT_50
ind_gene <- summary_RNA_expression[1:number_gene,c("Gene_Transcript","Entrez_Gene_id","Hugo_Entrez")]

### for the one without batch correction
## a subset containing only model genes
CPB4_exp_pretr_sub <- CPB4_exp_pretr %>% filter(hgnc_symbol %in% ind_gene$Gene_Transcript); cat("Any missing genes:",nrow(CPB4_exp_pretrB_sub)==length(ind_gene))
## annotate with model gene name
table(paste(CPB4_exp_pretr_sub$hgnc_symbol,CPB4_exp_pretr_sub$entrezgene_id,sep="_") %in% model$coefnames)  ## check whether Hugo_Entrez is the same as model
CPB4_exp_pretr_sub <- CPB4_exp_pretr_sub %>% left_join(ind_gene[,c("Gene_Transcript","Hugo_Entrez")],by=c(hgnc_symbol="Gene_Transcript")) %>% relocate("Hugo_Entrez")
## generate workable input table
CPB4_exp_pretr_sub <- CPB4_exp_pretr_sub %>% remove_rownames() %>% column_to_rownames("Hugo_Entrez") ## gene name to rowname
CPB4_exp_pretr_sub <- CPB4_exp_pretr_sub[,-c(1:3)] %>% t() %>% as.data.frame()
CPB4_exp_pretr_sub <- CPB4_exp_pretr_sub %>% relocate(ind_gene$Hugo_Entrez)
## rescale
processCenter <- preProcess(CPB4_exp_pretr_sub, method = c("scale", "center"))
CPB4_exp_pretr_sub <- predict(processCenter,CPB4_exp_pretr_sub); rm(processCenter)
## prediction
predicitons_ExpDT <- predict(model, CPB4_exp_pretr_sub)
names(predicitons_ExpDT) <- rownames(CPB4_exp_pretr_sub)  ## marked with patient id
predicitons_ExpDT <- data.frame(predicted_DT = predicitons_ExpDT)  ## convert it into data frame
predicitons_ExpDT$DEPMAPID <- rownames(predicitons_ExpDT)

## fix patient id format
SKCM_NIR_patient_list_fullcohort_sub <- SKCM_NIR_patient_list_fullcohort
SKCM_NIR_patient_list_fullcohort_sub$rna_id <- gsub("-","_",SKCM_NIR_patient_list_fullcohort_sub$rna_id) ## change all "-" into "_"
## fix patient id format
predicitons_ExpDT$DEPMAPID_fix <- predicitons_ExpDT$DEPMAPID
predicitons_ExpDT$DEPMAPID_fix <- gsub("\\.","_",predicitons_ExpDT$DEPMAPID_fix) ## change all "." into "_"
predicitons_ExpDT$DEPMAPID_fix <- gsub("-","_",predicitons_ExpDT$DEPMAPID_fix) ## change all "-" into "_"

predicitons_ExpDT <- predicitons_ExpDT %>% left_join(SKCM_NIR_patient_list_fullcohort_sub[,c("rna_id","cohort","clinical_rna_id")], by=c(DEPMAPID_fix="rna_id"))

## a separated table to keep the data for each cohort
predicitons_ExpDT_cohortSplit <- predicitons_ExpDT %>% split(predicitons_ExpDT$cohort)
names(predicitons_ExpDT_cohortSplit)[4] <- "Allen"

for(i in 1:length(predicitons_ExpDT_cohortSplit)){
  current_cohort <- names(predicitons_ExpDT_cohortSplit)[i]
  predicitons_ExpDT_cohortSplit[[i]] <- predicitons_ExpDT_cohortSplit[[i]] %>% 
    left_join(SKCM_NIR_patient_4cohort[[current_cohort]][,c(5:12)], by = c(clinical_rna_id="Sample_id"))
}
predicitons_ExpDT_cohortSplit$Riaz$Response[predicitons_ExpDT_cohortSplit$Riaz$Response=="NA"] <- NA

## final result
predicitons_ExpDT_v2 <- rbind(predicitons_ExpDT_cohortSplit[[1]],
                              predicitons_ExpDT_cohortSplit[[2]],
                              predicitons_ExpDT_cohortSplit[[3]],
                              predicitons_ExpDT_cohortSplit[[4]]
)

## output
predicitons_ExpDT_CPB4 <- predicitons_ExpDT
predicitons_ExpDT_v2_CPB4 <- predicitons_ExpDT_v2
predicitons_ExpDT_cohortSplit_CPB4 <- predicitons_ExpDT_cohortSplit


## extract Riaz data
Riaz_pre_t <- predicitons_ExpDT_v2_CPB4 %>% filter(cohort=="Riaz")

## self defined labels
label_1 <- c("With Prior CTLA4","Without Prior CTLA4")
names(label_1) <- c("TRUE","FALSE")
label_2 <- c("Post Anti-PD1","Pre Anti-PD1")
names(label_2) <- c("Post","Pre")

## This is Figure 8C ##
ggplot(Riaz_pre_t %>% filter(!is.na(Response)), aes(Response,predicted_DT)) + geom_violin() +
  geom_jitter(shape=1) +
  geom_boxplot(width=0.1,colour="black",outlier.colour = NA)+
  stat_compare_means(family = "Times New Roman",size=5, label.y=127, label.x=1.5, )+
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5)) +
  guides(color="none")+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = "Response", y = "pDT hr", col = "none", shape = "Dataset",
       title = "pDT on Therapy Response from Riaz Cohort") +
  facet_grid(c("Prior_CTLA4","Timing"), labeller = labeller(Prior_CTLA4=label_1,Timing=label_2)) + labs(subtitle = "Anti-PD1")


## MKI67 level for 4 cohorts from Freeman ##
## make a single data frame for all patient data
SKCM_NIR_patient_4cohort_df <- rbind(SKCM_NIR_patient_4cohort[[1]][,-13],
                                     SKCM_NIR_patient_4cohort[[2]],
                                     SKCM_NIR_patient_4cohort[[3]],
                                     SKCM_NIR_patient_4cohort[[4]])

working_table <- CPB4_exp_pretr

working_table <- working_table[which(working_table$hgnc_symbol=="MKI67"),]  ## RNA-seq

working_table <- working_table[,-1:-3] %>% t() %>% as.data.frame()
colnames(working_table) <- "gene"
working_table$Sample_id <- gsub("-","_",rownames(working_table))

working_table <- working_table %>% left_join(SKCM_NIR_patient_list_fullcohort_sub[,c("rna_id","cohort","clinical_rna_id")], 
                                             by=c(Sample_id="rna_id")) %>% 
  left_join(SKCM_NIR_patient_4cohort_df[,c(5:12)], by = c(clinical_rna_id="Sample_id"))
working_table$Response[working_table$Response=="NA"] <- NA

## output
CPB4_exp_pretr_MKI67 <- working_table; rm(working_table)

##------- NR/R ratio data for volcano plot for 4 cohorts ------
## pick one of the following to construct the working_table_collection: pDT, or MKI67
## pDT ratio between R and NR
working_table_collection <- predicitons_ExpDT_cohortSplit_CPB4

working_table_collection$Riaz <- working_table_collection$Riaz %>% filter(Overall_survival!="NA")
working_table_collection$`Riaz-Pre` <- working_table_collection$Riaz[grep("_Pre_",working_table_collection$Riaz$DEPMAPID_fix),]
working_table_collection$`Riaz-Pre`$cohort <- "Riaz-Pre"
working_table_collection$`Riaz-Post` <- working_table_collection$Riaz[grep("_On_",working_table_collection$Riaz$DEPMAPID_fix),]
working_table_collection$`Riaz-Post`$cohort <- "Riaz-Post"
working_table_collection$`MGH-Pre` <- working_table_collection$MGH[working_table_collection$MGH$Timing=="Pre",]
working_table_collection$`MGH-Pre`$cohort <- "MGH-Pre"
working_table_collection$`MGH-Post` <- working_table_collection$MGH[working_table_collection$MGH$Timing=="Post",]
working_table_collection$`MGH-Post`$cohort <- "MGH-Post"

## MKI67 ratio between R and NR
working_table_collection <- CPB4_exp_pretr_MKI67 %>% split(CPB4_exp_pretr_MKI67$cohort)
names(working_table_collection)[4] <- "Allen"

working_table_collection$Riaz <- working_table_collection$Riaz %>% filter(Overall_survival!="NA")
working_table_collection$`Riaz-Pre` <- working_table_collection$Riaz[grep("_Pre_",working_table_collection$Riaz$Sample_id),]
working_table_collection$`Riaz-Pre`$cohort <- "Riaz-Pre"
working_table_collection$`Riaz-Post` <- working_table_collection$Riaz[grep("_On_",working_table_collection$Riaz$Sample_id),]
working_table_collection$`Riaz-Post`$cohort <- "Riaz-Post"
working_table_collection$`MGH-Pre` <- working_table_collection$MGH[working_table_collection$MGH$Timing=="Pre",]
working_table_collection$`MGH-Pre`$cohort <- "MGH-Pre"
working_table_collection$`MGH-Post` <- working_table_collection$MGH[working_table_collection$MGH$Timing=="Post",]
working_table_collection$`MGH-Post`$cohort <- "MGH-Post"

## then run the following:
result <- data.frame(matrix(NA,1,13))
colnames(result) <- c("cohort","mean_NR","mean_R","median_NR","median_R",
                      "mean_ratio","median_ratio","wilcox_p","n_NR","n_R","types","drugs","targets")
for(i in 1:length(working_table_collection)){
  working_table <- working_table_collection[[i]]
  working_table_NR <- working_table %>% filter(Response=="NR")
  working_table_R <- working_table %>% filter(Response=="R")
  ## cohort
  result[i,"cohort"] <- names(working_table_collection)[i]
  ## mean
  result[i,"mean_NR"] <- working_table_NR[,1] %>% mean()
  result[i,"mean_R"] <- working_table_R[,1] %>% mean()
  ## median
  result[i,"median_NR"] <- working_table_NR[,1] %>% median()
  result[i,"median_R"] <- working_table_R[,1] %>% median()
  
  ## n
  result[i,"n_NR"] <- working_table_NR %>% nrow()
  result[i,"n_R"] <- working_table_R %>% nrow()
  
  ## ratio
  result[i,"mean_ratio"] <- result$mean_NR[i]/result$mean_R[i]
  result[i,"median_ratio"] <- result$median_NR[i]/result$median_R[i]
  ## wilcox p-value
  result[i,"wilcox_p"] <- wilcox.test(working_table_NR[,1], working_table_R[,1])$p.value
  
  result[i,"types"] <- paste(unique(working_table$types), collapse = ",")
  result[i,"drugs"] <- paste(unique(working_table$drugs), collapse = ",")
  result[i,"targets"] <- paste(unique(working_table$targets), collapse = ",")
}; rm(i, working_table_collection,working_table_NR,working_table_R, working_table)
## output: pick one of the following
## pDT
summary_response_pDTratio_CPB4 <- result; rm(result)
## MKI67
summary_response_pDTratio_CPB4_MKI67 <- result; rm(result)


##------- therapy analysis (Lee, Immunotherapy, targeted drug therapy): data treatment and prediction -----------
## apply model and make prediction (pDT index)
## pick one of the following:
## targeted drug
model_collection <- model_ExpDT_targeteddrug_collection; working_table_collection <- Targeted_10; list_dataset <- names(model_collection); list_gene_collection <- list_gene_targeteddrug_10_DT; list_patient_collection<-list_patient_targeteddrug_collection
## immuno
model_collection <- model_ExpDT_immuno_collection; working_table_collection <- Immuno_13_opt; list_dataset <- names(model_collection); list_gene_collection <- list_gene_immuno_13; list_patient_collection <- list_patient_immuno_collection

## then run this: (no code silence is needed)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")  ## for gene length in TPM
predicitons_ExpDT_collection <- c()
for(i in 1:length(list_dataset)){
  if(i==1){result_test <- as.data.frame(matrix(NA, length(list_dataset),2)); colnames(result_test)<-c("Wilcox","t.test")}
  current_dataset <- list_dataset[i]
  model <- model_collection[[current_dataset]]$model_ExpDT_50
  
  working_table <- working_table_collection[[current_dataset]]$mRNA
  rownames(working_table) <- gsub("\\.","-",rownames(working_table))
  if(current_dataset=="VanAllen"){
    colnames(working_table)=sub("MEL-IPI_","",colnames(working_table))
    colnames(working_table)=sub("-Tumor-SM.*","",colnames(working_table))
  }
  ## normalized by TPM
  if(current_dataset %in% c("Cho","Riaz")){
    ## calculate gene length
    gene_info <- getBM(attributes = c("hgnc_symbol", "start_position","end_position"),
                       filters = 'hgnc_symbol',
                       values = rownames(working_table),
                       mart = ensembl)
    gene_info$gene_length <- gene_info$end_position - gene_info$start_position
    working_table$gene_length <- gene_info$gene_length[match(rownames(working_table),gene_info$hgnc_symbol)]
    ## TPM
    working_table <- apply(working_table[,c(-ncol(working_table))], 2,function(X) X/(working_table$gene_length/1000))
    working_table <- apply(working_table, 2, function(X) X/sum(as.numeric(na.omit(X)))*10^6) %>% as.data.frame()
    working_table <- log2(working_table+1) ## Pseudo log2
  }
  
  list_gene <- list_gene_collection[[current_dataset]]
  ## in case Hugo_Entrez is missing which is the format of predictor name used by model:
  # list_gene <- list_gene %>% left_join(summary_RNA_expression[,c("Gene_Transcript","Hugo_Entrez")],by=c(Name="Gene_Transcript"))
  ind_gene <- match(sub("_.*","",gsub("`","",model$coefnames)), rownames(working_table))
  ## in case furthur process on gene name is needed:
  # ind_gene <- match(sub("_.*","",gsub("`","",model$coefnames)), gsub("\\.","-",rownames(working_table)))
  if(length(ind_gene)!=length(model$coefnames)){stop("gene is missing")}
  working_table <- working_table[ind_gene,] %>% t() %>% as.data.frame() 
  
  ## apply model
  ## rescale
  processCenter <- preProcess(working_table, method = c("scale", "center"))
  working_table <- predict(processCenter,working_table); rm(processCenter)
  ## prediction
  predicitons_ExpDT <- predict(model, working_table)
  names(predicitons_ExpDT) <- rownames(working_table)  ## marked with patient id
  predicitons_ExpDT <- data.frame(predicted_DT = predicitons_ExpDT)  ## convert it into data frame
  predicitons_ExpDT$DEPMAPID <- rownames(predicitons_ExpDT)
  
  predicitons_ExpDT <- predicitons_ExpDT %>% left_join(list_patient_collection[[current_dataset]], by=c(DEPMAPID="Patient_id"))
  
  if(current_dataset %in% c("Liu","Riaz")){
    predicitons_ExpDT$Response2 <- predicitons_ExpDT$Response
    ind <- grep("R",list_patient_collection[[current_dataset]][["Response"]])
    indd <- grep("D",list_patient_collection[[current_dataset]][["Response"]])
    predicitons_ExpDT$Response[ind] <- "Responder"
    predicitons_ExpDT$Response[indd] <- "Non-responder"
  }
  predicitons_ExpDT_collection[[current_dataset]] <- predicitons_ExpDT
  
  
  result_test$Wilcox[i] <- wilcox.test(na.omit(predicitons_ExpDT$predicted_DT[predicitons_ExpDT$Response=="Responder"]),
                                       na.omit(predicitons_ExpDT$predicted_DT[predicitons_ExpDT$Response=="Non-responder"]))$p.value
  result_test$t.test[i] <- t.test(na.omit(predicitons_ExpDT$predicted_DT[predicitons_ExpDT$Response=="Responder"]),
                                  na.omit(predicitons_ExpDT$predicted_DT[predicitons_ExpDT$Response=="Non-responder"]))$p.value
  rownames(result_test)[i] <- current_dataset
}
result_test$Wilcox.FDR <- p.adjust(result_test$Wilcox, method = "BH")
result_test$t.test.FDR <- p.adjust(result_test$t.test, method = "BH")

## output - pick one of the following, then run the specific annotation:
## targeted drug
predicitons_ExpDT_targeteddrug_collection <- predicitons_ExpDT_collection; rm(predicitons_ExpDT_collection)
result_test_targeteddrug <- result_test; rm(result_test)
## further annotation
predicitons_ExpDT_targeteddrug_collection$GSE109211$types <- Targeted_10$GSE109211$types; predicitons_ExpDT_targeteddrug_collection$GSE109211$drugs <-Targeted_10$GSE109211$drugs; predicitons_ExpDT_targeteddrug_collection$GSE109211$targets <- paste(Targeted_10$GSE109211$targets,collapse=", ")
predicitons_ExpDT_targeteddrug_collection$GSE16391$types <- Targeted_10$GSE16391$types; predicitons_ExpDT_targeteddrug_collection$GSE16391$drugs <-Targeted_10$GSE16391$drugs; predicitons_ExpDT_targeteddrug_collection$GSE16391$targets <- paste(Targeted_10$GSE16391$targets,collapse=", ")
predicitons_ExpDT_targeteddrug_collection$GSE50509$types <- Targeted_10$GSE50509$types; predicitons_ExpDT_targeteddrug_collection$GSE50509$drugs <-Targeted_10$GSE50509$drugs; predicitons_ExpDT_targeteddrug_collection$GSE50509$targets <- paste(Targeted_10$GSE50509$targets,collapse=", ")
predicitons_ExpDT_targeteddrug_collection$GSE65185$types <- Targeted_10$GSE65185$types; predicitons_ExpDT_targeteddrug_collection$GSE65185$drugs <-Targeted_10$GSE65185$drugs; predicitons_ExpDT_targeteddrug_collection$GSE65185$targets <- paste(Targeted_10$GSE65185$targets,collapse=", ")
predicitons_ExpDT_targeteddrug_collection$GSE66399$types <- Targeted_10$GSE66399$types; predicitons_ExpDT_targeteddrug_collection$GSE66399$drugs <-Targeted_10$GSE66399$drugs; predicitons_ExpDT_targeteddrug_collection$GSE66399$targets <- paste(Targeted_10$GSE66399$targets,collapse=", ")
predicitons_ExpDT_targeteddrug_collection$GSE68871$types <- Targeted_10$GSE68871$types; predicitons_ExpDT_targeteddrug_collection$GSE68871$drugs <-Targeted_10$GSE68871$drugs; predicitons_ExpDT_targeteddrug_collection$GSE68871$targets <- paste(Targeted_10$GSE68871$targets,collapse=", ")
predicitons_ExpDT_targeteddrug_collection$GSE99898$types <- Targeted_10$GSE99898$types; predicitons_ExpDT_targeteddrug_collection$GSE99898$drugs <-Targeted_10$GSE99898$drugs; predicitons_ExpDT_targeteddrug_collection$GSE99898$targets <- paste(Targeted_10$GSE99898$targets,collapse=", ")

#### immunotherapy
predicitons_ExpDT_immuno_collection <- predicitons_ExpDT_collection; rm(predicitons_ExpDT_collection)
result_test_immuno <- result_test; rm(result_test)
## further annotation
predicitons_ExpDT_immuno_collection$Liu <- predicitons_ExpDT_immuno_collection$Liu %>% left_join(Immuno_13_opt$Liu$clin[,c("samples","Primary_Type")], by=c(DEPMAPID="samples"))
predicitons_ExpDT_immuno_collection$Cho$types <- "NSCLC"; predicitons_ExpDT_immuno_collection$Cho$drugs <-"Anti-PD1"; predicitons_ExpDT_immuno_collection$Cho$targets <- "PD1"
predicitons_ExpDT_immuno_collection$Gide$types <- NA; predicitons_ExpDT_immuno_collection$Gide$drugs <- NA; predicitons_ExpDT_immuno_collection$Gide$targets <- NA
predicitons_ExpDT_immuno_collection$Liu$types <- "SKCM"; predicitons_ExpDT_immuno_collection$Liu$drugs <-"Anti-PD1"; predicitons_ExpDT_immuno_collection$Liu$targets <- "PD1"
predicitons_ExpDT_immuno_collection$Miao$types <- "KIRC"
## it is assumed that the order/id is the same, so make sure the order is correct before doing this
Immuno_13_opt$Miao[["samples"]]%>% as.numeric() == lapply(predicitons_ExpDT_immuno_collection$Miao$DEPMAPID %>% str_extract_all("[:digit:]"), function(X) paste(X, collapse = "")) %>% unlist() %>% as.numeric()
predicitons_ExpDT_immuno_collection$Miao$drugs <- Immuno_13_opt$Miao$drug ## it is assumed that the order/id is the same, so make sure the order is correct before doing this
predicitons_ExpDT_immuno_collection$Miao$targets <- "PD1"
predicitons_ExpDT_immuno_collection$Nathanson_pre$types <- "SKCM"; predicitons_ExpDT_immuno_collection$Nathanson_pre$drugs <-"Anti-CTLA4"; predicitons_ExpDT_immuno_collection$Nathanson_pre$targets <- "CTLA4"
predicitons_ExpDT_immuno_collection$Riaz$types <- NA; predicitons_ExpDT_immuno_collection$Riaz$drugs <- NA; predicitons_ExpDT_immuno_collection$Riaz$targets <- NA
predicitons_ExpDT_immuno_collection$VanAllen$types <- NA; predicitons_ExpDT_immuno_collection$VanAllen$drugs <- NA; predicitons_ExpDT_immuno_collection$VanAllen$targets <- NA


## MKI67
## pick one of the following:
## targeted drug
working_table_collection <- Targeted_10
list_patient_collection <- list_patient_targeteddrug_collection
## Immuno
working_table_collection <- Immuno_13_opt
list_patient_collection <- list_patient_immuno_collection

## then run this: (no code silence is needed)
gene = "MKI67"
output <- list()
for(i in 1:length(working_table_collection)){
  current_dataset <- names(working_table_collection)[i]
  working_clinic_df <- list_patient_collection[[current_dataset]]
  
  if(current_dataset=="Chen"){output[[current_dataset]] <- NA; next}
  working_table <- working_table_collection[[current_dataset]]$mRNA
  working_table <- working_table[rownames(working_table)==gene,]
  rownames(working_table) <- "gene"
  
  if(current_dataset=="VanAllen"){
    colnames(working_table)=sub("MEL-IPI_","",colnames(working_table))
    colnames(working_table)=sub("-Tumor-SM.*","",colnames(working_table))
  }
  
  working_table <- working_table[,-1] %>% t() %>% as.data.frame() %>% filter(gene!="NaN",gene!="NA",!is.na(gene),)
  working_table$gene <- as.numeric(working_table$gene)
  working_table$Sample_id <- rownames(working_table)
  
  
  working_table <- working_table %>% left_join(working_clinic_df, by=c(Sample_id="Patient_id"))
  if(current_dataset %in% c("Liu","Riaz")){
    working_table$Response2 <- working_table$Response
    ## in case external patient information is needed:
    # ind <- grep("R",list_patient_collection[[current_dataset]][["Response"]])
    # indd <- grep("D",list_patient_collection[[current_dataset]][["Response"]])
    ind <- grep("R",working_table$Response2)
    indd <- grep("D",working_table$Response2)
    working_table$Response[ind] <- "Responder"
    working_table$Response[indd] <- "Non-responder"
    rm(ind,indd)
  }
  
  output[[current_dataset]] <- working_table; rm(working_table)
}; rm(i, current_dataset, working_clinic_df, working_table_collection)
View(output)
## output - pick one of the following, then run the specific annotation:
Targeted_10_data_MKI67 <- output; rm(output)
Immuno_13_opt_data_MKI67 <- output; rm(output)

## further annotate
Targeted_10_data_MKI67$GSE109211$types <- Targeted_10$GSE109211$types; Targeted_10_data_MKI67$GSE109211$drugs <-Targeted_10$GSE109211$drugs; Targeted_10_data_MKI67$GSE109211$targets <- paste(Targeted_10$GSE109211$targets,collapse=", ")
Targeted_10_data_MKI67$GSE16391$types <- Targeted_10$GSE16391$types; Targeted_10_data_MKI67$GSE16391$drugs <-Targeted_10$GSE16391$drugs; Targeted_10_data_MKI67$GSE16391$targets <- paste(Targeted_10$GSE16391$targets,collapse=", ")
Targeted_10_data_MKI67$GSE50509$types <- Targeted_10$GSE50509$types; Targeted_10_data_MKI67$GSE50509$drugs <-Targeted_10$GSE50509$drugs; Targeted_10_data_MKI67$GSE50509$targets <- paste(Targeted_10$GSE50509$targets,collapse=", ")
Targeted_10_data_MKI67$GSE65185$types <- Targeted_10$GSE65185$types; Targeted_10_data_MKI67$GSE65185$drugs <-Targeted_10$GSE65185$drugs; Targeted_10_data_MKI67$GSE65185$targets <- paste(Targeted_10$GSE65185$targets,collapse=", ")
Targeted_10_data_MKI67$GSE66399$types <- Targeted_10$GSE66399$types; Targeted_10_data_MKI67$GSE66399$drugs <-Targeted_10$GSE66399$drugs; Targeted_10_data_MKI67$GSE66399$targets <- paste(Targeted_10$GSE66399$targets,collapse=", ")
Targeted_10_data_MKI67$GSE68871$types <- Targeted_10$GSE68871$types; Targeted_10_data_MKI67$GSE68871$drugs <-Targeted_10$GSE68871$drugs; Targeted_10_data_MKI67$GSE68871$targets <- paste(Targeted_10$GSE68871$targets,collapse=", ")
Targeted_10_data_MKI67$GSE99898$types <- Targeted_10$GSE99898$types; Targeted_10_data_MKI67$GSE99898$drugs <-Targeted_10$GSE99898$drugs; Targeted_10_data_MKI67$GSE99898$targets <- paste(Targeted_10$GSE99898$targets,collapse=", ")


Immuno_13_opt_data_MKI67$Liu <- Immuno_13_opt_data_MKI67$Liu %>% left_join(Immuno_13_opt$Liu$clin[,c("samples","Primary_Type")], by=c(Sample_id="samples"))
Immuno_13_opt_data_MKI67$Cho$types <- "NSCLC"; Immuno_13_opt_data_MKI67$Cho$drugs <-"Anti-PD1"; Immuno_13_opt_data_MKI67$Cho$targets <- "PD1"
Immuno_13_opt_data_MKI67$Gide$types <- NA; Immuno_13_opt_data_MKI67$Gide$drugs <- NA; Immuno_13_opt_data_MKI67$Gide$targets <- NA
Immuno_13_opt_data_MKI67$Liu$types <- "SKCM"; Immuno_13_opt_data_MKI67$Liu$drugs <-"Anti-PD1"; Immuno_13_opt_data_MKI67$Liu$targets <- "PD1"
Immuno_13_opt_data_MKI67$Miao$types <- "KIRC"
## it is assumed that the order/id is the same, so make sure the order is correct before doing this
Immuno_13_opt$Miao[["samples"]][-1]%>% as.numeric() == lapply(Immuno_13_opt_data_MKI67[["Miao"]][["Sample_id"]] %>% str_extract_all("[:digit:]"), function(X) paste(X, collapse = "")) %>% unlist() %>% as.numeric()
Immuno_13_opt_data_MKI67$Miao$drugs <- Immuno_13_opt$Miao$drug[-1]
Immuno_13_opt_data_MKI67$Miao$targets <- "PD1"
Immuno_13_opt_data_MKI67$Nathanson_pre$types <- "SKCM"; Immuno_13_opt_data_MKI67$Nathanson_pre$drugs <-"Anti-CTLA4"; Immuno_13_opt_data_MKI67$Nathanson_pre$targets <- "CTLA4"
Immuno_13_opt_data_MKI67$Riaz$types <- NA; Immuno_13_opt_data_MKI67$Riaz$drugs <- NA; Immuno_13_opt_data_MKI67$Riaz$targets <- NA
Immuno_13_opt_data_MKI67$VanAllen$types <- NA; Immuno_13_opt_data_MKI67$VanAllen$drugs <- NA; Immuno_13_opt_data_MKI67$VanAllen$targets <- NA


##------- NR/R ratio data for volcano plot for Lee ------
## pick one of the following:
## pDT ratio between R and NR
working_table_collection <- predicitons_ExpDT_targeteddrug_collection
working_table_collection <- predicitons_ExpDT_immuno_collection

## MKI67
## GSE32603, GSE3964, GSE50509, GSE8465 are in diff. Expression level - exclude
working_table_collection <- Targeted_10_data_MKI67[-which(names(Targeted_10_data_MKI67)%in%c("GSE32603", "GSE3964", "GSE50509", "GSE8465"))]
## Nathanson_pre is in diff. Expression level - exclude
working_table_collection <- Immuno_13_opt_data_MKI67[-1]
working_table_collection <- working_table_collection[-which(names(working_table_collection)%in%c("Nathanson_pre"))]

## then run the following: (code silence is not needed)
result <- data.frame(matrix(NA,1,13))
colnames(result) <- c("cohort","mean_NR","mean_R","median_NR","median_R",
                      "mean_ratio","median_ratio","wilcox_p","n_NR","n_R","types","drugs","targets")
for(i in 1:length(working_table_collection)){
  working_table <- working_table_collection[[i]]
  if(any(colnames(working_table)%in%"gene")){colnames(working_table)[colnames(working_table)=="gene"]="predicted_DT"}
  working_table_NR <- working_table %>% filter(Response=="Non-responder")
  working_table_R <- working_table %>% filter(Response=="Responder")
  ## cohort
  result[i,"cohort"] <- names(working_table_collection)[i]
  ## mean
  result[i,"mean_NR"] <- working_table_NR$predicted_DT %>% mean()
  result[i,"mean_R"] <- working_table_R$predicted_DT %>% mean()
  ## median
  result[i,"median_NR"] <- working_table_NR$predicted_DT %>% median()
  result[i,"median_R"] <- working_table_R$predicted_DT %>% median()
  
  ## n
  result[i,"n_NR"] <- working_table_NR %>% nrow()
  result[i,"n_R"] <- working_table_R %>% nrow()
  
  ## ratio
  result[i,"mean_ratio"] <- result$mean_NR[i]/result$mean_R[i]
  result[i,"median_ratio"] <- result$median_NR[i]/result$median_R[i]
  # wilcox p-value
  if(try(wilcox.test(working_table_NR$predicted_DT, working_table_R$predicted_DT)$p.value,T) %>% grep(pattern="Error") %>% length() != 1) {
    result[i,"wilcox_p"] <- wilcox.test(working_table_NR$predicted_DT, working_table_R$predicted_DT)$p.value
  } else {result[i,"wilcox_p"] <- NA}
  
  result[i,"types"] <- paste(unique(working_table$types), collapse = " | ")
  result[i,"drugs"] <- paste(unique(working_table$drugs), collapse = " | ")
  result[i,"targets"] <- paste(unique(working_table$targets), collapse = " | ")
}; #rm(i, working_table_collection,working_table_NR,working_table_R, working_table)
## output - pick one of the following:
summary_response_pDTratio_targeteddrug <- result %>% filter(!is.na(wilcox_p)); rm(result)
summary_response_pDTratio_immuno <- result %>% filter(!is.na(wilcox_p)); rm(result)

summary_response_pDTratio_targeteddrug_MKI67 <- result %>% filter(!is.na(wilcox_p)); rm(result)
summary_response_pDTratio_immuno_MKI67 <- result %>% filter(!is.na(wilcox_p)); rm(result)

##------- therapy analysis (CTR database): data treatment and prediction -----------
## for single data frame method:
## apply model - pick one of the followings:
## RNA_seq
model <- model_ExpDT_reg$model_ExpDT_50
working_table = RNAseq_data_tpm_df
working_clinic_df <- RNAseq_clinic_df 
## microarray
## MC20161
model <- model_ExpDT_ctr_MC20161$model_ExpDT_50
working_table = MC20161_data_tpm_df
working_clinic_df <- MC_clinic_df 
## MC12399
model <- model_ExpDT_ctr_MC12399$model_ExpDT_50
working_table = MC12399_data_tpm_df
working_clinic_df <- MC_clinic_df 

## then run the following shared code:
rownames(working_table) <- working_table$X
working_table <- working_table[,-1]

ind_gene <- match(sub("_.*","",gsub("`","",model$coefnames)), rownames(working_table))
working_table <- working_table[ind_gene,] %>% t() %>% as.data.frame() 
colnames(working_table) <- gsub("`","",model$coefnames)[match(colnames(working_table), sub("_.*","",gsub("`","",model$coefnames)))]  ## change into model gene name

## rescale
processCenter <- preProcess(working_table, method = c("scale", "center"))
working_table <- predict(processCenter,working_table); rm(processCenter)
## prediction
predicitons_ExpDT <- predict(model, working_table)
names(predicitons_ExpDT) <- rownames(working_table)  ## marked with patient id
predicitons_ExpDT <- data.frame(predicted_DT = predicitons_ExpDT)  ## convert it into data frame
predicitons_ExpDT$DEPMAPID <- rownames(predicitons_ExpDT)

## annotation
predicitons_ExpDT <- predicitons_ExpDT %>% left_join(working_clinic_df[,c("Sample_id","Cancer_type_level1","Original_response","Response","Response_standard",
                                                                          "Drug_list", "Additional_information_about_the_therapeutic_regimen","CTR_cohort_id")], 
                                                     by=c(DEPMAPID="Sample_id"))

## output - pick the followings:
predicitons_ExpDT_RNASeq <- predicitons_ExpDT; rm(predicitons_ExpDT)
predicitons_ExpDT_MC20161 <- predicitons_ExpDT; rm(predicitons_ExpDT)
predicitons_ExpDT_MC12399 <- predicitons_ExpDT; rm(predicitons_ExpDT)


## use MKI67 level instead of pDT index from the model

## extracting MKI67 expression level
## pick one of the following:
## RNA-seq
working_table <- RNAseq_data_tpm_df
working_clinic_df <- RNAseq_clinic_df
working_table <- working_table[rownames(working_table)=="MKI67",]  ## RNA-seq
## MC20161
working_table <- MC20161_data_tpm_df
working_clinic_df <- MC_clinic_df
working_table <- working_table[working_table$X=="MKI67",]  ## microarray
## MC12399
working_table <- MC12399_data_tpm_df
working_clinic_df <- MC_clinic_df
working_table <- working_table[working_table$X=="MKI67",]  ## microarray

## then run the following:
rownames(working_table) <- "gene" ## consistent name
working_table <- working_table[,-1] %>% t() %>% as.data.frame()
working_table$Sample_id <- rownames(working_table)

working_table <- working_table %>% left_join(working_clinic_df[,c("Sample_id","Cancer_type_level1","Original_response","Response","Response_standard",
                                                                  "Drug_list", "Additional_information_about_the_therapeutic_regimen","CTR_cohort_id")], 
                                             by=c(Sample_id="Sample_id"))

## output - pick one of the following:
## RNA-seq
RNAseq_data_MKI67 <- working_table
## MC20161
MC20161_data_MKI67 <- working_table
## MC12399
MC12399_data_MKI67 <- working_table


## single cohort plot
## This is Figure 8A ##
plot_table <- predicitons_ExpDT_MC12399 %>% filter(Drug_list=="Anthracycline+Taxane")
## wilcox q-value annotation
numb <- summary_MC12399_drug_cancer$wilcox_qvalue[summary_MC12399_drug_cancer$drug == "Anthracycline+Taxane"]
input_bquote <- bquote(paste('Wilcox q-value =' ~ .(sub("e.*","",signif(numb,digits=4))%>%as.numeric())
                             %*% 10^.(sub(".*e","",signif(numb,digits=4))%>%as.numeric())))

ggplot(plot_table, aes(Response, predicted_DT)) + geom_violin(aes(fill = Response)) + geom_boxplot(aes(fill = Response),width=0.1, alpha = 0.8) + geom_jitter(shape = 1, alpha = 0.5, width = 0.2)+
  theme_classic() + 
  scale_x_discrete(labels=c(Non_response="Nonresponder", Response="Responder")) +
  scale_y_continuous(limits = c(min(plot_table$predicted_DT), max(plot_table$predicted_DT)*1.05))+
  theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5), 
        axis.text = element_text(size = 20), legend.position = "none") +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  annotate("text",family = "Times New Roman", x=1.8, y=133, size = 6, 
           label = input_bquote) +  ## qvalue 
  labs(x="",y = "pDT hr", fill = "Response", title = paste("pDT on Therapy Response on Breast Cancer"), subtitle = "Anthracycline + Taxane")


##------- NR/R ratio data for volcano plot for CTR database ------
## pick one of the following:
## pDT ratio between R and NR
working_table_collection <- predicitons_ExpDT_RNASeq; current_DS <- "RNASeq"
working_table_collection <- predicitons_ExpDT_MC20161; current_DS <- "MC20161"
working_table_collection <- predicitons_ExpDT_MC12399; current_DS <- "MC12399"
## MKI67 ratio between R and NR
working_table_collection <- RNAseq_data_MKI67; current_DS <- "RNASeq"
working_table_collection <- MC20161_data_MKI67; current_DS <- "MC20161"
working_table_collection <- MC12399_data_MKI67; current_DS <- "MC12399"

## then run the following:
colnames(working_table_collection)[colnames(working_table_collection)=="gene"] = "predicted_DT" ## have a consistent template

result <- data.frame(matrix(NA,1,15))
colnames(result) <- c("cohort","mean_NR","mean_R","median_NR","median_R",
                      "mean_ratio","median_ratio","wilcox_p","n_NR","n_R","types","drugs","targets",
                      "cancer","drug")
## pick one of the following within this shared code:
## split samples according to cancer type
working_table_collection <- working_table_collection %>% split(working_table_collection$Cancer_type_level1); split_type = "cancer"
## split samples according to drug
working_table_collection <- working_table_collection %>% split(working_table_collection$Drug_list); split_type = "drug"
## split samples according to drug and cancer
working_table_collection <- working_table_collection %>% split(paste(working_table_collection$Drug_list,working_table_collection$Cancer_type_level1,sep=" on ")); split_type = "drug_cancer"
## then run the loop for each split type
for(i in 1:length(working_table_collection)){
  working_table <- working_table_collection[[i]]
  if(any(colnames(working_table)%in%"gene")){colnames(working_table)[colnames(working_table)=="gene"]="predicted_DT"}
  working_table_NR <- working_table %>% filter(Response=="Non_response")
  working_table_R <- working_table %>% filter(Response=="Response")
  ## cohort
  result[i,"cohort"] <- paste(names(working_table_collection)[i], current_DS, sep="_")
  ## mean
  result[i,"mean_NR"] <- working_table_NR$predicted_DT %>% mean()
  result[i,"mean_R"] <- working_table_R$predicted_DT %>% mean()
  ## median
  result[i,"median_NR"] <- working_table_NR$predicted_DT %>% median()
  result[i,"median_R"] <- working_table_R$predicted_DT %>% median()
  
  ## n
  result[i,"n_NR"] <- working_table_NR %>% nrow()
  result[i,"n_R"] <- working_table_R %>% nrow()
  
  ## ratio
  result[i,"mean_ratio"] <- result$mean_NR[i]/result$mean_R[i]
  result[i,"median_ratio"] <- result$median_NR[i]/result$median_R[i]
  ## wilcox p-value
  if(try(wilcox.test(working_table_NR$predicted_DT, working_table_R$predicted_DT)$p.value,T) %>% grep(pattern="Error") %>% length() != 1) {
    result[i,"wilcox_p"] <- wilcox.test(working_table_NR$predicted_DT, working_table_R$predicted_DT)$p.value
  } else {result[i,"wilcox_p"] <- NA}
  
  result[i,"types"] <- paste(unique(working_table$types), collapse = ",")
  result[i,"drugs"] <- paste(unique(working_table$drugs), collapse = ",")
  result[i,"targets"] <- paste(unique(working_table$targets), collapse = ",")
  
  if(result[i,"mean_NR"]=="NaN"){result[i,"mean_NR"]=NA}
  if(result[i,"mean_R"]=="NaN"){result[i,"mean_R"]=NA}
  if(result[i,"mean_ratio"]=="NaN"){result[i,"mean_ratio"]=NA}
  
  if(split_type=="cancer"){
    result[i,"cancer"] <- unique(working_table$Cancer_type_level1)
  } else if(split_type=="drug"){
    result[i,"drug"] <- unique(working_table$Drug_list)
  } else if(split_type=="drug_cancer"){
    result[i,"cancer"] <- unique(working_table$Cancer_type_level1)
    result[i,"drug"] <- unique(working_table$Drug_list)
  }
  
}; rm(i, working_table_collection,working_table_NR,working_table_R, working_table)
## output according to each case, and repeat the whole process until done
## output - pDT
summary_response_pDTratio_RNAseq <- list()
summary_response_pDTratio_RNAseq[["cancer"]] <- result; rm(result)
summary_response_pDTratio_RNAseq[["drug"]] <- result; rm(result)
summary_response_pDTratio_RNAseq[["cancer_drug"]] <- result; rm(result)

summary_response_pDTratio_MC20161 <- list()
summary_response_pDTratio_MC20161[["cancer"]] <- result; rm(result)
summary_response_pDTratio_MC20161[["drug"]] <- result; rm(result)
summary_response_pDTratio_MC20161[["cancer_drug"]] <- result; rm(result)

summary_response_pDTratio_MC12399 <- list()
summary_response_pDTratio_MC12399[["cancer"]] <- result; rm(result)
summary_response_pDTratio_MC12399[["drug"]] <- result; rm(result)
summary_response_pDTratio_MC12399[["cancer_drug"]] <- result; rm(result)

## output - MKI67
summary_response_pDTratio_RNAseq_MKI67 <- list()
summary_response_pDTratio_RNAseq_MKI67[["cancer"]] <- result; rm(result)
summary_response_pDTratio_RNAseq_MKI67[["drug"]] <- result; rm(result)
summary_response_pDTratio_RNAseq_MKI67[["cancer_drug"]] <- result; rm(result)

summary_response_pDTratio_MC20161_MKI67 <- list()
summary_response_pDTratio_MC20161_MKI67[["cancer"]] <- result; rm(result)
summary_response_pDTratio_MC20161_MKI67[["drug"]] <- result; rm(result)
summary_response_pDTratio_MC20161_MKI67[["cancer_drug"]] <- result; rm(result)

summary_response_pDTratio_MC12399_MKI67 <- list()
summary_response_pDTratio_MC12399_MKI67[["cancer"]] <- result; rm(result)
summary_response_pDTratio_MC12399_MKI67[["drug"]] <- result; rm(result)
summary_response_pDTratio_MC12399_MKI67[["cancer_drug"]] <- result; rm(result)








##------- summarize NR/R ratio and vocalno plot for therapy response prediction -------
## pDT index
plot_table <- rbind(data.frame(summary_response_pDTratio_targeteddrug,DB="Lee-TD"),
                    data.frame(summary_response_pDTratio_immuno%>%filter(!(cohort%in%c("Riaz","VanAllen"))),DB="Lee-Immuno"),
                    data.frame(summary_response_pDTratio_CPB4,DB="Freeman"),
                    data.frame(summary_response_pDTratio_RNAseq$drug[,1:13],DB="CTR-RNASeq"),
                    data.frame(summary_response_pDTratio_MC20161$drug[,1:13],DB="CTR-MC20K"),
                    data.frame(summary_response_pDTratio_MC12399$drug[,1:13],DB="CTR-MC12K"))
plot_table$cohort <- gsub("_MC20161","",plot_table$cohort)
plot_table$cohort <- gsub("_MC12399","",plot_table$cohort)
plot_table$cohort <- gsub("_RNASeq","",plot_table$cohort)
# ratio_cutoff <- c(-0.1,0.1)
ratio_cutoff <- c(0,0)
plot_table <- plot_table %>% mutate(related_ratio = case_when(log2(mean_ratio) >= ratio_cutoff[2] & wilcox_p <= 0.05 ~ "R-lowerDT",
                                                              log2(mean_ratio) <= ratio_cutoff[1] & wilcox_p <= 0.05 ~ "NR-lowerDT",
                                                              TRUE ~ "Unchanged") )
plot_table_pDT <- plot_table
plot_table_pDT_auc <- plot_table_pDT %>% left_join(summary_auc_ALL[,-2], by=c("cohort", "DB")) %>% relocate(auc, .after=wilcox_p)

## MKI67
plot_table <- rbind(data.frame(summary_response_pDTratio_targeteddrug_MKI67,DB="Lee-TD"),
                    data.frame(summary_response_pDTratio_immuno_MKI67%>%filter(!(cohort%in%c("Riaz","VanAllen"))),DB="Lee-Immuno"),
                    data.frame(summary_response_pDTratio_CPB4_MKI67,DB="Freeman"),
                    data.frame(summary_response_pDTratio_RNAseq_MKI67$drug[,1:13],DB="CTR-RNASeq"),
                    data.frame(summary_response_pDTratio_MC20161_MKI67$drug[,1:13],DB="CTR-MC20K"),
                    data.frame(summary_response_pDTratio_MC12399_MKI67$drug[,1:13],DB="CTR-MC12K"))
plot_table$cohort <- gsub("_MC20161","",plot_table$cohort)
plot_table$cohort <- gsub("_MC12399","",plot_table$cohort)
plot_table$cohort <- gsub("_RNASeq","",plot_table$cohort)
# ratio_cutoff <- c(-0.1,0.1)
ratio_cutoff <- c(0,0)
plot_table <- plot_table %>% mutate(related_ratio = case_when(log2(mean_ratio) >= ratio_cutoff[2] & wilcox_p <= 0.05 ~ "R-lowerDT",
                                                              log2(mean_ratio) <= ratio_cutoff[1] & wilcox_p <= 0.05 ~ "NR-lowerDT",
                                                              TRUE ~ "Unchanged") )
plot_table_MKI67 <- plot_table


## integrate pDT index & MKI67 into one table
plot_table_pDT_MKI67 <- plot_table_pDT %>% full_join(plot_table_MKI67, by=c("cohort","DB"),suffix=c(".pDT",".MKI67"))
plot_table_pDT_MKI67_sub <- plot_table_pDT_MKI67[,c("cohort","mean_ratio.pDT","mean_ratio.MKI67","wilcox_p.pDT","wilcox_p.MKI67",
                                                    "mean_NR.pDT","mean_NR.MKI67","mean_R.pDT","mean_R.MKI67")]
plot_table_pDT_MKI67_sub$mean_ratio.pDT <- log2(plot_table_pDT_MKI67_sub$mean_ratio.pDT)
plot_table_pDT_MKI67_sub$mean_ratio.MKI67 <- log2(plot_table_pDT_MKI67_sub$mean_ratio.MKI67)


## volcano plot - pick one of the following:
plot_table <- plot_table_pDT
plot_table <- plot_table_MKI67
## then run the following:
plot_table$wilcox_q <- p.adjust(plot_table$wilcox_p,method="BH")
plot_table$related_ratio[which(plot_table$wilcox_q>0.1 & plot_table$wilcox_p<0.05)] <- "Unchanged_pvalue"

## These are the Figure 8B and S2 ##
## silence and un-silence the code for different cases
ggplot(data=plot_table, aes(x=log2(mean_ratio), y=-log10(wilcox_p), col=related_ratio, label=cohort, shape = DB)) +
  geom_vline(xintercept=ratio_cutoff, col="gray") + geom_hline(yintercept=-log10(0.05), col="gray") +
  geom_point(alpha=0.7, size = 2) + 
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5)) +
  guides(color="none")+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('log' ['2'] ('NR/R Mean pDT Ratio')), y = bquote('log' ['10'] ('Wilcoxon p-value')), col = "none", shape = "Dataset",
       title = "pDT on Therapy Response") +
  # labs(x = bquote('log' ['2'] ('NR/R Mean KI-67 Ratio')), y = bquote('log' ['10'] ('Wilcoxon p-value')), shape = "Dataset",
  #      title = "MKI67 on Therapy Response") +
  ## pDT
  geom_label_repel(data=plot_table%>%filter(related_ratio=="NR-lowerDT"),seed=22,point.padding = 0.5,nudge_x = -0.1, nudge_y = .5,
                   alpha = 0.8,max.overlaps = 15) +
  geom_label_repel(data=plot_table%>%filter(related_ratio=="R-lowerDT"),seed=20,point.padding = 0.5,nudge_x = 0.3, nudge_y = .6,
                   alpha = 0.8,max.overlaps = 15) +
  geom_label_repel(data=plot_table%>%filter(related_ratio=="Unchanged_pvalue", mean_ratio<1),seed=20,point.padding = 0.5,nudge_x = -.1, nudge_y = .5,
                   alpha = 0.8,max.overlaps = 15) +
  geom_label_repel(data=plot_table%>%filter(related_ratio=="Unchanged_pvalue", mean_ratio>1),seed=20,point.padding = 0.5,nudge_x = .6, nudge_y = .5,
                   alpha = 0.8,max.overlaps = 15) +
  ## MKI67
  # geom_label_repel(data=plot_table%>%filter(related_ratio=="NR-lowerDT",log10(wilcox_p)< -5),seed=22,point.padding = 0.5,nudge_x = -0.1, nudge_y = -.6, 
  #                  alpha = 0.8,max.overlaps = 5) +
  # geom_label_repel(data=plot_table%>%filter(related_ratio=="NR-lowerDT",log10(wilcox_p)< -3.1,log10(wilcox_p)> -5),seed=22,point.padding = 0.5,nudge_x = -0.5, nudge_y = 3.5, 
  #                  alpha = 0.8,max.overlaps = 5) +
  # geom_label_repel(data=plot_table%>%filter(related_ratio=="NR-lowerDT",log10(wilcox_p)> -3.1),seed=22,point.padding = 0.5,nudge_x = -0.1, nudge_y = 1.1, 
  #                  alpha = 0.8,max.overlaps = 5) +
  # geom_label_repel(data=plot_table%>%filter(related_ratio=="R-lowerDT"),seed=20,point.padding = 0.5,nudge_x = 0.6, nudge_y = -.6,
  #                  alpha = 0.8,max.overlaps = 15) +
  # geom_label_repel(data=plot_table%>%filter(related_ratio=="Unchanged_pvalue"),seed=20,point.padding = 0.5,nudge_x = 0.6, nudge_y = -1.2,
  #                  alpha = 0.8,max.overlaps = 15) +
# geom_label_repel(seed=20,point.padding = 0.5,nudge_x = 0.3, alpha = 0.8,
#                   nudge_y = .6, max.overlaps = 15) +

scale_color_manual(values=c(`NR-lowerDT`="blue",`Unchanged`= "grey40", `Unchanged_pvalue`="black", `R-lowerDT`="red"))


