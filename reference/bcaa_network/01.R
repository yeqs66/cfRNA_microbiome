library(latex2exp)
library(grid)
library(openxlsx)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(ggplot2)
library(gridExtra)
library(fgsea)
library(ggord)
library(FactoMineR)

##### load cfRNA #####
load("__PATH_001__")
#pathways <- pathways$transciptome_ENTREZID[grep("(KEGG)",names(pathways$transciptome_ENTREZID))]
pathways <- do.call(rbind, lapply(seq_along(pathways$transciptome_ENTREZID), function(i) {
  data.frame(name=rep(names(pathways$transciptome_ENTREZID[i]),length(pathways$transciptome_ENTREZID[[i]])), value=pathways$transciptome_ENTREZID[[i]], stringsAsFactors=FALSE)
})) %>%
  dplyr::filter(grepl("(KEGG)",name))
symbol_ensg<-read.table("__PATH_002__")

select_pathway_rna<-symbol_ensg %>%
  dplyr::filter(V1 %in% pathways$value[which(grepl("Valine, leucine and isoleucine degradation",pathways$name))]) %>%
  dplyr::pull(V2)

coldata_all<-read.csv("__PATH_003__",row.names = 1) %>%
  mutate(visit=case_when(trimester=="1"~"early",
                         trimester=="2"~"middle",
                         trimester=="3"~"late"),
         ID=as.character(ID),
         sign=paste(ID,visit,sep="_"))

tpm_1562<-read.csv("__PATH_004__",row.names = 1,header = T) %>%
  dplyr::select(coldata_all$sample) %>%
  dplyr::rename(!!!setNames(coldata_all$sample, coldata_all$sign)) %>%
  dplyr::filter(rownames(.) %in% select_pathway_rna) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("sign") %>%
  dplyr::mutate(across(where(is.numeric), ~log2(.x + 1)))

## pca
pca_tpm<-tpm_1562 %>%
  column_to_rownames("sign") %>%
  PCA(graph = FALSE)

pca_tpm<-data.frame(PC1=pca_tpm$ind$coord[, 1],PC2=pca_tpm$ind$coord[, 2]) %>%
  dplyr::rename_with(~gsub("^","cfRNA_",.x)) %>%
  rownames_to_column("sign")




##### load metabolome #####
load("__PATH_005__")
load("__PATH_006__")
load("__PATH_007__")

coldata_349 <- read.csv("__PATH_003__",row.names = 1) %>%
  dplyr::mutate(sample=gsub("S","",sample),
                ID=as.character(ID),
                trimester=case_when(trimester==1 ~ "early",
                                    trimester==2 ~ "middle",
                                    trimester==3 ~ "late"),
                sign=paste(ID,trimester,sep = "_")
                ) %>%
  dplyr::rename(sample_cfRNA=sample) %>%
  dplyr::select(-batch) ## 不需要cfRNA的batch变量

# 匹配 349 例 cfRNA 样本，并一次性带入 5 个协变量
metabolics_matrix <- metabolics_matrix %>%
  dplyr::mutate(sign=paste(ID,trimester,sep = "_"),
                ID=as.character(ID)) %>%
  left_join(coldata_349[,c("sign","sample_cfRNA")],by="sign") %>%
  dplyr::filter(!is.na(sample_cfRNA))

select_meta <- read.csv("__PATH_008__")
select_meta_SFCA <- read.csv("__PATH_009__")
bcaa_id <- unique(c(mergeinfo$Vname[which(mergeinfo$Compound %in% select_meta$BCAA)],
    "V34",   # L-Valine
    "V54",   # Isoleucine
    "V55",   # L-Leucine
    "V415",  # alpha-Ketoisovaleric acid
    "V181",  # Ketoleucine
    "V182",  # 3-Methyl-2-oxovaleric acid，构型说明见下
    "V39",   # 3-Hydroxyisovaleric acid
    "V170",  # Acetoacetic acid
    "V337"))
scfa_id <- mergeinfo$Vname[which(mergeinfo$Vname %in% select_meta_SFCA$Vname)]

trans_meta_name <- mergeinfo %>%
  dplyr::filter(Vname %in% c(bcaa_id)) %>%
  dplyr::select(Vname,Compound,类别) %>%
  dplyr::rename(category=类别) %>%
  mutate(category=case_when(category=="Fatty acids" ~ "SCFA",
                            category=="Amino acids, peptides, and analogues" ~ "BCAA"))

metabolics_matrix_val <- metabolics_matrix %>%
  dplyr::select(sample_cfRNA,ID,trimester,group,bcaa_id) %>%
  dplyr::rename(!!!setNames(trans_meta_name$Vname,trans_meta_name$Compound)) %>%
  dplyr::mutate(trimester=factor(trimester,levels=c("early","middle","late")),
                sign=paste(ID,trimester,sep = "_")) %>%
  dplyr::select(-c(sample_cfRNA,ID,trimester,group))




##### load microbiome #####
ori_wd<-getwd()
setwd("__PATH_010__")

load("kogmt.rda")
load("klist_to_module_20240325_1524_updated.rda")
gmt <- merge(kogmt,klist_to_module,by=1) %>%
  dplyr::select(mod,V2)

# gmt for pathway
pathgmt <- read.table("ko00001.forR.KO_only.gmt", header=FALSE)
setwd(ori_wd)

load("__PATH_011__")
load("microbiome ko gsea复现/microbiome_sampleinfo_GDM1999s_20240318.rda")

select_pathway_micro<-pathgmt %>%
  dplyr::filter(V1=="00280_Valine,_leucine_and_isoleucine_degradation_[PATH:ko00280]") %>%
  dplyr::pull(V2)

matrix_meta <- read.csv("microbiome ko gsea复现/TJGDM1962s_ko_tri_95_rpk.csv", row.names = 1,header=T) %>% 
  dplyr::filter(rownames(.) %in% select_pathway_micro)

microbiome_sampleinfo<-microbiome_sampleinfo %>%
  mutate(sign=paste(ID,trimester,sep="_")) %>%
  dplyr::filter(sample %in% colnames(matrix_meta))

matrix_meta<-matrix_meta %>%
  dplyr::rename(!!!setNames(microbiome_sampleinfo$sample, microbiome_sampleinfo$sign)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("sign")
matrix_meta_log<-matrix_meta %>%
  dplyr::mutate(across(where(is.numeric), ~log2(.x + 1))) %>%
  left_join(tpm_1562,by="sign") %>%
  left_join(metabolics_matrix_val,by="sign")

## pca
pca_micro<-matrix_meta %>%
  dplyr::mutate(across(where(is.numeric), ~log2(.x + 1))) %>%
  column_to_rownames("sign")
dat.pca_micro <- PCA(pca_micro, graph = FALSE)

pca_micro<-data.frame(PC1=dat.pca_micro$ind$coord[, 1],PC2=dat.pca_micro$ind$coord[, 2]) %>%
  dplyr::rename_with(~gsub("^","micro_",.x)) %>%
  rownames_to_column("sign")



##### PCA correlation #####
matrix_pca <- metabolics_matrix %>%
  dplyr::select(sample_cfRNA,ID,trimester,group,bcaa_id) %>%
  dplyr::rename(!!!setNames(trans_meta_name$Vname,trans_meta_name$Compound)) %>%
  dplyr::mutate(trimester=factor(trimester,levels=c("early","middle","late")),
                sign=paste(ID,trimester,sep = "_")) %>%
  dplyr::select(-c(sample_cfRNA,ID,trimester,group)) %>%
  left_join(pca_tpm,by="sign") %>%
  left_join(pca_micro,by="sign") 


