library(tidyverse)
library(openxlsx)
library(readxl)

##arrange DEG of mRNA
tmp<-read.csv("limma_349_batch_only.csv") %>%
  dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("1","1st",trimester),
         trimester=gsub("2","2nd",trimester),
         trimester=gsub("3","3rd",trimester),
         trimester=gsub("all","overall",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 1_Differential_analysis_results_for_cf-mRNA.xlsx")

###arrange DEG of miRNA
tmp<-data.frame()
for (i in 1:4) {
  temp<-read_excel("limma_349_miRNA_batch_only.xlsx",sheet = i)
  tmp<-rbind(tmp,temp)
}

tmp<-tmp %>%
  dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("1","1st",trimester),
         trimester=gsub("2","2nd",trimester),
         trimester=gsub("3","3rd",trimester),
         trimester=gsub("all","overall",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 2_Differential_analysis_results_for_cf-miRNA.xlsx")

##arrange DEG of lncRNA
tmp<-data.frame()
for (i in 1:4) {
  temp<-read_excel("limma_349_lncRNA_batch_only.xlsx",sheet = i)
  tmp<-rbind(tmp,temp)
}

tmp<-tmp %>%
  dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("1","1st",trimester),
         trimester=gsub("2","2nd",trimester),
         trimester=gsub("3","3rd",trimester),
         trimester=gsub("all","overall",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 3_Differential_analysis_results_for_cf-lncRNA.xlsx")

##arrange DEG of microbiome KO
file<-dir("microbiome ko gsea复现/new_deg/")
tmp<-data.frame()
for (i in file) {
  trimester<-gsub(".*zinb_","",i)
  trimester<-gsub("_rpk.*","",trimester)
  temp<-read.csv(paste0("microbiome ko gsea复现/new_deg/",i)) %>%
    mutate(trimester=trimester)
  tmp<-rbind(tmp,temp)
}

tmp<-tmp %>%
  #dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("early","1st",trimester),
         trimester=gsub("middle","2nd",trimester),
         trimester=gsub("late","3rd",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 11_Differential_analysis_results_for_microbiome_derived_KO.xlsx")

##arrange GSEA of cfRNA
tmp<-read_excel("gsea_deg_vs_kegg_251212.xlsx") %>%
  dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("1","1st",trimester),
         trimester=gsub("2","2nd",trimester),
         trimester=gsub("3","3rd",trimester),
         trimester=gsub("all","overall",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 4_Gene_Set_Enrichment_Analysis_of_KEGG_for_cfRNA.xlsx")

##arrange GO of cfRNA
tmp<-read_excel("gsea_deg_go_251230.xlsx") %>%
  dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("1","1st",trimester),
         trimester=gsub("2","2nd",trimester),
         trimester=gsub("3","3rd",trimester),
         trimester=gsub("all","overall",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 5_Gene_Set_Enrichment_Analysis_of_GO_for_cfRNA.xlsx")

##arrange GSEA of KO
load("microbiome ko gsea复现/gsea_deg_KO_260212.rda")
tmp<-gsea %>%
  dplyr::rename(trimester=period) %>%
  mutate(trimester=gsub("early","1st",trimester),
         trimester=gsub("middle","2nd",trimester),
         trimester=gsub("late","3rd",trimester))
write.xlsx(tmp,"manu_v1.0/Supplementary Table 12_Gene_Set_Enrichment_Analysis_of_microbiome_derived_KO.xlsx")


##### WGCNA module #####
module_gene<-read.csv("投稿期刊/cell reports/material/Supplementary Table 6_WGCNA_module_gene_set.csv",row.names=1) %>%
  mutate(moduleLabels=moduleLabels+1,moduleLabels=as.character(moduleLabels),
         moduleColors=gsub("^","ME",moduleColors)) %>%
  dplyr::select(moduleLabels)
module_cor_r<-read.csv("投稿期刊/cell reports/material/Supplementary Table 7_WGCNA_moduleTraitCorrelation.csv",row.names=1) %>%
  dplyr::rename(!!!setNames(module_gene$moduleColors, module_gene$moduleLabels))
module_cor_p<-read.csv("投稿期刊/cell reports/material/Supplementary Table 8_WGCNA_moduleTraitCorPvalue.csv",row.names=1) %>%
  dplyr::rename(!!!setNames(module_gene$moduleColors, module_gene$moduleLabels))
module_cor_fdr<-read.csv("投稿期刊/cell reports/material/Supplementary Table 9_WGCNA_moduleTraitCorFDR.csv",row.names=1) %>%
  dplyr::rename(!!!setNames(module_gene$moduleColors, module_gene$moduleLabels))

write.csv(module_gene,"投稿期刊/cell reports/material/simple_scirep/Supplementary Table 6.csv",row.names = T)
write.csv(module_cor_r,file="投稿期刊/cell reports/material/simple_scirep/Supplementary Table 7.csv",row.names = T)
write.csv(module_cor_p,file="投稿期刊/cell reports/material/simple_scirep/Supplementary Table 8.csv",row.names = T)
write.csv(module_cor_fdr,file="投稿期刊/cell reports/material/simple_scirep/Supplementary Table 9.csv",row.names = T)

