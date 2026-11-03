library(tidyverse)
library(readxl)
library(clusterProfiler)

#load("__PATH_001__")
#load("__PATH_002__")
load("__PATH_003__")
gmt_microbiome<-read.table("__PATH_004__",header = T,sep = "\t")

deg_micro_ko<-read.csv("__PATH_005__") 
deg_micro_ko_flt<-deg_micro_ko %>%
  dplyr::filter(p.adj<0.05)

enrich_kegg<-data.frame()
gsea_ko_kegg<-data.frame()
for (i in unique(deg_micro_ko_flt$contrast)) {
  tmp<-deg_micro_ko_flt %>%
    dplyr::filter(contrast==i)
  KO<-tmp$feature
  
  enrich_tmp<-enricher(KO,
                       pvalueCutoff = 1, 
                       pAdjustMethod = "BH", 
                       TERM2GENE = gmt_microbiome, 
                       universe = unique(deg_micro_ko$feature))
  
  enrich_tmp<-enrich_tmp@result %>%
    mutate(contrast=i)
  enrich_kegg<-rbind(enrich_kegg,enrich_tmp)
  
  tmp<-deg_micro_ko %>%
    dplyr::filter(contrast==i) %>%
    ##log2MedianRatio
    dplyr::filter(!(is.na(log2fc) | log2fc=="Inf" | log2fc=="-Inf")) %>%
    arrange(desc(log2fc))
    ##log10padj
    #mutate(log10padj=-log10(p.adj)) %>%
    #dplyr::filter(!is.na(log10padj)) %>%
    #arrange(desc(log10padj))
    
  KO<-tmp$log2fc
  names(KO)<-tmp$feature
  
  #set.seed(3490)
  gsea_ko_tmp<-GSEA(KO,pvalueCutoff = 1,pAdjustMethod = "BH",TERM2GENE = gmt_microbiome,seed = F)
  gsea_ko_tmp<-gsea_ko_tmp@result %>%
    mutate(contrast=i)
  gsea_ko_kegg<-rbind(gsea_ko_kegg,gsea_ko_tmp)
}







