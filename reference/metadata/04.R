library(tidyverse)
library(readxl)

rpm_meta<-read.csv("__PATH_001__",header = T,row.names = 1)
sample_meta<-read_excel("__PATH_002__",sheet = 10) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_meta=sampleid)
sample_cfRNA<-read_excel("__PATH_002__",sheet = 8) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_cfRNA=sampleid) %>%
  left_join(sample_meta[,c("sign","sampleid_meta")])

load("__PATH_003__")
coldata_1562<-read.csv("__PATH_004__",row.names = 1) %>% 
  mutate(ID=as.character(ID),sample=gsub("^S","",sample))

coldata_1476<-coldata_1562 %>%
  dplyr::filter(sample %in% intersect(sample_cfRNA$sampleid_cfRNA,sample)) %>%
  mutate(sample=gsub("^","S",sample),
         trimester=gsub("1","early",trimester),
         trimester=gsub("2","middle",trimester),
         trimester=gsub("3","late",trimester),
         sign=paste0(ID,trimester,sep=""))

limma_voom_multiRNA_1476<-limma_voom_multiRNA %>%
  dplyr::select(coldata_1476$sample) %>%
  dplyr::rename(!!!setNames(coldata_1476$sample,coldata_1476$sign))
limma_voom_multiRNA_1476<-limma_voom_multiRNA_1476[,coldata_1476$sign]


rpm_meta<-read.csv("__PATH_001__",header = T,row.names = 1)
rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  rownames_to_column("Kname") %>%
  mutate(Kname=gsub("^","Kname_",Kname)) %>%
  column_to_rownames("Kname")
rpm_meta<-rpm_meta[,coldata_1476$sign]
rpm_meta<-as.data.frame(log2(rpm_meta+1))
all(colnames(rpm_meta)==coldata_1476$sign)

deg_ko_all<-read.csv("microbiome ko gsea复现/new_deg/TJGDM1962s_ko_tri_95.zinb_overall_rpk.csv") %>%
  dplyr::mutate(period="overall")
deg_ko_early<-read.csv("microbiome ko gsea复现/new_deg/TJGDM1962s_ko_tri_95.zinb_early_rpk.csv") %>%
  dplyr::mutate(period="early")
deg_ko_middle<-read.csv("microbiome ko gsea复现/new_deg/TJGDM1962s_ko_tri_95.zinb_middle_rpk.csv") %>%
  dplyr::mutate(period="middle")
deg_ko_late<-read.csv("microbiome ko gsea复现/new_deg/TJGDM1962s_ko_tri_95.zinb_late_rpk.csv") %>%
  dplyr::mutate(period="late")
deg_ko_merge<-rbind(deg_ko_early,deg_ko_middle,deg_ko_late,deg_ko_all)
ko_sig<-deg_ko_merge$feature[which(deg_ko_merge$p.adj<0.05)]
deg_ko_merge<-deg_ko_merge  %>%
  dplyr::filter(feature %in% ko_sig) %>%
  mutate(feature=gsub("^","Kname_",feature)) %>%
  pivot_wider(id_cols=-c(contrast,difference,p.value,log2fc),names_from = period,values_from = c(log2MedianRatio,p.adj),values_fill = 0)


rpm_meta_sig<-rpm_meta[unique(deg_ko_merge$feature),]

save()



