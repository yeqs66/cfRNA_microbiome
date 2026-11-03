library(WGCNA)
library(tidyr)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(sva)
library(openxlsx)
library(readxl)
library(matrixStats)

setwd("__PATH_001__")

load("coldata_all_meta_260119.rda")
load("__PATH_002__")
load("__PATH_003__")

sample_meta<-read_excel("__PATH_004__",sheet = 10) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_meta=sampleid)
sample_cfRNA<-read_excel("__PATH_004__",sheet = 8) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_cfRNA=sampleid) %>%
  left_join(sample_meta[,c("sign","sampleid_meta")])
#sample_cfRNA<-sample_cfRNA %>%
#  dplyr::filter(sampleid_meta %in% names(rpm_meta))

##### 50% low expression filter #####
coldata_1562<-read.csv("__PATH_005__",header = T,row.names = 1) %>%
  dplyr::rename(period=trimester) %>%
  mutate(period=gsub("1","early",period),period=gsub("2","middle",period),period=gsub("3","late",period),
         ID=as.character(ID),period_sign=paste(ID,period,sep = ""))
count_mrna<-read.table("__PATH_006__",header = T,row.names = 1)
count_mirna<-read.table("__PATH_007__",header = T,row.names = 1)
count_lncrna<-read.table("__PATH_008__",header = T,row.names = 1)

count_mrna_flt<-count_mrna[,coldata_1562$sample]
count_mrna_flt<-count_mrna_flt[apply(count_mrna_flt,1,function(x) sum(x>1) > floor(ncol(count_mrna_flt)/2)),] %>%
  dplyr::rename_with(~gsub("S","",.x))
count_mirna_flt<-count_mirna[,coldata_1562$sample] 
count_mirna_flt<-count_mirna_flt[apply(count_mirna_flt,1,function(x) sum(x>1) > floor(ncol(count_mirna_flt)/10)),] %>%
  dplyr::rename_with(~gsub("S","",.x))
count_lncrna_flt<-count_lncrna[,coldata_1562$sample]
count_lncrna_flt<-count_lncrna_flt[apply(count_lncrna_flt,1,function(x) sum(x>1) > floor(ncol(count_lncrna_flt)/10)),] %>%
  dplyr::rename_with(~gsub("S","",.x))

tpm_mrna_ori<-read.table("__PATH_009__",header = T,row.names = 1) %>%
  dplyr::rename_with(~gsub("S","",.x))
tpm_mirna_ori<-read.table("__PATH_010__",header = T,row.names = 1) %>%
  dplyr::rename_with(~gsub("S","",.x))
tpm_lncrna_ori<-read.table("__PATH_011__",header = T,row.names = 1) %>%
  dplyr::rename_with(~gsub("S","",.x))

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_cfRNA %in% names(count_lncrna_flt)) %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

tpm_mrna<-tpm_mrna_ori[rownames(count_mrna_flt),sample_cfRNA$sampleid_cfRNA] %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA,sample_cfRNA$sign))
tpm_mirna<-tpm_mirna_ori[rownames(count_mirna_flt),sample_cfRNA$sampleid_cfRNA] %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA,sample_cfRNA$sign))
tpm_lncrna<-tpm_lncrna_ori[rownames(count_lncrna_flt),sample_cfRNA$sampleid_cfRNA] %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA,sample_cfRNA$sign))

##filter top varience
vars <- rowVars(tpm_mrna)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(tpm_mrna)/2)]
tpm_mrna_flt <- tpm_mrna[top_rows, ]

vars <- rowVars(tpm_mirna)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(tpm_mirna)/2)]
tpm_mirna_flt <- tpm_mirna[top_rows, ]

vars <- rowVars(tpm_lncrna)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(tpm_lncrna)/2)]
tpm_lncrna_flt <- tpm_lncrna[top_rows, ]


tpm_mrna_flt<-log2(tpm_mrna_flt+1)
tpm_mirna_flt<-log2(tpm_mirna_flt+1)
tpm_lncrna_flt<-log2(tpm_lncrna_flt+1)

##combat
#common_sample<-Reduce(intersect(rownames(vst)))
coldata_1562_flt<-coldata_1562 %>%
  dplyr::filter(period_sign %in% sample_cfRNA$sign) %>%
  dplyr::select(period_sign,disease,batch)

##combat##
mod<-model.matrix(~disease,data=coldata_1562_flt)
batch<-as.factor(coldata_1562_flt$batch)
tpm_mrna_flt_combat<-t(ComBat((tpm_mrna_flt),batch=batch,mod=mod))
tpm_mirna_flt_combat<-t(ComBat((tpm_mirna_flt),batch=batch,mod=mod))
tpm_lncrna_flt_combat<-t(ComBat((tpm_lncrna_flt),batch=batch,mod=mod))
  

##### rpm_meta ko #####
rpm_meta<-read.csv("__PATH_012__",header = T,row.names = 1)
rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign))

rpm_meta<-na.omit(rpm_meta[,sample_cfRNA$sign])

vars <- rowVars(rpm_meta)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(rpm_meta)/2)]
rpm_meta_flt <- rpm_meta[top_rows, ]
rpm_meta_flt<-as.data.frame(t(log2(rpm_meta_flt+1)))


MOFAobject<-list(KO=t(rpm_meta_flt),
                 mRNA=t(tpm_mrna_flt_combat),
                 miRNA=t(tpm_mirna_flt_combat),
                 lncRNA=t(tpm_lncrna_flt_combat))

load("coldata_all_meta_260119.rda")
coldata<-coldata_all_meta[rownames(rpm_meta_flt),]

save(MOFAobject,coldata,file = "mofa_input_ko_m_mi_lncRNA_topVarience50_260120.rda")



##### no low expression filter #####
rpm_meta<-read.csv("__PATH_012__",header = T,row.names = 1)

coldata_1562<-read.csv("__PATH_005__",header = T,row.names = 1) %>%
  dplyr::rename(period=trimester) %>%
  mutate(period=gsub("1","early",period),period=gsub("2","middle",period),period=gsub("3","late",period),
         ID=as.character(ID),period_sign=paste(ID,period,sep = ""))
count_mrna<-read.table("__PATH_006__",header = T,row.names = 1)
count_mirna<-read.table("__PATH_007__",header = T,row.names = 1)
count_lncrna<-read.table("__PATH_008__",header = T,row.names = 1)

count_mrna_flt<-count_mrna[,coldata_1562$sample]
count_mrna_flt<-count_mrna_flt %>%
  dplyr::rename_with(~gsub("S","",.x))
count_mirna_flt<-count_mirna[,coldata_1562$sample] 
count_mirna_flt<-count_mirna_flt %>%
  dplyr::rename_with(~gsub("S","",.x))
count_lncrna_flt<-count_lncrna[,coldata_1562$sample]
count_lncrna_flt<-count_lncrna_flt %>%
  dplyr::rename_with(~gsub("S","",.x))

tpm_mrna_ori<-read.table("__PATH_009__",header = T,row.names = 1) %>%
  dplyr::rename_with(~gsub("S","",.x))
tpm_mirna_ori<-read.table("__PATH_010__",header = T,row.names = 1) %>%
  dplyr::rename_with(~gsub("S","",.x))
tpm_lncrna_ori<-read.table("__PATH_011__",header = T,row.names = 1) %>%
  dplyr::rename_with(~gsub("S","",.x))

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_cfRNA %in% names(count_lncrna_flt)) %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

tpm_mrna<-tpm_mrna_ori[rownames(count_mrna_flt),sample_cfRNA$sampleid_cfRNA] %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA,sample_cfRNA$sign))
tpm_mirna<-tpm_mirna_ori[rownames(count_mirna_flt),sample_cfRNA$sampleid_cfRNA] %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA,sample_cfRNA$sign))
tpm_lncrna<-tpm_lncrna_ori[rownames(count_lncrna_flt),sample_cfRNA$sampleid_cfRNA] %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA,sample_cfRNA$sign))

##filter top varience
vars <- rowVars(tpm_mrna)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(tpm_mrna)/4)]
tpm_mrna_flt <- tpm_mrna[top_rows, ]

vars <- rowVars(tpm_mirna)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(tpm_mirna)/4)]
tpm_mirna_flt <- tpm_mirna[top_rows, ]

vars <- rowVars(tpm_lncrna)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(tpm_lncrna)/4)]
tpm_lncrna_flt <- tpm_lncrna[top_rows, ]


tpm_mrna_flt<-log2(tpm_mrna_flt+1)
tpm_mirna_flt<-log2(tpm_mirna_flt+1)
tpm_lncrna_flt<-log2(tpm_lncrna_flt+1)

##combat
#common_sample<-Reduce(intersect(rownames(vst)))
coldata_1562_flt<-coldata_1562 %>%
  dplyr::filter(period_sign %in% sample_cfRNA$sign) %>%
  dplyr::select(period_sign,disease,batch)

##combat##
mod<-model.matrix(~disease,data=coldata_1562_flt)
batch<-as.factor(coldata_1562_flt$batch)
tpm_mrna_flt_combat<-t(ComBat((tpm_mrna_flt),batch=batch,mod=mod))
tpm_mirna_flt_combat<-t(ComBat((tpm_mirna_flt),batch=batch,mod=mod))
tpm_lncrna_flt_combat<-t(ComBat((tpm_lncrna_flt),batch=batch,mod=mod))


## rpm_meta ko
rpm_meta<-read.csv("__PATH_012__",header = T,row.names = 1)
rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign))

rpm_meta<-na.omit(rpm_meta[,sample_cfRNA$sign])

vars <- rowVars(rpm_meta)                # 比 apply 快一个量级
top_rows <- order(vars, decreasing = TRUE)[1:floor(nrow(rpm_meta)/4)]
rpm_meta_flt <- rpm_meta[top_rows, ]
rpm_meta_flt<-as.data.frame(t(log2(rpm_meta_flt+1)))


MOFAobject<-list(KO=t(rpm_meta_flt),
                 mRNA=t(tpm_mrna_flt_combat),
                 miRNA=t(tpm_mirna_flt_combat),
                 lncRNA=t(tpm_lncrna_flt_combat))

load("coldata_all_meta_260119.rda")
coldata<-coldata_all_meta[rownames(rpm_meta_flt),]

save(MOFAobject,coldata,file = "mofa_input_ko_m_mi_lncRNA_topVarience25_noLowExpressionFilter_260129.rda")
