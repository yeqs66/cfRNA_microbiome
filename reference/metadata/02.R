library(WGCNA)
library(tidyr)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(sva)
library(openxlsx)
library(fastDummies)
library(preprocessCore)
library(labelled)

rm(list = ls())

##### vst_cfRNA #####
###matrix-cfRNA
#load("__PATH_001__")
coldata_1562<-read.csv("__PATH_002__",header = T,row.names = 1) %>%
  dplyr::rename(period=trimester) %>%
  mutate(period=gsub("1","early",period),period=gsub("2","middle",period),period=gsub("3","late",period),
         ID=as.character(ID),period_sign=paste(ID,period,sep = ""))

matrix_rna<-read.table("__PATH_003__",row.names = 1,header = T) %>%
  #dplyr::rename_with(~ gsub("^.....","",.x)) %>%
  dplyr::select(coldata_1562$sample) %>%
  dplyr::rename(!!!setNames(coldata_1562$sample,coldata_1562$period_sign))

matrix_mirna<-read.table("__PATH_004__",row.names = 1,header = T) %>% 
  #rename_all(~ gsub("X", "S", .)) %>%
  dplyr::select(coldata_1562$sample) %>%
  dplyr::rename(!!!setNames(coldata_1562$sample,coldata_1562$period_sign))

matrix_lncrna<-read.table("__PATH_005__",row.names = 1,header = T) %>% 
  #rename_all(~ gsub("X", "S", .)) %>%
  dplyr::select(coldata_1562$sample) %>%
  dplyr::rename(!!!setNames(coldata_1562$sample,coldata_1562$period_sign))


##cfRNA
matrix_rna_flt<-matrix_rna%>%dplyr::select(coldata_1562$period_sign)
matrix_rna_flt<-matrix_rna_flt[apply(matrix_rna_flt,1,function(x) sum(x>1) > floor(ncol(matrix_rna_flt)/2)),]
matrix_rna_flt<-round(matrix_rna_flt)
vst_rna<-varianceStabilizingTransformation(as.matrix(matrix_rna_flt))
vst_rna<-t(vst_rna)

matrix_mirna_flt<-matrix_mirna%>%dplyr::select(coldata_1562$period_sign)
matrix_mirna_flt<-matrix_mirna_flt[apply(matrix_mirna_flt,1,function(x) sum(x>1) > floor(ncol(matrix_mirna_flt)/2)),]
matrix_mirna_flt<-round(matrix_mirna_flt)
vst_mirna<-varianceStabilizingTransformation(as.matrix(matrix_mirna_flt))
vst_mirna<-t(vst_mirna)

matrix_lncrna_flt<-matrix_lncrna%>%dplyr::select(coldata_1562$period_sign)
matrix_lncrna_flt<-matrix_lncrna_flt[apply(matrix_lncrna_flt,1,function(x) sum(x>1) > floor(ncol(matrix_lncrna_flt)/2)),]
matrix_lncrna_flt<-round(matrix_lncrna_flt)
vst_lncrna<-varianceStabilizingTransformation(as.matrix(matrix_lncrna_flt))
vst_lncrna<-t(vst_lncrna)

##聚类查看离群值
#mRNA
sampleTree = hclust(dist(vst_rna), method = "average")
sizeGrWindow(12,9) #视图sampleTree = hclust(dist(vst_count_mRNA_349), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 240, col = "red") #划定需要剪切的枝长
clust = cutreeStatic(sampleTree, cutHeight = 240, minSize = 10)
table(clust)   
keepSamples <- (clust==1)  
vst_rna_flt <- vst_rna[keepSamples, ] 

#miRNA
sampleTree = hclust(dist(vst_mirna), method = "average")
sizeGrWindow(12,9) #视图sampleTree = hclust(dist(vst_count_mRNA_349), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 70, col = "red") #划定需要剪切的枝长
clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
table(clust)   
keepSamples_mi <- (clust==1)  
vst_mirna_flt <- vst_mirna[keepSamples_mi, ] 

#lncrna
sampleTree = hclust(dist(vst_lncrna), method = "average")
sizeGrWindow(12,9) #视图sampleTree = hclust(dist(vst_count_mRNA_349), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 150, col = "red") #划定需要剪切的枝长
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)   
keepSamples_lnc <- (clust==1)  
vst_lncrna_flt <- vst_lncrna[keepSamples_lnc, ] 

##merge
common_sample<-Reduce(intersect,list(rownames(vst_rna_flt),rownames(vst_mirna_flt),rownames(vst_lncrna_flt)))
vst_rna_flt<-vst_rna_flt[which(rownames(vst_rna_flt) %in% common_sample),]
vst_mirna_flt<-vst_mirna_flt[which(rownames(vst_mirna_flt) %in% common_sample),]
vst_lncrna_flt<-vst_lncrna_flt[which(rownames(vst_lncrna_flt) %in% common_sample),]

all(rownames(vst_mirna_flt)==rownames(vst_rna_flt))
all(rownames(vst_mirna_flt)==rownames(vst_lncrna_flt))

##combat
#common_sample<-Reduce(intersect(rownames(vst)))
coldata_1562_flt<-coldata_1562 %>%
  dplyr::filter(period_sign %in% common_sample) %>%
  dplyr::select(period_sign,disease,batch)

##combat##
mod<-model.matrix(~disease,data=coldata_1562_flt)
batch<-as.factor(coldata_1562_flt$batch)
if (any(coldata_1562_flt$period_sign!=names(vst_rna_flt))) {
  stop("sample is inconsistant")
} else {
  vst_rna_flt_combat<-t(ComBat(t(vst_rna_flt),batch=batch,mod=mod))
  vst_mirna_flt_combat<-t(ComBat(t(vst_mirna_flt),batch=batch,mod=mod))
  vst_lncrna_flt_combat<-t(ComBat(t(vst_lncrna_flt),batch=batch,mod=mod))
  
  vst_cfRNA<-cbind(vst_rna_flt_combat,vst_mirna_flt_combat,vst_lncrna_flt_combat)
}

save(vst_cfRNA,file = "__PATH_006__")


##### coldata_wgcna #####
ori_wd<-getwd()
setwd("__PATH_007__")
load("__PATH_006__")
#load("__PATH_008__")

coldata_rna_1562<-read.table("../../../cfRNA_coldata_1562_batch.txt",header = T) %>%
  dplyr::rename(period=trimester) %>%
  mutate(ID=as.character(ID),period_sign=paste(ID,period,sep = ""))
cfDNA_seqff<-read_excel("../../../../GDM VS control/临床信息/cfDNA采样孕周.xlsx") %>%
  mutate(period_sign=gsub("F","early",period_sign),period_sign=gsub("S","middle",period_sign),period_sign=gsub("T","late",period_sign))
sample_period<-read_excel("../../../../GDM VS control/临床信息/样本编号-孕保手册号-采样时间.xlsx",col_names = F) %>%
  dplyr::rename("sample" = colnames(.)[1],"period_sign" = colnames(.)[2]) %>%
  mutate(period_sign=gsub("F","early",period_sign),period_sign=gsub("S","middle",period_sign),period_sign=gsub("T","late",period_sign))
link<-read.csv("../../../../GDM VS control/临床信息/样本编号-孕保手册号.csv",fileEncoding = "GBK",row.names = 1) %>% 
  dplyr::rename("孕保手册号" = colnames(.)[1], "sample" = colnames(.)[2]) %>%
  mutate(孕保手册号=as.character(孕保手册号)) %>%
  left_join(sample_period,by="sample") %>% 
  left_join(cfDNA_seqff,by="period_sign") %>%
  dplyr::rename(name_RNA=sample,batch_rna=batch)

#cfRNA leading
load("../../../../GDM VS control/临床信息/labelledclinical_20240223.rda") 
labelledclinical<-labelledclinical %>%
  mutate(ID=as.character(ID)) %>%
  mutate(blood_rountine_gws_1=as.numeric(gsub("\\+.*","",blood_rountine_gws)),blood_rountine_gws_2=as.numeric(gsub(".*\\+","",blood_rountine_gws)),tmp=as.numeric(7*blood_rountine_gws_1+blood_rountine_gws_2),blood_rountine_gws=tmp) %>%
  dplyr::select(-blood_rountine_gws_1,-blood_rountine_gws_2,-tmp) %>%
  dplyr::rename(blood_rountine_TransDays_gws=blood_rountine_gws) %>%
  dplyr::rename(first_sampletime=孕早期, middle_sampletime=孕中期, late_sampletime=孕晚期)
coldata_meta<-coldata_rna_1562[,c(1,2,3,6)] %>% 
  #left_join(cfDNA_seqff,by="period_sign") %>% 
  left_join(link,by="period_sign") %>%
  #left_join(info_meta[,c(3:4,7)],by="period_sign") %>%
  #dplyr::rename(ID=孕保手册号) %>%
  left_join(labelledclinical,by="ID") %>%
  dplyr::select(-sjid_kid) %>%
  #dplyr::filter(name_meta %in% metabolite_all$Name) %>%
  as.data.frame()
rownames(coldata_meta)<-coldata_meta$period_sign
coldata_meta<-remove_labels(coldata_meta)

coldata_meta<-coldata_meta %>% 
  #mutate_all(~str_replace_all(., "#NULL!", "")) %>%
  #mutate(late_sampletime=replace(late_sampletime,is.na(late_sampletime),"")) %>%
  #mutate_all(~ifelse(. == "", NA, .)) %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  mutate(birthday_mo=as.numeric(as.Date(birthday_mo)),
         birthday_fa=as.numeric(as.Date(birthday_fa)),
         birthday_kid=as.numeric(as.Date(birthday_kid)),
         last_menstrual_period.x=as.numeric(as.Date(last_menstrual_period.x)),
         delivery_date=as.numeric(as.Date(delivery_date)),
         last_menstrual_period.y=as.numeric(as.Date(last_menstrual_period.y)),
         first_sampletime=as.numeric(as.Date(first_sampletime)),
         middle_sampletime=as.numeric(as.Date(middle_sampletime)),
         late_sampletime=as.numeric(as.Date(late_sampletime))) %>%
  #dplyr::select(-c(ter_pre_gws,ter_pre_reasons,ter_pre_method,终止妊娠日期,stop_preg)) %>% 
  mutate(birth_defect=as.numeric(birth_defect)) %>%
  dplyr::select(-sample,-ID)


var_type<-read_excel("wgcna-性状变量类型.xlsx") %>%
  na.omit()

#

numeric_col<-c("period","family_income","sleep_quality_1year_beforepreg","sleep_quality_afterpreg","urine_occult_blood",
               "urine_protein","urine_glucose","urine_acetone_bodies","urine_white_cell","FT4_group",
               "tea_freq_prepreg","tea_freq_afterpreg","coffee_freq_prepreg","coffee_freq_afterpreg",
               "sodas_freq_prepreg","sodas_freq_afterpreg","autism_12month","autism_18month")

#rename_col<-data.frame(ori=names(coldata_meta),rep=gsub(" ","_",names(coldata_meta))) %>%
#  mutate(rep=gsub("\\+","plus",rep))

coldata_wgcna<-coldata_meta %>%
  mutate(period_sign=rownames(.),ID=gsub("[a-zA-Z]","",period_sign),
         gestational_wks_delivery=as.numeric(gestational_wks_delivery),
         DM_father_Familyhistory_mo=gsub("1","No",DM_father_Familyhistory_mo),
         DM_father_Familyhistory_mo=gsub("2","Yes",DM_father_Familyhistory_mo),
         DM_mother_Familyhistory_mo=gsub("1","No",DM_mother_Familyhistory_mo),
         DM_mother_Familyhistory_mo=gsub("2","Yes",DM_mother_Familyhistory_mo),
         DM_siblings_Familyhistory_mo=gsub("1","No",DM_siblings_Familyhistory_mo),
         DM_siblings_Familyhistory_mo=gsub("2","Yes",DM_siblings_Familyhistory_mo),
         birth_defect=gsub("1","No",birth_defect),
         birth_defect=gsub("2","Yes",birth_defect)) %>%
  mutate(across(everything(), ~ ifelse(. == "GDM", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Control", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "early", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "middle", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "late", 3, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "male", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "female", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "No", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Yes", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "bad", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "medium", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "good", 3, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "done 1st questionnaire", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "China", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Canada", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Han", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "minority", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Foreigner", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Above college or others", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Below senior high school", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Vaginal delivery", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Caesarean section", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "other_not_specified", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "A", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "AB", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "B", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "O", 3, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Positive", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Negative", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "-", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "±", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "+", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "++", 3, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "+++", 4, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Low", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Normal", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "High", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "lower than 50,000 CNY", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "50,000 to 99,999 CNY", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "100,000 to 199,999 CNY", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "200,000 TO 299,999 CNY", 3, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "300,000 to 399,000 CNY", 4, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "400,000 to 499,999 CNY", 5, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "greater than 500,000 CNY", 6, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "never", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "sometimes", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "once a week", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "2-3 times a week", 3, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "above 3 times a week", 4, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Term delivery", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "Preterm delivery", 1, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "exclusive formula fed", 2, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "exclusive breastfeeding", 0, .))) %>%
  mutate(across(everything(), ~ ifelse(. == "mixed feeding", 1, .))) %>%
  #mutate(blood_rountine_gws=gsub("\\+.*","",blood_rountine_gws)) %>%  ##已经更换为天数
  #mutate(batch_rna=gsub("149",0,batch_rna),batch_rna=gsub("200",1,batch_rna)) %>%
  mutate_at(vars(numeric_col),as.numeric) %>%
  #mutate(across(everything(), ~str_replace_all(., " ", "_"))) %>%
  dummy_cols(select_columns = c("ABO_type","Resident_address","feeding_type_28days"),ignore_na = T) %>%
  dplyr::select(-feeding_type_28days,-Resident_address,-ABO_type) %>%
  dplyr::select(-c(ID,孕保手册号,name_RNA))

coldata_wgcna[coldata_wgcna==""]<-NA

coldata_wgcna<-coldata_wgcna %>%
  column_to_rownames("period_sign")

coldata_wgcna<-coldata_wgcna[rownames(vst_cfRNA),]

coldata_wgcna[coldata_wgcna == "NA" ]<-""
tmp<-rownames(coldata_wgcna)
coldata_wgcna <- as.data.frame(sapply(coldata_wgcna, as.numeric))
rownames(coldata_wgcna)<-tmp


save(coldata_wgcna,file = "__PATH_001__")



##### coldata_wgcna add meta ko #####
load("__PATH_001__")
load("__PATH_009__")

sample_meta<-read_excel("__PATH_010__",sheet = 10) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_meta=sampleid)
sample_cfRNA<-read_excel("__PATH_010__",sheet = 8) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_cfRNA=sampleid) %>%
  left_join(sample_meta[,c("sign","sampleid_meta")])


rpm_meta<-read.csv("__PATH_011__",header = T,row.names = 1)

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  t() %>% as.data.frame()

rpm_meta<-na.omit(rpm_meta[rownames(coldata_wgcna),])
missing_sample <- setdiff(rownames(coldata_wgcna), rownames(rpm_meta))
rpm_meta[missing_sample, ] <- 0
rpm_meta<-rpm_meta[rownames(coldata_wgcna),]

##log2 transform
#rpm_meta<-as.data.frame(log2(rpm_meta+1))

##分位数标准化
#rpm_meta_norm<-as.data.frame(log2(rpm_meta+1))
#rpm_meta_norm <- normalize.quantiles(as.matrix(rpm_meta))
#rownames(rpm_meta_norm) <- rownames(rpm_meta)
#colnames(rpm_meta_norm) <- colnames(rpm_meta)
#rpm_meta<-rpm_meta_norm

all(rownames(rpm_meta)==rownames(coldata_wgcna))

coldata_wgcna_ko<-cbind(coldata_wgcna,rpm_meta)



##### coldata_wgcna add meta strain #####
rpm_meta<-read.csv("__PATH_012__",header = T,row.names = 1)

rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  t() %>% as.data.frame()

rpm_meta<-na.omit(rpm_meta[rownames(coldata_wgcna),])
missing_sample <- setdiff(rownames(coldata_wgcna), rownames(rpm_meta))
rpm_meta[missing_sample, ] <- 0
rpm_meta<-rpm_meta[rownames(coldata_wgcna),]

##log2 transform
#rpm_meta<-as.data.frame(log2(rpm_meta+1))

##分位数标准化
#rpm_meta_norm<-as.data.frame(log2(rpm_meta+1))
#rpm_meta_norm <- normalize.quantiles(as.matrix(rpm_meta))
#rownames(rpm_meta_norm) <- rownames(rpm_meta)
#colnames(rpm_meta_norm) <- colnames(rpm_meta)
#rpm_meta<-rpm_meta_norm

all(rownames(rpm_meta)==rownames(coldata_wgcna))

coldata_wgcna_strain<-cbind(coldata_wgcna,rpm_meta)



save(vst_cfRNA,coldata_wgcna,coldata_wgcna_ko,coldata_wgcna_strain,file = "__PATH_013__")

save(coldata_wgcna_ko,coldata_wgcna_strain,file="__PATH_014__")

###### sft #####
##sft
powers = c(c(1:10), seq(from = 10, to=19, by=2))
sft = pickSoftThreshold(vst_cfRNA, powerVector = powers, verbose = 5)

cairo_pdf(filename = "__PATH_015__",
          width = 10,height = 6,family = "Arial")
par(mfrow = c(1,2),
    mar = c(4.5, 4.5, 3, 1.5),  # 缩小边距，特别是右边距
    pin = c(3, 3))  
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8 ,col="red") ##截止值0.8

##平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$powerEstimate

dev.off()
