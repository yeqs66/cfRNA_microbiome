library(lmerTest)
library(lme4)
library(sjmisc)
library(readxl)
library(tidyverse)
library(labelled)
library(tidymodels)
library(dplyr)
library(fastDummies)
library(ggResidpanel)

rm(list = ls())
#load("coldata_wgcna_flt_select.rda")


#####  #####
ori_wd<-getwd()
setwd("__PATH_001__")
load("result/Result_diffsample/0604/wgcna_file.rda")
load("__PATH_002__")

##### new wgcna  #####
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
  left_join(cfDNA_seqff,by="period_sign") %>% 
  left_join(link,by="period_sign") %>%
  #left_join(info_meta[,c(3:4,7)],by="period_sign") %>%
  dplyr::rename(ID=孕保手册号) %>%
  left_join(labelledclinical,by="ID") %>%
  dplyr::select(-sjid_kid) %>%
  #dplyr::filter(name_meta %in% metabolite_all$Name) %>%
  as.data.frame(.)
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
  dplyr::select(-ID)

coldata_wgcna[coldata_wgcna==""]<-NA

coldata_wgcna<-coldata_wgcna %>%
  column_to_rownames("period_sign")

coldata_wgcna<-coldata_wgcna[rownames(vst_cfRNA),]

coldata_wgcna[coldata_wgcna == "NA" ]<-""
tmp<-rownames(coldata_wgcna)
coldata_wgcna <- as.data.frame(sapply(coldata_wgcna, as.numeric))
rownames(coldata_wgcna)<-tmp


save(coldata_wgcna,file = "__PATH_003__")

setwd(ori_wd)
