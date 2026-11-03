library(tidyverse)

load("test_prediction_260202.rda")
gene_in_module<-read.csv("wgcna/cfRNA_microbiome/m_mi_lnc/ko/wgcga_signed_signed_hybrid_pearson/module_gene-signed-signed_hybrid-pearson.csv")

limma_voom_module13<-limma_voom_multiRNA_1476 %>%
  dplyr::filter(rownames(.) %in% gene_in_module$X[which(gene_in_module$moduleColors=="tan")])



load("__PATH_001__")
load("__PATH_002__")

gender<-labelledclinical[,c("sjid_kid","gender_kid")] %>%
  mutate(gender_kid=gsub("^male$","M",gender_kid),
         gender_kid=gsub("^female$","F",gender_kid))

mom_child_id<-read_excel("__PATH_003__",sheet = 2)

bmi_child<-read_excel("__PATH_004__") %>%
  mutate(BMI=weight/(height/100)^2,date=as.Date(date),birth=as.Date(birth)) %>%
  dplyr::rename_with(~gsub("^","kid_",.x)) %>%
  dplyr::rename(sjid_kid=kid_IDchild)

coldata_349<-read.csv("__PATH_005__",header = T,row.names = 1) %>%
  dplyr::rename(sjid_mo=ID) %>%
  mutate(sjid_mo=as.character(sjid_mo),
         period=gsub("1","early",trimester),
         period=gsub("2","middle",period),
         period=gsub("3","late",period),
         sign=paste(sjid_mo,period,sep = "")) %>%
  dplyr::select(-period) %>%
  left_join(mom_child_id[,1:2]) %>%
  left_join(bmi_child)


coldata_kid_bmiz<-coldata_349 %>%
  left_join(gender) %>%
  #dplyr::select(-c(sample,trimester)) %>%
  unique() %>%
  na.omit() %>%
  mutate(disease_gender=paste(disease,gender_kid,sep = "_"),
         kid_date=as.Date(kid_date),kid_birth=as.Date(kid_birth),
         kid_date_month=as.numeric(kid_date-kid_birth)/30.4375) %>%
  mutate(
    bmi_z_kid = case_when(
      # 如果年龄 <= 60 个月，使用 anthro 包的 antho_zscores 函数
      kid_date_month <= 60 ~ anthro_zscores(
        sex = gender_kid, 
        age = kid_date_month, is_age_in_month = T, 
        weight = kid_weight, 
        lenhei = kid_height
      )$zbmi,
      
      # 如果年龄 > 60 个月，使用 anthroplus 包的 anthroplus_zscores 函数
      kid_date_month > 60 ~ anthroplus_zscores(
        sex = tolower(gender_kid), 
        age_in_months = kid_date_month,
        weight_in_kg = kid_weight, 
        height_in_cm = kid_height
      )$zbfa
    ))

coldata_1476_bmizkid<-coldata_kid_bmiz %>%
  dplyr::rename(ID=sjid_mo) %>%
  dplyr::select(names(coldata_1476),sjid_kid,kid_birth,kid_date_month,kid_age,bmi_z_kid) %>%
  dplyr::rename(sjid_mo=ID) %>%
  dplyr::filter(kid_date_month>=24)

save(coldata_1476_bmizkid,limma_voom_module13,file = "test_bmi_z_kid_prediction_260320.rda")
