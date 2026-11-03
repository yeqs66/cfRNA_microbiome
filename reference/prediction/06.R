library(tidyverse)
library(readxl)
library(openxlsx)
library(ggplot2)
library(splines)
library(anthro)
library(anthroplus)
library(lme4)
library(lmerTest)
library(ggsignif)
library(ggstats)
library(ggpubr)

rm(list=ls())

load("__PATH_001__")
load("__PATH_002__")

trio_ped<-read.table("__PATH_003__") %>%
  dplyr::rename_with(~c("FID","KID","DID","MID","V1","V2"))
gender<-labelledclinical[,c("sjid_kid","gender_kid")] %>%
  mutate(gender_kid=gsub("^male$","M",gender_kid),
         gender_kid=gsub("^female$","F",gender_kid))

mom_child_id<-read_excel("__PATH_004__",sheet = 2)

bmi_child<-read_excel("__PATH_005__") %>%
  mutate(BMI=weight/(height/100)^2,date=as.Date(date),birth=as.Date(birth)) %>%
  dplyr::rename_with(~gsub("^","kid_",.x)) %>%
  dplyr::rename(sjid_kid=kid_IDchild)

coldata_349<-read.csv("__PATH_006__",header = T,row.names = 1) %>%
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

  

#sp_pred <- smooth.spline(tmp$kid_age, tmp$kid_BMI)
#sp_df <- data.frame(
#  x = seq(min(tmp$kid_age), max(tmp$kid_age), length.out = 100),
#  y = predict(sp_pred, x = seq(min(tmp$kid_age), max(tmp$kid_age), length.out = 100))$y
#)

ggplot(data = coldata_kid_bmiz,aes(x=kid_age,y=bmi_z_kid,color=disease_gender))+
  geom_point(alpha=0.2)+
  #geom_line(data = sp_df, aes(x = x, y = y), color = "red", size = 1)+
  stat_smooth(method = "lm",
              formula = y~ns(x,df=3),
              aes(color = disease_gender),
              #fill = "lightblue",    # 置信区间填充色
              size = 1.5,
              level = 0.90)

#coldata_349<-coldata_349 %>%
#  left_join(coldata_kid_bmiz[,c("sjid_kid","bmi_z","kid_date_month","disease_gender")]) 
  

##### kid_age=0 #####
comparisons <- combn(na.omit(unique(coldata_kid_bmiz$disease_gender)), 2, simplify = FALSE)

ggplot(coldata_kid_bmiz,aes(x=disease_gender,y=bmi_z_kid,color=disease_gender))+ 
  geom_boxplot()+
  geom_signif(comparisons = comparisons, 
              y_position = c(6, 7, 8, 9,10,11,12),
              test = "t.test")+
  theme_bw()





##### LMM cfRNA #####
load("__PATH_007__")
load("__PATH_008__")
deg_40<-read_excel("manu_v1.1/Supplementary Table 1_Differential_analysis_results_for_cf-mRNA.xlsx") %>%
  dplyr::filter(adj.P.Val<0.1)
deg_40<-deg_40$SYMBOL

tmp<-limma_voom_mRNA %>%
  dplyr::filter(rownames(.) %in% deg_40) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var="sample")
  
coldata_349_cfRNA<-coldata_kid_bmiz %>%
  left_join(tmp)

res_cfRNA_bmi<-data.frame()
for (gene in deg_40) {
  res_tmp<-data.frame()
  for (period in unique(coldata_349_cfRNA$trimester)) {
    tmp<-coldata_349_cfRNA %>%
      dplyr::filter(trimester==period,kid_date_month>=24)
    
    formula_str <- paste0(
      #"bmi_z_kid ~ ns(kid_date_month, df = 3) * ", gene, 
      "bmi_z_kid ~ kid_date_month + ", gene, 
      " + disease + gender_kid + BMI_mo_Anthropometric + (1 | sjid_kid)"
    )
    
    model <- lmerTest::lmer(
      formula = as.formula(formula_str), 
      data = tmp,
      REML = TRUE
    )
    tmp<-as.data.frame(anova(model)) %>%
      rownames_to_column(var = "coef") %>%
      mutate(feature=gene,trimester=period)
    res_tmp<-rbind(res_tmp,tmp)
  }
  res_cfRNA_bmi<-rbind(res_cfRNA_bmi,res_tmp)
}

tmp<-res_cfRNA_bmi %>%
  dplyr::filter(str_detect(coef,"\\:"),
                trimester==1)
tmp$padj<-p.adjust(tmp$`Pr(>F)`,method = "BH")


##kid_age=0 & cfRNA
tmp<-coldata_349_cfRNA %>%
  dplyr::filter(kid_age==0)



##### LMM cfRNA module 7 #####
load("wgcna/cfRNA_microbiome/m_mi_lnc/ko/wgcga_signed_signed_hybrid_pearson/MEs.rda")
gene_in_module<-read.csv("wgcna/cfRNA_microbiome/m_mi_lnc/ko/wgcga_signed_signed_hybrid_pearson/module_gene-signed-signed_hybrid-pearson.csv") %>%
  dplyr::select(-X) %>% unique()
MEs<-MEs %>%
  rownames_to_column(var = "sign")

coldata_349_wgcna<-coldata_kid_bmiz %>%
  left_join(MEs)

##kid_age=0
tmp<-coldata_349_wgcna %>%
  dplyr::filter(kid_age==0) %>%
  na.omit() %>%
  unique()

ggplot(tmp,aes(x=bmi_z_kid,y=MEtan,color=disease_gender))+
  geom_point()+
  geom_smooth(method = "lm",
              #formula = y~ns(x,df=3),
              aes(color = disease_gender),
              #fill = "lightblue",    # 置信区间填充色
              size = 1.5,
              level = 0.90)

##### wgcna module cor #####
tmp<-coldata_349_wgcna %>%
  dplyr::filter(kid_date_month>=24) %>%
  mutate(period=gsub("1","early",trimester),
         period=gsub("2","middle",period),
         period=gsub("3","late",period)) %>%
  dplyr::select(sample,gender_kid,disease,kid_date_month,period,MEtan,MEred,MEcyan,MElightyellow,MEturquoise,MEsalmon,MEblack,MEmidnightblue,MEmagenta,MEblue,MEgreenyellow,MElightcyan,MEgreen,MEpink,MElightgreen,MEpurple,MEyellow,MEbrown,MEgrey60,MEgrey,bmi_z_kid,kid_date_month) %>%
  pivot_wider(
    # 如果每个性别-疾病组合只有一行，直接用 bmi_z_kid 作为标识列
    # 如果有多行（重复测量），需要添加 sample_id 作为 id_cols
    id_cols = c(sample,bmi_z_kid,kid_date_month),  
    names_from = c(gender_kid, disease,period),
    values_from = c(MEtan,MEred,MEcyan,MElightyellow,MEturquoise,MEsalmon,MEblack,MEmidnightblue,MEmagenta,MEblue,MEgreenyellow,MElightcyan,MEgreen,MEpink,MElightgreen,MEpurple,MEyellow,MEbrown,MEgrey60,MEgrey),
    names_glue = "{.value}_{gender_kid}_{disease}_{period}",  # 格式：MEtan_F_GDM
    values_fill = NA
  )

cor_wgcna_bmi<-data.frame(matrix(nrow = 1,ncol = 3))
for (me in names(tmp)[6:25]) {
  cor_p<-cor.test(tmp[[me]],tmp$bmi_z_kid,alternative="two.sided",method="pearson")
  cor<-cor_p$estimate
  cor_p<-cor_p$p.value
  cor_wgcna_bmi<-rbind(cor_wgcna_bmi,c(cor,cor_p,me))
}
cor_wgcna_bmi<-cor_wgcna_bmi %>%
  mutate(X1=as.numeric(X1),
         X2=as.numeric(X2)) %>%
  na.omit() %>%
  dplyr::rename_with(~c("cor","cor_p","ME")) %>%
  mutate(cor_padj=p.adjust(cor_p,method = "BH"))





##### LLM for wgcna #####
res_wgcna7_bmi<-data.frame()
for (me in colnames(MEs)[-1]) {
  res_tmp<-data.frame()
  #for (period in unique(coldata_349_wgcna$trimester)) {
    tmp<-coldata_349_wgcna %>%
      dplyr::filter(
        #trimester==period,
        kid_date_month>=24
        )
    
    formula_str <- paste0(
      "bmi_z_kid ~ kid_date_month + ", me, 
      " * disease + gender_kid + BMI_mo_Anthropometric + (1 | sjid_kid)"
    )
    
    model <- lmerTest::lmer(
      formula = as.formula(formula_str), 
      data = tmp,
      REML = TRUE
    )
    tmp<-summary(model)
    tmp<-as.data.frame(tmp$coefficients) %>%
      rownames_to_column(var = "coef") %>%
      mutate(
        feature=me,
        #trimester=period
        )
    res_tmp<-rbind(res_tmp,tmp)
  #}
  res_wgcna7_bmi<-rbind(res_wgcna7_bmi,res_tmp)
}


##### cor MOFA2 #####
model<-load_model("MOFA2/model_tpm_ko_topVarience50_260120.hdf5")
load("mofa_input_ko_m_mi_lncRNA_topVarience50_260120.rda")
coldata_mofa<-coldata %>%
  mutate(sample=rownames(.))
samples_metadata(model) <- coldata_mofa %>%
  dplyr::select(sample,disease,period,ALT_Liver,AST_Liver,Preterm_OffspringRelated) %>%
  mutate(group=period) %>%   ##group by period
  mutate(disease=gsub("0","control",disease),
         disease=gsub("1","GDM",disease))

factors <- as.data.frame(get_factors(model, factors = "all",groups = 1,as.data.frame = F, scale = T)) %>%   dplyr::rename_with(~gsub("^X1.","",.)) %>%
  rownames_to_column(var = "sign")


tmp<-coldata_kid_bmiz %>%
  dplyr::filter(kid_date_month>=24) %>%
  left_join(factors) %>%
  mutate(period=gsub("1","early",trimester),
         period=gsub("2","middle",period),
         period=gsub("3","late",period)) %>%
  dplyr::select(sample,gender_kid,disease,kid_date_month,period,Factor1,Factor2,Factor3,Factor4,Factor5,bmi_z_kid,kid_date_month) %>%
  pivot_wider(
    # 如果每个性别-疾病组合只有一行，直接用 bmi_z_kid 作为标识列
    # 如果有多行（重复测量），需要添加 sample_id 作为 id_cols
    id_cols = c(sample,bmi_z_kid,kid_date_month),  
    names_from = c(gender_kid, disease,period),
    values_from = c(Factor1,Factor2,Factor3,Factor4,Factor5),
    names_glue = "{.value}_{gender_kid}_{disease}_{period}",  # 格式：MEtan_F_GDM
    values_fill = NA
  )

cor_mofa_bmi<-data.frame(matrix(nrow = 1,ncol = 3))
for (me in names(tmp)[6:10]) {
  cor_p<-cor.test(tmp[[me]],tmp$bmi_z_kid,alternative="two.sided",method="pearson")
  cor<-cor_p$estimate
  cor_p<-cor_p$p.value
  cor_mofa_bmi<-rbind(cor_mofa_bmi,c(cor,cor_p,me))
}
cor_mofa_bmi<-cor_mofa_bmi %>%
  mutate(X1=as.numeric(X1),
         X2=as.numeric(X2)) %>%
  na.omit() %>%
  dplyr::rename_with(~c("cor","cor_p","Factor")) %>%
  mutate(cor_padj=p.adjust(cor_p,method = "BH"))











