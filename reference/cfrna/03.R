library(tidyverse)
library(limma)
library(openxlsx)
library(edgeR)
library(DESeq2)

#setwd("__PATH_001__")
#count_1562<-read.csv("cfRNA_count_1562_0603.csv",header = T,row.names = 1)
#count_1562<-round(count_1562)
coldata_1562<-read.csv("__PATH_002__",row.names = 1) %>%
  mutate(ID=as.character(ID),trimester=gsub("1","1st",trimester),trimester=gsub("2","2nd",trimester),trimester=gsub("3","3rd",trimester)) %>%
  mutate(batch=as.factor(batch)) %>%
  mutate(BMI_mo_Anthropometric=as.numeric(BMI_mo_Anthropometric))

count_lncRNA_1562<-read.table("__PATH_003__",row.names = 1,header = T)[,coldata_1562$sample]
count_lncRNA_1562<-count_lncRNA_1562[apply(count_lncRNA_1562,1,function(x) sum(x>1) > floor(ncol(count_lncRNA_1562)/10)),]

# BMI分层：正常 18.5–23.9，超重 24–27.9，肥胖 >=28
coldata_1562$BMI_group <- cut(
  coldata_1562$BMI_mo_Anthropometric,
  breaks = c(18.5, 24, 28, Inf),
  right = FALSE,
  labels = c("normal", "overweight", "obese")
)

# 仅保留可归入三个 BMI 组的样本，并按 BMI 组生成 coldata
coldata_bmi_list <- split(
  coldata_1562[!is.na(coldata_1562$BMI_group), , drop = FALSE],
  coldata_1562$BMI_group[!is.na(coldata_1562$BMI_group)],
  drop = TRUE
)

for (bmi_name in names(coldata_bmi_list)) {
  print(bmi_name)
  coldata_bmi <- coldata_bmi_list[[bmi_name]]
  count_lncRNA_bmi <- count_lncRNA_1562[, coldata_bmi$sample, drop = FALSE]

  #dds <- DESeqDataSetFromMatrix(countData=count_lncRNA_bmi, colData=coldata_bmi, design=~ disease + trimester + disease:trimester + batch)
  #dds <- DESeq(dds,quiet=TRUE)
  design_test<-model.matrix(~trimester*disease+batch,data = coldata_bmi) # set coefficient for linear model
  colnames(design_test) = make.names(colnames(design_test))
  #cont.matrix<-cbind(de_1=c(0,0,0,1,0,0,0,0,0,0),de_2=c(0,0,0,1,0,0,0,0,1,0),
  #                   de_3=c(0,0,0,1,0,0,0,0,0,1),de_all=c(0,0,0,1,0,0,0,0,1/3,1/3))

  cont.matrix<-cbind(de_1=c(0,0,0,1,0,0,0),de_2=c(0,0,0,1,0,1,0),
                     de_3=c(0,0,0,1,0,1,1),de_all=c(0,0,0,1,0,1/3,1/3))

  limma_matrix<-count_lncRNA_bmi %>%
    DGEList() # filter genes with low counts
  print(paste0(bmi_name,": 1"))
  limma_matrix<-calcNormFactors(limma_matrix,method = "TMM")
  print(paste0(bmi_name,": 2"))
  limma_res_weight_tmp<-voomWithQualityWeights(limma_matrix, design_test,normalize.method = "none",method = "genebygene",trace=T)
  print(paste0(bmi_name,": 3"))
  corfit_weight<-duplicateCorrelation(limma_res_weight_tmp,design_test,block=coldata_bmi$ID)
  print(paste0(bmi_name,": 4"))
  limma_res_weight<-voomWithQualityWeights(limma_matrix, design_test,normalize.method = "none",method = "genebygene",correlation =corfit_weight$consensus,trace=T)
  print(paste0(bmi_name,": 5"))
  fit <- lmFit(limma_res_weight, design_test, block=coldata_bmi$ID, correlation=corfit_weight$consensus)
  print(paste0(bmi_name,": 6"))
  fit2<- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  de_1<-topTable(fit2,coef="de_1",adjust.method="BH",n=Inf) %>% mutate(period="1",SYMBOL=rownames(.),BMI_group=bmi_name)
  de_1<-de_1[,c(8,1:7)]
  de_2<-topTable(fit2,coef="de_2",adjust.method="BH",n=Inf) %>% mutate(period="2",SYMBOL=rownames(.),BMI_group=bmi_name)
  de_2<-de_2[,c(8,1:7)]
  de_3<-topTable(fit2,coef="de_3",adjust.method="BH",n=Inf) %>% mutate(period="3",SYMBOL=rownames(.),BMI_group=bmi_name)
  de_3<-de_3[,c(8,1:7)]
  de_all<-topTable(fit2,coef="de_all",adjust.method="BH",n=Inf) %>% mutate(period="all",SYMBOL=rownames(.),BMI_group=bmi_name)
  de_all<-de_all[,c(8,1:7)]

  de_merge <- createWorkbook()
  addWorksheet(de_merge, sheetName = "deg_1")
  writeData(de_merge, sheet = "deg_1", x = de_1)
  addWorksheet(de_merge, sheetName = "deg_2")
  writeData(de_merge, sheet = "deg_2", x = de_2)
  addWorksheet(de_merge, sheetName = "deg_3")
  writeData(de_merge, sheet = "deg_3", x = de_3)
  addWorksheet(de_merge, sheetName = "deg_all")
  writeData(de_merge, sheet = "deg_all", x = de_all)
  saveWorkbook(de_merge, file = paste0("limma_349_lncRNA_BMI_", bmi_name, ".xlsx"), overwrite = TRUE)

  de_merge<-rbind(de_1,de_2,de_3,de_all)
  write.csv(de_merge,file=paste0("limma_349_lncRNA_BMI_",bmi_name,".csv"),row.names = F,quote = F)
}
