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
  mutate(batch=as.factor(batch))
count_lncRNA_1562<-read.table("__PATH_003__",row.names = 1,header = T)[,coldata_1562$sample]
count_lncRNA_1562<-count_lncRNA_1562[apply(count_lncRNA_1562,1,function(x) sum(x>1) > floor(ncol(count_lncRNA_1562)/10)),]

#dds <-  DESeqDataSetFromMatrix(countData=count_1562, colData=coldata_1562, design= ~ disease + trimester + disease:trimester + batch + smoking_mo + drinking_mo_18topreg + ethnicity_mo)
#dds <- DESeq(dds,quiet=TRUE)
design_test<-model.matrix(~trimester*disease+batch,data = coldata_1562) # set coefficient for linear model
colnames(design_test) = make.names(colnames(design_test))

#cont.matrix<-cbind(de_1=c(0,0,0,1,0,0,0,0,0,0),de_2=c(0,0,0,1,0,0,0,0,1,0),
#                   de_3=c(0,0,0,1,0,0,0,0,0,1),de_all=c(0,0,0,1,0,0,0,0,1/3,1/3))

cont.matrix<-cbind(de_1=c(0,0,0,1,0,0,0),de_2=c(0,0,0,1,0,1,0),
                   de_3=c(0,0,0,1,0,1,1),de_all=c(0,0,0,1,0,1/3,1/3))


limma_matrix<-count_lncRNA_1562 %>%
  DGEList() # filter genes with low counts
print("1")
limma_matrix<-calcNormFactors(limma_matrix,method = "TMM")
print("2")
limma_res_weight_tmp<-voomWithQualityWeights(limma_matrix, design_test,normalize.method = "none",method = "genebygene",trace=T)
print("3")
corfit_weight<-duplicateCorrelation(limma_res_weight_tmp,design_test,block=coldata_1562$ID)
print("4")
limma_res_weight<-voomWithQualityWeights(limma_matrix, design_test,normalize.method = "none",method = "genebygene",correlation =corfit_weight$consensus,trace=T)
print("5")
fit <- lmFit(limma_res_weight, design_test, block=coldata_1562$ID, correlation=corfit_weight$consensus)
print("6")
save.image("limma_lncRNA.RData")
fit2<- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

save.image()

de_1<-topTable(fit2,coef="de_1",adjust.method="BH",n=Inf) %>% mutate(period="1",SYMBOL=rownames(.))
de_1<-de_1[,c(8,1:7)]
de_2<-topTable(fit2,coef="de_2",adjust.method="BH",n=Inf) %>% mutate(period="2",SYMBOL=rownames(.))
de_2<-de_2[,c(8,1:7)]
de_3<-topTable(fit2,coef="de_3",adjust.method="BH",n=Inf) %>% mutate(period="3",SYMBOL=rownames(.))
de_3<-de_3[,c(8,1:7)]
de_all<-topTable(fit2,coef="de_all",adjust.method="BH",n=Inf) %>% mutate(period="all",SYMBOL=rownames(.))
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
saveWorkbook(de_merge, file = "limma_349_lncRNA_batch_only.xlsx", overwrite = TRUE)
