library(tidyverse)
library(readxl)
library(WGCNA)
#library(sva)
library(openxlsx)

load("__PATH_001__")

vst_cfRNA<-read.csv("__PATH_002__",row.names = 1)

sample_meta<-read_excel("__PATH_003__",sheet = 10) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_meta=sampleid)
sample_cfRNA<-read_excel("__PATH_003__",sheet = 8) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_cfRNA=sampleid) %>%
  left_join(sample_meta[,c("sign","sampleid_meta")])


##### meta pca #####
load("../代谢组/wgcna/wgcna merge-cfRNA &  metabolite/coldata_wgcna.rda")
rpm_meta<-read.csv("__PATH_004__",header = T,row.names = 1)

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign))

coldata_wgcna<-coldata_wgcna %>%
  dplyr::filter(rownames(.) %in% intersect(rownames(.),names(rpm_meta)))

rpm_meta<-rpm_meta[,rownames(coldata_wgcna)]
rpm_meta<-log2(rpm_meta+1)

dat <- as.data.frame(t(rpm_meta))
dat.pca <- PCA(dat, graph = FALSE)
group_list_datail<-factor(coldata_wgcna$disease)
ggord(dat.pca,group_list_datail,coord_fix = FALSE,ellipse=F,
      xlims=c(-270,270),ylims=c(-270,270),
      arrow = 0,vec_ext = 0,txt = NULL)

##### add meta ko #####
load("__PATH_005__")
rpm_meta<-read.csv("__PATH_004__",header = T,row.names = 1)

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  t() %>% as.data.frame()

## log2 trans
#rpm_meta<-log2(rpm_meta+1)

#rpm_meta<-rpm_meta %>%
#  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA[which(sample_cfRNA$sampleid_meta %in% colnames(.))],sample_cfRNA$sampleid_meta[which(sample_cfRNA$sampleid_meta %in% colnames(.))])) %>%
#  t() %>% as.data.frame()

rpm_meta<-na.omit(rpm_meta[rownames(coldata_wgcna),])
missing_sample <- setdiff(rownames(coldata_wgcna), rownames(rpm_meta))
rpm_meta[missing_sample, ] <- 0
rpm_meta<-rpm_meta[rownames(coldata_wgcna),]
rpm_meta<-as.data.frame(log2(rpm_meta+1))

all(rownames(rpm_meta)==rownames(coldata_wgcna))

coldata_wgcna_ko<-cbind(coldata_wgcna,rpm_meta)




#meta_ko<-read.csv("__PATH_006__") %>%
#  dplyr::rename(K_name=feature) %>%
#  left_join(unique(kegg_ko_gmt_micro[,c("V1","K_name")])) %>%
#  dplyr::rename(name=V1) %>%
#  dplyr::select(1:2,13,3:12) %>%
#  mutate(name=gsub("D.K..... ","",name))



##### meta 菌种 #####
load("__PATH_005__")

rpm_meta<-read.csv("__PATH_007__",header = T,row.names = 1)

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  t() %>% as.data.frame()

#rpm_meta<-log2(rpm_meta+1)

#rpm_meta<-rpm_meta %>%
#  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_cfRNA[which(sample_cfRNA$sampleid_meta %in% colnames(.))],sample_cfRNA$sampleid_meta[which(sample_cfRNA$sampleid_meta %in% colnames(.))])) %>%
#  t() %>% as.data.frame()

rpm_meta<-na.omit(rpm_meta[rownames(coldata_wgcna),])
missing_sample <- setdiff(rownames(coldata_wgcna), rownames(rpm_meta))
rpm_meta[missing_sample, ] <- 0
rpm_meta<-rpm_meta[rownames(coldata_wgcna),]
rpm_meta<-as.data.frame(log2(rpm_meta+1))

all(rownames(rpm_meta)==rownames(coldata_wgcna))

coldata_wgcna_strain<-cbind(coldata_wgcna,rpm_meta)



save(vst_cfRNA,coldata_wgcna_ko,coldata_wgcna_strain,
     file="__PATH_008__")

##### meta WGCNA seperate #####
rpm_meta<-read.csv("__PATH_007__",header = T,row.names = 1)

sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))

rpm_meta<-rpm_meta[apply(rpm_meta,1,function(x) sum(x>0) > floor(ncol(rpm_meta)/10)),]

rpm_meta_flt<-rpm_meta %>%
  #dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  t() %>% as.data.frame()

#rpm_meta_flt<-as.data.frame(log2(rpm_meta_flt+1))

##sample filter
sampleTree = hclust(dist(rpm_meta_flt), method = "average")
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
rpm_meta_flt <- rpm_meta_flt[keepSamples, ] 

##sft
powers = c(c(1:10), seq(from = 12, to=20, by=1))
sft = pickSoftThreshold(rpm_meta_flt, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
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
