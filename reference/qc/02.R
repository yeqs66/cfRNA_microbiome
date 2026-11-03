library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library("vegan")
library(FactoMineR)
library(ggord)
library(patchwork)
library(readxl)

master_key<-read.csv("master_key.csv")
matrix_meta<-read.csv("microbiome ko gsea复现/new_deg/TJGDM1962s_ko_tri_95_rpk.csv",header = T,sep = ",",row.names = 1)
load("microbiome ko gsea复现/microbiome_sampleinfo_GDM1999s_20240318.rda")

#sequence_info_meta<-read.csv("microbiome ko gsea复现/TJ_meta.QC.all.update20260809.csv") %>%
#  mutate(QC__RawReadsCount=as.numeric(QC__RawReadsCount),
#         QC__HostRate=as.numeric(QC__HostRate),
#         QC__GCContentPercent=as.numeric(QC__GCContentPercent),
#         QC__Q20Percent=as.numeric(QC__Q20Percent),
#         QC__Q30Percent=as.numeric(QC__Q30Percent))
#sequence_info_meta<-sequence_info_meta[match(colnames(matrix_meta), sequence_info_meta$sampleID),]

sequence_info_meta<-read_excel("microbiome ko gsea复现/1962_meta_sampleid_sjidmo_addTASK20260810.xlsx") %>%
  dplyr::select("样品名称","任务单号") %>%
  dplyr::rename(sampleID="样品名称",TASK="任务单号") %>%
  unique() %>%
  left_join(microbiome_sampleinfo[,c("sample","group","trimester")], by=c("sampleID"="sample"))


## profile of sequencing depth & GC & batch
sequence_info_meta_tmp<-read_excel("microbiome ko gsea复现/1962_meta_sampleid_sjidmo_addTASK20260810.xlsx") %>%
  dplyr::filter(样品名称 %in% colnames(matrix_meta),状态=="已交付") %>%
  dplyr::select("样品名称","任务单号","GC%","ReadNum","BaseNum",`>Q30%(total)`) %>%
  mutate(
    `GC%` = as.numeric(`GC%`),
    ReadNum = as.numeric(ReadNum),
    BaseNum = as.numeric(BaseNum),
    `>Q30%(total)` = as.numeric(`>Q30%(total)`)
  ) %>%
  unique() %>%
  group_by(样品名称) %>%
  summarise(
    任务单号 = first(任务单号),                 # 同一样本任务单号取第一个
    `GC%`    = mean(`GC%`,    na.rm = TRUE),    # 同名样本 GC% 取平均
    ReadNum  = mean(ReadNum,  na.rm = TRUE),    # ReadNum 取平均
    BaseNum  = mean(BaseNum,  na.rm = TRUE),    # BaseNum 取平均
    `>Q30%(total)` = mean(`>Q30%(total)`, na.rm = TRUE), # >Q30%(total) 取平均
    .groups  = "drop"
  )

summary(sequence_info_meta_tmp$`GC%`)
summary(sequence_info_meta_tmp$ReadNum)
summary(sequence_info_meta_tmp$BaseNum)
summary(sequence_info_meta_tmp$`>Q30%(total)`)


## PERMANOVA
sequence_info_meta_tmp<-sequence_info_meta %>%
  dplyr::filter(!is.na(QC__RawReadsCount))
data<-matrix_meta[,sequence_info_meta_tmp$sampleID]
dist_mat <- vegdist(t(data), method = "euclidean")
adonis2(dist_mat ~ QC__RawReadsCount+QC__HostRate+QC__GCContentPercent+QC__Q20Percent+QC__Q30Percent, data = data.frame(sequence_info_meta_tmp), method = "euclidean",permutations = 999)

ord <- metaMDS(dist_mat, k = 2)
sc <- as.data.frame(scores(ord))
sc$cov <- sequence_info_meta_tmp$covar[rownames(sc)]

ggplot(sc, aes(x = NMDS1, y = NMDS2, color = cov)) +
  geom_point(size = 2) +
  scale_color_gradientn(colors = c("navy","cyan","yellow","red"),
                        values = seq(0,1,length.out=4)) +
  stat_ellipse(aes(group = cut(cov,3)), linetype = 2) +
  theme_minimal()


## PCA
#color_rawreadcount<-
exp <- log2(matrix_meta+1)
dat <- as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
group_list_group<-factor(sequence_info_meta$group)
group_list_batch<-factor(sequence_info_meta$TASK)
#group_list_ethnicity<-as.factor(sequence_info_meta$ethnicity)
group_list_trimester<-as.factor(sequence_info_meta$trimester)
col_group<-c("GDM"="#E18E6D","Control"="#62B197")
col_batch<-c("W220914003"="#2f5688","W220512037"="#CC0000")
#col_ethnicity<-c("Han"="#E18E6D","minority"="#62B197")
col_trimester<-c("early"="#AED185","middle"="#F89FA8","late"="#AECDD7")

test<-data.frame(sample=rownames(dat.pca$ind$coord),pc1=dat.pca$ind$coord[,1],pc2=dat.pca$ind$coord[,2],
                 Batch=as.factor(sequence_info_meta$TASK),Group=as.factor(sequence_info_meta$group),Trimester=as.factor(sequence_info_meta$trimester))

#p_microbiome_qc<-list()
#for (i in c("QC__RawReadsCount","QC__HostRate","QC__GCContentPercent","QC__Q20Percent","QC__Q30Percent")) {
#  p_microbiome_qc[[i]]<-ggplot(coord, aes_string(x = "Dim.1", y = "Dim.2", color = i)) +
#    geom_point(size = 2) +
#    scale_color_viridis_c() +            # 连续渐变色
#    labs(color = i) +
#    theme_minimal()
#}

pc1_var <- dat.pca$eig[1, "percentage of variance"]
pc2_var <- dat.pca$eig[2, "percentage of variance"]

p<-ggplot(test,aes(x=pc1,y=pc2,color=Batch,shape=Batch))+
  geom_point(size=2,alpha=0.8)+
  stat_ellipse(aes(color=Batch),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_task)+
  scale_shape_manual(values = c("W220914003"=16,"W220512037"=18))+
  #scale_color_manual(values = c(NA,NA))+
  xlab(sprintf("PC1 (%.2f%%)", pc1_var))+
  ylab(sprintf("PC2 (%.2f%%)", pc2_var))+
  ggtitle("microbiome KO")+
  #labs(fill = "Batch")+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = 'gray'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(color="black",face="bold"),
        axis.title = element_text(color = "black",face="bold"),
        plot.title = element_text(color="black",face="bold",size=15,hjust = 0.5),
        legend.text = element_text(face="bold")
  )





