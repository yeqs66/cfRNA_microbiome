library(openxlsx)
library(readxl)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggbeeswarm)
library(ggseqlogo)

rm(list = ls())

##### stat #####
stat_149<-read_excel("__PATH_001__",sheet = 1)
stat_200<-read_excel("__PATH_001__",sheet = 2)
stat_9<-read_excel("__PATH_001__",sheet = 3)
stat_1980<-read.

gene_dect<-read_excel("__PATH_002__",sheet = 1) %>%
  dplyr::rename(Sample=...1,'mRNA(Count>0)'=mRNA,'lncRNA(Count>0)'=lncRNA,'miRNA(Count>0)'=miRNA,
                'tRNA(Count>0)'=tRNA,'piRNA(Count>0)'=piRNA)

stat_149<-stat_149 %>%
  dplyr::filter(!Sample %in% stat_9$Sample) %>%
  rbind(stat_9)

stat_all<-rbind(stat_149,stat_200) %>%
  #na.omit() %>%
  dplyr::select(-Uniq_Genome,-Mult_Genome) %>%
  dplyr::rename(first_clean_reads=clean_reads) %>%
  mutate(first_filter= (Total - first_clean_reads) / Total * 100) %>%
  mutate(rRNA=as.numeric(rRNA),miRNA=as.numeric(miRNA),tpiRNA=as.numeric(tpiRNA),
         mRNA=as.numeric(mRNA),otherRNA=as.numeric(otherRNA),Sec_flt_reads=as.numeric(Sec_flt_reads)) %>%
  mutate(second_filter=Sec_flt_reads / first_clean_reads * 100, 
         final_clean_reads=first_clean_reads - Sec_flt_reads,
         final_filter= (Total - final_clean_reads) / Total * 100) %>%
  mutate("r/Y/vtRNA"=rRNA / final_clean_reads * 100, miRNA=miRNA / final_clean_reads * 100, 
         "t/piRNA"=tpiRNA / final_clean_reads * 100, "m/lncRNA"=mRNA / final_clean_reads * 100,
         otherRNA=otherRNA / final_clean_reads * 100) %>%
  dplyr::select(-c(rRNA,tpiRNA,mRNA)) %>%
  left_join(gene_dect,by="Sample")

name_order_stat<-names(stat_all)
name_order_stat<-name_order_stat[c(1,2,4,9,5,10,11,12,13,7,14,15,16:20)]

stat_all<-stat_all %>%
  dplyr::select(name_order_stat)

##### pca-cfRNA mRNA #####
coldata_all<-read.csv("__PATH_003__",row.names = 1)
#coldata_all<-read.csv("coldata")
#tpm_1562<-read.csv("../cfRNA_tpm_1562.csv",row.names = 1) %>%
#  dplyr::select(coldata_all$sample)

tpm_1562<-read.csv("__PATH_004__",row.names = 1,header = T) %>%
  dplyr::select(coldata_all$sample)

library(ggord)
library(FactoMineR)

exp <- log2(tpm_1562+1)
dat <- as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
group_list_group<-factor(coldata_all$GDMgroup2023)
group_list_batch<-factor(coldata_all$batch)
group_list_ethnicity<-as.factor(coldata_all$ethnicity)
group_list_trimester<-as.factor(coldata_all$trimester)
col_group<-c("GDM"="#E18E6D","control"="#62B197")
col_batch<-c("149"="#2f5688","200"="#CC0000")
col_ethnicity<-c("Han"="#E18E6D","minority"="#62B197")
col_trimester<-c("1"="#AED185","2"="#F89FA8","3"="#AECDD7")
#shape_group<-c("GDM"="#E18E6D","control"="#62B197")
#shape_batch<-c("149"="#2f5688","200"="#CC0000")

test<-data.frame(sample=rownames(dat.pca$ind$coord),pc1=dat.pca$ind$coord[,1],pc2=dat.pca$ind$coord[,2],
                 Batch=as.factor(coldata_all$batch),Group=as.factor(coldata_all$disease),Trimester=as.factor(coldata_all$trimester))
#test_mean<-aggregate(cbind(pc1,pc2)~Batch,data=test,FUN=mean)
#test<-merge(test,setNames(test_mean,c("Batch","pc1_mean","pc2_mean")))

p_batch<-ggplot(test,aes(x=pc1,y=pc2,color=Batch,shape=Batch))+
  geom_point(size=2,alpha=0.8)+
  stat_ellipse(aes(color=Batch),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_batch)+
  scale_shape_manual(values = c("149"=16,"200"=18))+
  #scale_color_manual(values = c(NA,NA))+
  xlab("P1 (31.83%)")+
  ylab("P2 (8.74%)")+
  ggtitle("mRNA")+
  labs(fill = "Batch")+
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

p_trimester<-ggplot(test,aes(x=pc1,y=pc2,color=Trimester,shape=Trimester))+
  geom_point(size=2,alpha=1)+
  stat_ellipse(aes(color=Trimester),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_trimester)+
  scale_shape_manual(values = c("1"=16,"2"=15,"3"=18))+
  #scale_color_manual(values = c(NA,NA,NA))+
  xlab("P1 (31.83%)")+
  ylab("P2 (8.74%)")+
  ggtitle("mRNA")+
  labs(fill = "Trimester")+
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

p_group<-ggplot(test,aes(x=pc1,y=pc2,color=Group,shape=Group))+
  geom_point(size=2,alpha=1)+
  stat_ellipse(aes(color=Group),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_group)+
  scale_shape_manual(values = c("GDM"=16,"control"=18))+
  xlab("P1 (31.83%)")+
  ylab("P2 (8.74%)")+
  ggtitle("mRNA")+
  labs(fill = "Group")+
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


p_mRNA<-p_group+p_batch+p_trimester

##### pca-cfRNA lncRNA #####
coldata_all<-read.csv("__PATH_003__",row.names = 1)
#tpm_1562<-read.csv("../cfRNA_tpm_1562.csv",row.names = 1) %>%
#  dplyr::select(coldata_all$sample)

#tpm_1562<-read.table("__PATH_005__",row.names = 1,header = T) %>%
#  dplyr::select(coldata_all$sample)
load("__PATH_006__")
tpm_lnc_1562_6473<-tpm_lnc_1562_6473 %>%
  dplyr::select(coldata_all$sample)

library(ggord)
library(FactoMineR)

exp <- log2(tpm_lnc_1562_6473+1)
dat <- as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
group_list_group<-factor(coldata_all$GDMgroup2023)
group_list_batch<-factor(coldata_all$batch)
group_list_ethnicity<-as.factor(coldata_all$ethnicity)
group_list_trimester<-as.factor(coldata_all$trimester)
col_group<-c("GDM"="#E18E6D","control"="#62B197")
col_batch<-c("149"="#2f5688","200"="#CC0000")
col_ethnicity<-c("Han"="#E18E6D","minority"="#62B197")
col_trimester<-c("1"="#AED185","2"="#F89FA8","3"="#AECDD7")
#shape_group<-c("GDM"="#E18E6D","control"="#62B197")
#shape_batch<-c("149"="#2f5688","200"="#CC0000")

test<-data.frame(sample=rownames(dat.pca$ind$coord),pc1=dat.pca$ind$coord[,1],pc2=dat.pca$ind$coord[,2],
                 Batch=as.factor(coldata_all$batch),Group=as.factor(coldata_all$disease),Trimester=as.factor(coldata_all$trimester))
#test_mean<-aggregate(cbind(pc1,pc2)~Batch,data=test,FUN=mean)
#test<-merge(test,setNames(test_mean,c("Batch","pc1_mean","pc2_mean")))

p_batch<-ggplot(test,aes(x=pc1,y=pc2,color=Batch,shape=Batch))+
  geom_point(size=2,alpha=0.8)+
  stat_ellipse(aes(color=Batch),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_batch)+
  scale_shape_manual(values = c("149"=16,"200"=18))+
  #scale_color_manual(values = c(NA,NA))+
  xlab("P1 (26.64%)")+
  ylab("P2 (2.43%)")+
  ggtitle("lncRNA")+
  #xlim(c(-120,200))+
  #ylim(c(-10,10))+
  labs(fill = "Batch")+
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

p_trimester<-ggplot(test,aes(x=pc1,y=pc2,color=Trimester,shape=Trimester))+
  geom_point(size=2,alpha=1)+
  stat_ellipse(aes(color=Trimester),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_trimester)+
  scale_shape_manual(values = c("1"=16,"2"=15,"3"=18))+
  #scale_color_manual(values = c(NA,NA,NA))+
  xlab("P1 (26.64%)")+
  ylab("P2 (2.43%)")+
  ggtitle("lncRNA")+
  #xlim(c(-120,200))+
  #ylim(c(-10,10))+
  labs(fill = "Trimester")+
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

p_group<-ggplot(test,aes(x=pc1,y=pc2,color=Group,shape=Group))+
  geom_point(size=2,alpha=1)+
  stat_ellipse(aes(color=Group),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_group)+
  scale_shape_manual(values = c("GDM"=16,"control"=18))+
  xlab("P1 (26.64%)")+
  ylab("P2 (2.43%)")+
  ggtitle("lncRNA")+
  #xlim(c(-120,200))+
  #ylim(c(-10,10))+
  labs(fill = "Group")+
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

#cairo_pdf(filename = "cfRNA_m_mi_lncRNA_PCA质控.pdf",width = 12,height = 5,family = "Arial")

p_lncRNA<-p_group+p_batch+p_trimester

#dev.off()


##### pca-cfRNA miRNA #####
coldata_all<-read.csv("__PATH_003__",row.names = 1)
#fpkm_1562<-read.csv("../cfRNA_fpkm_1562.csv",row.names = 1) %>%
#  dplyr::select(coldata_all$sample)
fpkm_1562<-read.csv("__PATH_007__",row.names = 1,header = T) %>%
  dplyr::select(coldata_all$sample)


library(ggord)
library(FactoMineR)

exp <- log2(fpkm_1562+1)
dat <- as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
group_list_group<-factor(coldata_all$GDMgroup2023)
group_list_batch<-factor(coldata_all$batch)
group_list_ethnicity<-as.factor(coldata_all$ethnicity)
group_list_trimester<-as.factor(coldata_all$trimester)
col_group<-c("GDM"="#E18E6D","control"="#62B197")
col_batch<-c("149"="#2f5688","200"="#CC0000")
col_ethnicity<-c("Han"="#E18E6D","minority"="#62B197")
col_trimester<-c("1"="#AED185","2"="#F89FA8","3"="#AECDD7")

test<-data.frame(sample=rownames(dat.pca$ind$coord),pc1=dat.pca$ind$coord[,1],pc2=dat.pca$ind$coord[,2],
                 Batch=as.factor(coldata_all$batch),Group=as.factor(coldata_all$disease),Trimester=as.factor(coldata_all$trimester))
#test_mean<-aggregate(cbind(pc1,pc2)~Batch,data=test,FUN=mean)
#test<-merge(test,setNames(test_mean,c("Batch","pc1_mean","pc2_mean")))

p_batch<-ggplot(test,aes(x=pc1,y=pc2,color=Batch,shape=Batch))+
  geom_point(size=2,alpha=0.8)+
  stat_ellipse(aes(color=Batch),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_batch)+
  scale_shape_manual(values = c("149"=16,"200"=18))+
  #scale_color_manual(values = c(NA,NA))+
  xlab("P1 (23.06%)")+
  ylab("P2 (9.08%)")+
  ggtitle("miRNA")+
  labs(fill = "Batch")+
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

p_trimester<-ggplot(test,aes(x=pc1,y=pc2,color=Trimester,shape=Trimester))+
  geom_point(size=2,alpha=1)+
  stat_ellipse(aes(color=Trimester),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_trimester)+
  scale_shape_manual(values = c("1"=16,"2"=15,"3"=18))+
  #scale_color_manual(values = c(NA,NA,NA))+
  xlab("P1 (23.06%)")+
  ylab("P2 (9.08%)")+
  ggtitle("miRNA")+
  labs(fill = "Trimester")+
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

p_group<-ggplot(test,aes(x=pc1,y=pc2,color=Group,shape=Group))+
  geom_point(size=2,alpha=1)+
  stat_ellipse(aes(color=Group),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_group)+
  scale_shape_manual(values = c("GDM"=16,"control"=18))+
  xlab("P1 (23.06%)")+
  ylab("P2 (9.08%)")+
  ggtitle("miRNA")+
  labs(fill = "Group")+
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

p_miRNA<-p_group+p_batch+p_trimester



cairo_pdf(filename = "cfRNA_m_mi_lncRNA_PCA质控.pdf",height = 15,width = 12,family = "Arial")

p_mRNA/p_lncRNA/p_miRNA

dev.off()

##### motif #####
motif_mRNA<-read.csv("__PATH_008__",row.names = 1)
motif_lncRNA<-read.csv("__PATH_009__",row.names = 1)
motif_position<-data.frame(motif=names(motif_lncRNA)) %>%
  separate(motif, into = c("NA","p1","p2","p3","p4"), sep = "",remove = F) %>%
  dplyr::select(-"NA") %>%
  mutate(mRNA=colSums(motif_mRNA),lncRNA=colSums(motif_lncRNA))
#motif_mRNA_sum<-data.frame(mRNA=colSums(motif_mRNA),lncRNA=colSums(motif_lncRNA))

plot_motif<-list(
  lncRNA=as.matrix(data.frame(p1=c("A"=sum(motif_position$lncRNA[which(motif_position$p1=="A")]),
                                   "C"=sum(motif_position$lncRNA[which(motif_position$p1=="C")]),
                                   "G"=sum(motif_position$lncRNA[which(motif_position$p1=="G")]),
                                   "T"=sum(motif_position$lncRNA[which(motif_position$p1=="T")])
  ),
  p2=c("A"=sum(motif_position$lncRNA[which(motif_position$p2=="A")]),
       "C"=sum(motif_position$lncRNA[which(motif_position$p2=="C")]),
       "G"=sum(motif_position$lncRNA[which(motif_position$p2=="G")]),
       "T"=sum(motif_position$lncRNA[which(motif_position$p2=="T")])
  ),
  p3=c("A"=sum(motif_position$lncRNA[which(motif_position$p3=="A")]),
       "C"=sum(motif_position$lncRNA[which(motif_position$p3=="C")]),
       "G"=sum(motif_position$lncRNA[which(motif_position$p3=="G")]),
       "T"=sum(motif_position$lncRNA[which(motif_position$p3=="T")])
  ),
  p4=c("A"=sum(motif_position$lncRNA[which(motif_position$p4=="A")]),
       "C"=sum(motif_position$lncRNA[which(motif_position$p4=="C")]),
       "G"=sum(motif_position$lncRNA[which(motif_position$p4=="G")]),
       "T"=sum(motif_position$lncRNA[which(motif_position$p4=="T")])
  ))),
  mRNA=as.matrix(data.frame(p1=c("A"=sum(motif_position$mRNA[which(motif_position$p1=="A")]),
                                 "C"=sum(motif_position$mRNA[which(motif_position$p1=="C")]),
                                 "G"=sum(motif_position$mRNA[which(motif_position$p1=="G")]),
                                 "T"=sum(motif_position$mRNA[which(motif_position$p1=="T")])
  ),
  p2=c("A"=sum(motif_position$mRNA[which(motif_position$p2=="A")]),
       "C"=sum(motif_position$mRNA[which(motif_position$p2=="C")]),
       "G"=sum(motif_position$mRNA[which(motif_position$p2=="G")]),
       "T"=sum(motif_position$mRNA[which(motif_position$p2=="T")])
  ),
  p3=c("A"=sum(motif_position$mRNA[which(motif_position$p3=="A")]),
       "C"=sum(motif_position$mRNA[which(motif_position$p3=="C")]),
       "G"=sum(motif_position$mRNA[which(motif_position$p3=="G")]),
       "T"=sum(motif_position$mRNA[which(motif_position$p3=="T")])
  ),
  p4=c("A"=sum(motif_position$mRNA[which(motif_position$p4=="A")]),
       "C"=sum(motif_position$mRNA[which(motif_position$p4=="C")]),
       "G"=sum(motif_position$mRNA[which(motif_position$p4=="G")]),
       "T"=sum(motif_position$mRNA[which(motif_position$p4=="T")])
  )))
)

p_motif_lncRNA<-ggseqlogo(plot_motif$lncRNA,method="probability")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        #panel.border = element_rect(fill=NA),
        axis.text = element_text(color="black",face="bold"),
        axis.title = element_text(color = "black",face="bold"),
        plot.title = element_text(color="black",face="bold",size=15,hjust = 0.5),
        legend.text = element_text(face="bold"))

p_motif_mRNA<-ggseqlogo(plot_motif$mRNA,method="probability")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        #panel.border = element_rect(fill=NA),
        axis.text = element_text(color="black",face="bold"),
        axis.title = element_text(color = "black",face="bold"),
        plot.title = element_text(color="black",face="bold",size=15,hjust = 0.5),
        legend.text = element_text(face="bold"))

#ggsave(filename = "cfRNA_motif_mRNA.pdf",plot = p_motif_mRNA,width = 6,height = 6)
cairo_pdf(filename = "cfRNA_m_lncRNA_motif.pdf",width = 8,height = 4,family = "Arial")
p_motif_mRNA+p_motif_lncRNA
dev.off()

##### gene num detected #####
gene_num<-read_excel("__PATH_010__") %>%
  dplyr::rename(sample="...1") %>%
  pivot_longer(cols = -sample,names_to = "class",values_to = "count")
rna_level=factor(c("mRNA","lncRNA","miRNA","tRNA","piRNA"))
plot_violin<-summary(gene_num[,c(2,3)])

ggplot(data = gene_num,aes(x=class,y=count,fill=class))+
  geom_boxplot()+
  scale_x_discrete(limits=rna_level)+
  scale_fill_manual(values = c("mRNA"="#E18E6D","lncRNA"="#62B197","miRNA"="#F2BB6B","tRNA"="#8DA0CB","piRNA"="#E78AC3"))+
  xlab("RNA type")+
  ylab("Number")+
  labs(fill="RNA type")+
  ggtitle("Gene")
theme(panel.background = element_blank(),
      axis.line = element_line(),
      #panel.border = element_rect(fill=NA),
      axis.text = element_text(color="black",face="bold"),
      axis.title = element_text(color = "black",face="bold"),
      plot.title = element_text(color="black",face="bold",size=15,hjust = 0.5),
      legend.text = element_text(face="bold"))



