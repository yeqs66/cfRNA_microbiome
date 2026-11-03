library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(extrafont)

coldata_rna<-read.table("__PATH_001__",header = T) %>%
  mutate(ID=as.character(ID))
res_deconformer<-read.table("细胞溯源/res_deconv_349.txt") %>%
  mutate(sample=rownames(.)) %>%
  left_join(coldata_rna[]) %>%
  pivot_longer(cols = !c(ID,sample,trimester,GDMgroup2023,batch),names_to = "cell_type",values_to = "score")

for (i in unique(res_deconformer$trimester)) {
  res_tmp<-subset(res_deconformer,res_deconformer$trimester==i)
  ggplot(data = res_tmp,aes(x=cell_type,y=score,fill=GDMgroup2023))+
    geom_boxplot(outlier.shape = NA)+
    #coord_cartesian(ylim = c(0, 0.1))+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    #panel.background = element_rect(fill = "white"),
    #panel.grid.major = element_line(color = "gray"))+
    theme_bw()+
    #geom_signif(comparisons = list(c("GDM","control")),test = "wilcox.test")
    stat_compare_means(aes(group = GDMgroup2023),method = "wilcox.test",label.y = 0.9,label = "p.signif")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text = element_text(face="bold"),
          #legend.position = "bottom",
          panel.border = element_blank())
}

selected_cell<-c("cell_of_skeletal_muscle","pancreatic_ductal_cell","pancreatic_acinar_cell","pancreatic_stellate_cell",
                 "pancreatic_alpha_cell__pancreatic_beta_cell","hepatocyte","early_VCT","late_VCT",
                 "early_SCT","late_SCT","early_EVT","late_EVT")


for (i in unique(res_deconformer$trimester)) {
  for (j in unique(res_deconformer$GDMgroup2023)) {
 res_tmp<-res_deconformer %>%
  dplyr::filter(GDMgroup2023==j,trimester==i) %>%
  group_by(cell_type) %>%
  summarise(mean_score = mean(score, na.rm = TRUE)) %>%
  mutate(type=ifelse(mean_score>0.01 | cell_type %in% selected_cell,cell_type,"other_tmp")) %>%
  arrange(desc(mean_score))
res_tmp<-res_tmp %>%
  filter(type == "other_tmp") %>%
  summarise(mean_score = sum(mean_score, na.rm = TRUE), type = "other") %>%
  # 将新生成的行添加到原始数据框的末尾
  bind_rows(res_tmp, .) %>%
  dplyr::filter(!type=="other_tmp") %>%
  mutate(cell_type=type) %>%
  mutate(cell_type=gsub("_"," ",cell_type),cell_type=gsub("  "," ",cell_type)) %>%
  dplyr::rename("Cell type"=cell_type) %>%
  mutate(`Cell type`=factor(`Cell type`,levels=`Cell type`),score_percent=paste0(round(mean_score,digits = 4)*100,"%"),label="")

#cairo_pdf(paste0("__PATH_002__","_",i,"_",j,".pdf"),width = 10,height = 6,family = "Arial")
p<-ggdonutchart(res_tmp, "mean_score", 
              label = "label",
              fill = "Cell type",
              lab.pos = c("out"),
              lab.adjust = 0,
              #palette =colorRampPalette(pal_frontiers("default",alpha = 0.5)(10))(length(res_tmp$`Cell type`)),
              lab.font = 2,
             color="grey70",font.family = "Arial",size=0.05)+
  scale_fill_manual(values = alpha(colorRampPalette(pal_locuszoom("default")(7))(length(res_tmp$`Cell type`)),alpha=0.6),labels=paste(res_tmp$`Cell type`,res_tmp$score_percent,sep = " "))+
  ggtitle(paste(i,j,sep=" "))+
  guides(fill=guide_legend(ncol=2))+
  theme(legend.position = "right",
        plot.title = element_text(vjust = 1.5,hjust = 0.5),
        legend.key = element_blank())

#dev.off()
    ggsave(filename = paste0("__PATH_003__",i,"_",j,"_donut.pdf"),p,width = 10,height = 6)
  }
}


## overall GDM vs control
res_tmp<-res_deconformer %>%
  dplyr::filter(GDMgroup2023=="control") %>%
  group_by(cell_type) %>%
  summarise(mean_score = mean(score, na.rm = TRUE)) %>%
  mutate(type=ifelse(mean_score>0.01 | cell_type %in% selected_cell,cell_type,"other_tmp")) %>%
  arrange(desc(mean_score))
res_tmp<-res_tmp %>%
  filter(type == "other_tmp") %>%
  summarise(mean_score = sum(mean_score, na.rm = TRUE), type = "other") %>%
  # 将新生成的行添加到原始数据框的末尾
  bind_rows(res_tmp, .) %>%
  dplyr::filter(!type=="other_tmp") %>%
  mutate(cell_type=type) %>%
  mutate(cell_type=gsub("_"," ",cell_type),cell_type=gsub("  "," ",cell_type)) %>%
  dplyr::rename("Cell type"=cell_type) %>%
  mutate(`Cell type`=factor(`Cell type`,levels=`Cell type`),score_percent=paste0(round(mean_score,digits = 4)*100,"%"),label="")

#cairo_pdf(paste0("__PATH_002__","_",i,"_",j,".pdf"),width = 10,height = 6,family = "Arial")
p<-ggdonutchart(res_tmp, "mean_score", 
              label = "label",
              fill = "Cell type",
              lab.pos = c("out"),
              lab.adjust = 0,
              #palette =colorRampPalette(pal_frontiers("default",alpha = 0.5)(10))(length(res_tmp$`Cell type`)),
              lab.font = 2,
             color="grey70",font.family = "Arial",size=0.05)+
  scale_fill_manual(values = alpha(colorRampPalette(pal_locuszoom("default")(7))(length(res_tmp$`Cell type`)),alpha=0.6),labels=paste(res_tmp$`Cell type`,res_tmp$score_percent,sep = " "))+
  ggtitle("Control")+
  guides(fill=guide_legend(ncol=2))+
  theme(legend.position = "right",
        plot.title = element_text(vjust = 1.5,hjust = 0.5),
        legend.key = element_blank())

#dev.off()
    ggsave(filename = "__PATH_004__",p,width = 10,height = 6)
