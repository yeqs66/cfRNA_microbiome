library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)

#coldata_rna<-read.table("../cfRNA_coldata_1562_batch.txt",header = T) %>%
#  mutate(ID=as.character(ID))
coldata_rna<-read.csv("__PATH_001__",,row.names = 1) %>%
  mutate(ID=as.character(ID)) %>%
  dplyr::select(sample,ID,trimester,disease,batch)
res_deconformer<-read.table("细胞溯源/res_deconv_349.txt") %>%
  mutate(sample=rownames(.)) %>%
  left_join(coldata_rna[]) %>%
  pivot_longer(cols = !c(ID,sample,trimester,disease,batch),names_to = "cell_type",values_to = "score")

res_deconformer_stat_wilcox<-data.frame()
res_deconformer_stat_ttest<-data.frame()
for (i in unique(res_deconformer$cell_type)) {
  tmp_1<-res_deconformer[which(res_deconformer$cell_type==i & res_deconformer$trimester==1),]
  tmp_2<-res_deconformer[which(res_deconformer$cell_type==i & res_deconformer$trimester==2),]
  tmp_3<-res_deconformer[which(res_deconformer$cell_type==i & res_deconformer$trimester==3),]
  tmp_all<-res_deconformer[which(res_deconformer$cell_type==i),]
  
  tmp1<-wilcox.test(tmp_1$score[which(tmp_1$disease=="GDM")],tmp_1$score[which(tmp_1$disease=="control")])
  tmp2<-wilcox.test(tmp_2$score[which(tmp_2$disease=="GDM")],tmp_2$score[which(tmp_2$disease=="control")])
  tmp3<-wilcox.test(tmp_3$score[which(tmp_3$disease=="GDM")],tmp_3$score[which(tmp_3$disease=="control")])
  tmpall<-wilcox.test(tmp_all$score[which(tmp_all$disease=="GDM")],tmp_all$score[which(tmp_all$disease=="control")])
  res_deconformer_stat_wilcox<-rbind(res_deconformer_stat_wilcox,c(i,tmp1$p.value,tmp2$p.value,tmp3$p.value,tmpall$p.value))
  
  tmp1<-t.test(tmp_1$score[which(tmp_1$disease=="GDM")],tmp_1$score[which(tmp_1$disease=="control")])
  tmp2<-t.test(tmp_2$score[which(tmp_2$disease=="GDM")],tmp_2$score[which(tmp_2$disease=="control")])
  tmp3<-t.test(tmp_3$score[which(tmp_3$disease=="GDM")],tmp_3$score[which(tmp_3$disease=="control")])
  tmpall<-t.test(tmp_all$score[which(tmp_all$disease=="GDM")],tmp_all$score[which(tmp_all$disease=="control")])
  res_deconformer_stat_ttest<-rbind(res_deconformer_stat_ttest,c(i,tmp1$p.value,tmp2$p.value,tmp3$p.value,tmpall$p.value))
}
rm(tmp_1,tmp_2,tmp_3,tmp_all,tmp1,tmp2,tmp3,tmpall)
names(res_deconformer_stat_wilcox)<-c("cell_type","1st Trimester","2nd Trimester","3rd Trimester","All_period")
names(res_deconformer_stat_ttest)<-c("cell_type","1st Trimester","2nd Trimester","3rd Trimester","All_period")

##折线图
data_plot<-res_deconformer %>%
  dplyr::filter(cell_type=="pancreatic_alpha_cell__pancreatic_beta_cell") %>%
  #mutate(trimester=factor(trimester,levels=c("1st","2nd","3rd")))
  mutate(color=disease,color=gsub("GDM","#E41A1C",color),color=gsub("control","#377EB8",color))
ggplot(data_plot,aes(x=trimester,y=score,fill=disease,group=ID)) +
  geom_line(aes(color=color,alpha=0.4))+
  geom_point(size=2,aes(color=color,alpha=0.5)) +
  geom_smooth(aes(group = disease, color = color),
              method = "lm", se = TRUE,
              linewidth = 1.2, alpha = 0.15)+
  scale_color_identity()+
  ylim(0,0.01)+
  theme(legend.position = "none")


##柱状图
res_deconformer_noblood<-res_deconformer %>%
  dplyr::filter(!cell_type=="erythrocyte__erythroid_progenitor",
                trimester==3)
ggplot(res_deconformer_noblood,
       aes(x = disease, y = score, color = disease)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(comparisons = list(c("GDM", "control")),
              test = "wilcox.test",
              map_signif_level = F,   # 自动 */**/***
              y_position  = 0.05,
              textsize = 3,
              color="grey50") +       # 括号高度
  #ylim(c(0,0.1))+
  coord_cartesian(ylim = c(0, 0.1),clip = "off") +
  facet_wrap(~cell_type) +
  theme(plot.margin = margin(t = 20, unit = "pt"))+
  theme_bw()


for (i in 1:4) {
  period<-c("1st Trimester","2nd Trimester","3rd Trimester","All_period")[i]
  tmp<-res_deconformer_stat_wilcox$cell_type[which(res_deconformer_stat_wilcox[[period]]<0.05)]

  for (j in tmp) {
    
    if (period=="All_period") {
      data_plot<-res_deconformer[which(res_deconformer$cell_type==j),]
    } else {
      data_plot<-res_deconformer[which(res_deconformer$cell_type==j & res_deconformer$trimester==i),]
    }
    
    p<-ggplot(data_plot,
           aes(x = disease, y = score, color = disease)) +
      geom_boxplot(outlier.shape = NA) +
      scale_color_manual(values=c("GDM"="#E41A1C","control"="#377EB8"))+
      geom_signif(comparisons = list(c("GDM", "control")),
                  test = "t.test",
                  map_signif_level = F,   # 自动 */**/***
                  #y_position  = 0.05,
                  textsize = 3,
                  color="grey50") +       # 括号高度
      #ylim(c(0,0.1))+
      ggtitle(period,subtitle = j)+
      #coord_cartesian(ylim = c(0, 0.1),clip = "off") +
      #facet_wrap(~cell_type) +
      theme(plot.margin = margin(t = 20, unit = "pt"))+
      theme_bw()
    ggsave(filename = paste0("__PATH_002__",period,"_",j,"_wilcox.pdf"),p,width = 4.5,height = 4)
  }
  
}


##### 去除血细胞 #####
##柱状图
blood_cells <- c(
  "innate_lymphoid_cell",
  "neutrophil",
  "basophil",
  "b_cell",
  "monocyte",
  "plasmablast",
  "plasma_cell",
  "erythrocyte__erythroid_progenitor",
  "mature_conventional_dendritic_cell",
  "myeloid_progenitor",
  "mast_cell",
  "platelet",
  "macrophage",
  "hematopoietic_stem_cell",
  "t_cell",
  "thymocyte",
  "nk_cell"
)

res_deconformer_noblood<-res_deconformer %>%
  mutate(is_blood= cell_type %in% blood_cells) %>%
  #dplyr::filter(trimester=="3") %>%
  group_by(sample) %>%
  mutate(sum_noblood=sum(score[!is_blood])) %>%
  ungroup() %>%
  filter(!is_blood) %>%                    # 如果你不想保留血细胞，直接去掉
  mutate(score_noblood = score / sum_noblood,
         cell_type=gsub("_"," ",cell_type),cell_type=gsub(" "," ",cell_type),cell_type=gsub("  "," ",cell_type))

ggplot(res_deconformer_noblood,
       aes(x = disease, y = score_noblood, color = disease)) +
  geom_boxplot(outlier.shape = NA) +
  geom_signif(comparisons = list(c("GDM", "control")),
              test = "wilcox.test",
              map_signif_level = F,   # 自动 */**/***
              y_position  = 0.1,
              textsize = 3,
              color="grey50") +       # 括号高度
  #ylim(c(0,0.1))+
  coord_cartesian(ylim = c(0, 0.4),clip = "off") +
  facet_wrap(~cell_type) +
  theme(plot.margin = margin(t = 20, unit = "pt"))+
  theme_bw()

res_deconformer_noblood_stat_wilcox<-data.frame()
for (i in unique(res_deconformer_noblood$cell_type)) {
  tmp_1<-res_deconformer_noblood[which(res_deconformer_noblood$cell_type==i & res_deconformer_noblood$trimester==1),]
  tmp_2<-res_deconformer_noblood[which(res_deconformer_noblood$cell_type==i & res_deconformer_noblood$trimester==2),]
  tmp_3<-res_deconformer_noblood[which(res_deconformer_noblood$cell_type==i & res_deconformer_noblood$trimester==3),]
  tmp_all<-res_deconformer_noblood[which(res_deconformer_noblood$cell_type==i),]
  
  tmp1<-wilcox.test(tmp_1$score_noblood[which(tmp_1$disease=="GDM")],tmp_1$score_noblood[which(tmp_1$disease=="control")])
  tmp2<-wilcox.test(tmp_2$score_noblood[which(tmp_2$disease=="GDM")],tmp_2$score_noblood[which(tmp_2$disease=="control")])
  tmp3<-wilcox.test(tmp_3$score_noblood[which(tmp_3$disease=="GDM")],tmp_3$score_noblood[which(tmp_3$disease=="control")])
  tmpall<-wilcox.test(tmp_all$score_noblood[which(tmp_all$disease=="GDM")],tmp_all$score_noblood[which(tmp_all$disease=="control")])
  res_deconformer_noblood_stat_wilcox<-rbind(res_deconformer_noblood_stat_wilcox,c(i,tmp1$p.value,tmp2$p.value,tmp3$p.value,tmpall$p.value))
}
rm(tmp_1,tmp_2,tmp_3,tmp_all,tmp1,tmp2,tmp3,tmpall)
names(res_deconformer_noblood_stat_wilcox)<-c("cell_type","1st Trimester","2nd Trimester","3rd Trimester","All_period")


preg_cell<-c("early_VCT","late_VCT","early_SCT","late_SCT","early_EVT","late_EVT")

for (i in 1:4) {
  for (cell in preg_cell){
    period<-c("1st Trimester","2nd Trimester","3rd Trimester","All_period")[i]
    tmp<-res_deconformer_noblood_stat_wilcox %>%
      dplyr::filter(cell_type==cell)
    
    for (j in tmp) {
      
      if (period=="All_period") {
        data_plot<-res_deconformer_noblood[which(res_deconformer_noblood$cell_type==j),]
      } else {
        data_plot<-res_deconformer_noblood[which(res_deconformer_noblood$cell_type==j & res_deconformer_noblood$trimester==i),]
      }
      
      p<-ggplot(data_plot,
                aes(x = disease, y = score_noblood, color = disease)) +
        geom_boxplot(outlier.shape = NA) +
        scale_color_manual(values=c("GDM"="#E41A1C","control"="#377EB8"))+
        geom_signif(comparisons = list(c("GDM", "control")),
                    test = "wilcox.test",
                    map_signif_level = F,   # 自动 */**/***
                    y_position  = 0.05,
                    textsize = 3,
                    color="grey50") +       # 括号高度
        #ylim(c(0,0.1))+
        ggtitle(period,subtitle = j)+
        #coord_cartesian(ylim = c(0, 0.1),clip = "off") +
        #facet_wrap(~cell_type) +
        #theme(plot.margin = margin(t = 20, unit = "pt"))+
        theme_bw()
      ggsave(filename = paste0("__PATH_003__",period,"_",j,"_wilcox.pdf"),p,width = 4.5,height = 4)
    }

  }
  
}


##### special cell type #####
preg_cell<-data.frame(ori_cell=c("early VCT","late VCT","early SCT","late SCT","early EVT","late EVT"),
                      new_cell=c("The proportion of Early Villous cytotrophoblast (VCT)","The proportion of Late Villous cytotrophoblast (VCT)",
                      "The proportion of Early Syncytiotrophoblast (SCT)","The proportion of Late Syncytiotrophoblast (SCT)",
                      "The proportion of Early Extra-villous trophoblast (EVT)","The proportion of Late Extra-villous trophoblast (EVT)"))
tri<-NA

for (cell in preg_cell$ori_cell) {
  if (is.na(tri)) {
      data_plot<-res_deconformer_noblood %>%
        left_join(preg_cell,by=c("cell_type"="ori_cell")) %>%
        dplyr::filter(cell_type==cell)
  } else {
    data_plot<-res_deconformer_noblood %>%
      dplyr::filter(trimester==tri) %>%
      left_join(preg_cell,by=c("cell_type"="ori_cell")) %>%
        dplyr::filter(cell_type==cell)
  }

  cell<-preg_cell$new_cell[which(preg_cell$ori_cell==cell)]

  p<-ggplot(data_plot,
        aes(x = disease, y = score_noblood, color = disease)) +
        geom_boxplot(outlier.shape = NA) +
        scale_color_manual(values=c("GDM"="#E41A1C","control"="#377EB8"))+
        geom_signif(comparisons = list(c("GDM", "control")),
                    test = "wilcox.test",
                    map_signif_level = F,   # 自动 */**/***
                    y_position  = 0.04,
                    textsize = 3,
                    margin_top=0.005,
                    tip_length = 0.001,
                    color="grey50") +       # 括号高度
        #ylim(c(0,0.1))+
        ggtitle("Overall trimester",subtitle=paste0(cell))+
        #ggtitle("3rd Trimester",subtitle = "intrahepatic cholangiocyte")+
        coord_cartesian(ylim = c(0, 0.06),clip = "off") +
        #facet_wrap(~cell_type) +
        #theme(plot.margin = margin(t = 20, unit = "pt"))+
        theme_bw()
ggsave(plot=p,filename = paste0("__PATH_004__","_",cell,"_wilcox.pdf"),width = 4.5,height = 4)

}


## late intrahepatic_cholangiocyte
p<-ggplot(data_plot,
       aes(x = disease, y = score_noblood, color = disease)) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values=c("GDM"="#E41A1C","control"="#377EB8"))+
  geom_signif(comparisons = list(c("GDM", "control")),
              test = "wilcox.test",
              map_signif_level = F,   # 自动 */**/***
              y_position  = 0.004,
              textsize = 3,
              margin_top=0.005,
              tip_length = 0.001,
              color="grey50") +       # 括号高度
  #ylim(c(0,0.1))+
  ggtitle("3rd Trimester",subtitle = "intrahepatic cholangiocyte")+
  coord_cartesian(ylim = c(0, 0.006),clip = "off") +
  #facet_wrap(~cell_type) +
  #theme(plot.margin = margin(t = 20, unit = "pt"))+
  theme_bw()
ggsave(filename = "__PATH_005__",p,width = 4.5,height = 4)

