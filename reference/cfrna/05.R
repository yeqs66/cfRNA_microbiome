library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

target_genes <- c("ALDOB", "CPT1A", "FERMT3")

deg_all <- read.csv("limma_349_batch_only.csv") %>% 
  dplyr::filter(period=="all") %>%
  mutate(logp=-1*log10(adj.P.Val),sign=rep("up",nrow(.)),label="") %>%
  arrange(desc(logp)) %>%
  mutate(markedlabel = case_when(adj.P.Val < 0.1 & abs(logFC) > 0.15 ~ "marked",
                                 TRUE ~"unmarked"),
         highlight_fill = ifelse(SYMBOL %in% target_genes, "Target", "Normal"),)
deg_all$sign[which(deg_all$logFC<0)]<-"down"
deg_all$group<-"Non-significant"
deg_all$group[which(deg_all$logFC>0.15 & deg_all$logp>1)]<-"Up-regulated"
deg_all$group[which(deg_all$logFC< -0.15 & deg_all$logp>1)]<-"Down-regulated"
upgene<-head(deg_all$SYMBOL[which(deg_all$group=="Up-regulated")],10)
downgene<-head(deg_all$SYMBOL[which(deg_all$group=="Down-regulated")],10)
deg_sign<-c(as.character(upgene),as.character(downgene))
deg_all$label[match(deg_sign,deg_all$SYMBOL)]<-deg_sign




ggscatter(deg_all,x="logFC",y="logp",
          color = "group",
          palette = c("#2f5688","#BBBBBB","#CC0000"),
          size = 1,
          label = deg_all$label,
          font.label = 8,
          repel = T,
          xlab = "log2FC",
          ylab = "-log10(p.adj)") +
  #theme_base() +
  geom_hline(yintercept = 1,linetype = "dashed") +
  geom_vline(xintercept = c(-0.15,0.15),linetype = "dashed")

p<-ggplot(data = deg_all, aes(x=logFC, y=logp,color=group))+
  geom_point(alpha=0.6,size=4.5)+
  geom_text_repel(data=subset(deg_all, deg_all$group=="Up-regulated"), 
                  aes(x = logFC, y = logp, label = SYMBOL),
                  seed = 233,                           # 设置随机数种子，用于确定标签位置
                  size = 3.5,                           # 设置标签的字体大小
                  color = 'black',                      # 设置标签的字体颜色
                  min.segment.length = 0,               # 始终为标签添加指引线段
                  force = 3,                            # 设置标签重叠时的排斥力
                  force_pull = 2,                       # 设置标签与数据点之间的吸引力
                  box.padding = 0.4,                    # 设置标签周围的填充量
                  max.overlaps = Inf,                   # 设置排斥重叠过多标签的阈值为无穷大，保持始终显示所有标签
                  segment.linetype = 3,                 # 设置线段类型为虚线
                  segment.color = 'black',              # 设置线段颜色
                  segment.alpha = 0.5,                  # 设置线段不透明度
                  fontface = "italic",
                  nudge_x = 2 - subset(deg_all, deg_all$group=="Up-regulated")$logFC,          # 调整标签x轴起始位置
                  direction = "y",                      # 按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 1)+
  geom_text_repel(data=subset(deg_all, deg_all$group=="Down-regulated"), 
                  aes(x = logFC, y = logp, label = SYMBOL),
                  seed = 233,                           # 设置随机数种子，用于确定标签位置
                  size = 3.5,                           # 设置标签的字体大小
                  color = 'black',                      # 设置标签的字体颜色
                  min.segment.length = 0,               # 始终为标签添加指引线段
                  force = 4,                            # 设置标签重叠时的排斥力
                  force_pull = 2,                       # 设置标签与数据点之间的吸引力
                  box.padding = 0.4,                    # 设置标签周围的填充量
                  max.overlaps = Inf,                   # 设置排斥重叠过多标签的阈值为无穷大，保持始终显示所有标签
                  segment.linetype = 3,                 # 设置线段类型为虚线
                  segment.color = 'black',              # 设置线段颜色
                  segment.alpha = 0.5,                  # 设置线段不透明度
                  fontface = "italic",
                  nudge_x = -2 - subset(deg_all, deg_all$group=="Down-regulated")$logFC,          # 调整标签x轴起始位置
                  direction = "y",                      # 按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 1)+
  scale_color_manual(values=c("Down-regulated"='#62B197',"Up-regulated"='#E18E6D',"Non-significant"="grey80")) +
  scale_fill_manual((values=c("Target"='#FFEB99', "Normal"='white')))+
  #guides(size = guide_legend("-log10(FDR)"))+
  xlab("Log2 fold change") + 
  ylab("-Log10 FDR") + 
  #scale_x_discrete(labels = c("Control 2 vs Control 1", "GDM 2 vs 1","Control 3 vs 2", "GDM 3 vs 2","GDM vs Control"))+
  theme_classic() +
  ggtitle("GDM VS Control") +
  geom_hline(yintercept = 1,linetype = "dashed") +
  geom_vline(xintercept = c(-0.15,0.15),linetype = "dashed")+
  theme(axis.text=element_text(size=14,face="bold", color = "black"), 
        axis.title=element_text(size=14,face="bold"), 
        legend.text = element_text(size=14,face="bold"), 
        legend.title = element_blank(), 
        plot.title = element_text(size=14, face="bold", hjust = 0.5), 
        axis.text.x = element_text(size=14,face="bold"))

ggsave(filename = "volcano_deg_mRNA_all.pdf", plot = p, width = 8, height = 10, dpi = 600)
