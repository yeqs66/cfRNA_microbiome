library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(readxl)
library(patchwork)
library(scales)

## load data
module7<-read.csv("wgcna/cfRNA_microbiome/m_mi_lnc/ko/wgcga_signed_signed_hybrid_spearman_ori/TissueSpecificGeneEnrichment_module_red.tsv",sep = "\t") %>%
  dplyr::rename(log10padj="X.Log10PValue") %>%
  mutate(log10padj=as.numeric(log10padj),
         Tissue=factor(Tissue,levels=Tissue))

module13<-read.csv("wgcna/cfRNA_microbiome/m_mi_lnc/ko/wgcga_signed_signed_hybrid_spearman_ori/TissueSpecificGeneEnrichment_module_tan.tsv",sep = "\t") %>%
  dplyr::rename(log10padj="X.Log10PValue") %>%
  mutate(log10padj=as.numeric(log10padj),
         Tissue=factor(Tissue,levels=Tissue))

nonzero_tissues <- as.character(module7$Tissue[module7$log10padj > 0])
nz_cols <- hue_pal()(length(nonzero_tissues))
names(nz_cols) <- nonzero_tissues

all_cols <- c(nz_cols, setNames(rep("grey80", sum(module7$log10padj == 0)),
                                as.character(module7$Tissue[module7$log10padj == 0])))
p1<-ggplot(module7, aes(x = Tissue, y = log10padj, fill = Tissue)) +
  geom_col(width = 0.8) +
  scale_fill_nejm(palette = "default") +
  #labs(x = "Tissue", y = "Value", title = "Tissue Distribution") +
  theme_bw(base_size = 12) +
  xlab("Tissue")+
  ylab("-Log10(P.Adjusted)")+
  ggtitle("TissueEnrich results of Module 7")+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )




p2<-ggplot(module13, aes(x = Tissue, y = log10padj, fill = Tissue)) +
  geom_col(width = 0.8) +
  scale_fill_npg(palette = "nrc") +
  #labs(x = "Tissue", y = "Value", title = "Tissue Distribution") +
  theme_bw(base_size = 12) +
  xlab("Tissue")+
  ylab("-Log10(P.Adjusted)")+
  ggtitle("TissueEnrich results of Module 13")+
  theme(
    axis.text.x = element_text(family = "Arial",angle = 60, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


cairo_pdf(filename = "wgcna/cfRNA_microbiome/m_mi_lnc/ko/wgcga_signed_signed_hybrid_pearson/TissueEnrich.pdf",width = 6,height = 8,family = "Arial")
p<-p1/p2
p

dev.off()
