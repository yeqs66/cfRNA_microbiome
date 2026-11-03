library(latex2exp)
library(grid)
library(openxlsx)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(ggplot2)
library(gridExtra)
library(fgsea)
library(org.Hs.eg.db)


load("__PATH_001__")
pathways_kegg <- pathways$transcriptome[grep("(KEGG)",names(pathways$transcriptome))]
tmp <- do.call(rbind, lapply(seq_along(pathways$transcriptome), function(i) {
  data.frame(name=rep(names(pathways$transcriptome[i]),length(pathways$transcriptome[[i]])), value=pathways$transcriptome[[i]], stringsAsFactors=FALSE)
})) %>%
  dplyr::filter(grepl("\\(KEGG\\)",name))
symbol_ensg<-read.table("__PATH_002__") %>%
  dplyr::rename_with(~c("ENSG","SYMBOL")) %>%
  mutate(ENSG=gsub("\\..*","",ENSG))

##### vs deg #####
deg_res_349_nocoef_addbatch<-read.csv("__PATH_003__")

gsea_deg_kegg<-data.frame()
for (i in unique(deg_res_349_nocoef_addbatch$period)) {
  temp<-deg_res_349_nocoef_addbatch %>%
    dplyr::filter(period==i) %>%
    arrange(desc(logFC))
  gene_fc<-temp$logFC
  names(gene_fc)<-temp$SYMBOL
  
  #symbol_ensg_flt<-symbol_ensg[match(names(gene_fc),symbol_ensg$V2),]
  #all(symbol_ensg_flt$V2==names(gene_fc))
  #names(gene_fc)<-symbol_ensg_flt$V1
  
  set.seed(3490)
  gsea_other<-fgsea(pathways_kegg, stats=gene_fc) %>%
    mutate(period=i)
  
  #gsea_other<-GSEA(geneList = gene_fc,TERM2GENE = pathways,pvalueCutoff = 1,pAdjustMethod = "BH")
  
  #gsea_other<-as.data.frame(gsea_other@result) %>%
  #  mutate(type="other",period=i)
  
  gsea_deg_kegg<-rbind(gsea_deg_kegg,gsea_other)
}

gsea_deg_kegg<-gsea_deg_kegg %>%
  arrange(padj)
#save(gsea_deg_kegg,file = "gsea_deg_")

tmp_gsea<-gsea_deg_kegg %>%
  arrange(padj) %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::mutate(log10padj= -log10(padj)) %>%
  mutate(pathway=gsub("\\(KEGG\\)","",pathway),
         period=gsub("1","1st",period),period=gsub("2","2nd",period),
         period=gsub("3","3rd",period),period=gsub("all","All Period",period))

tmp_gsea_other<-gsea_deg_kegg %>%
  mutate(pathway=gsub("\\(KEGG\\)","",pathway),
         period=gsub("1","1st",period),period=gsub("2","2nd",period),
         period=gsub("3","3rd",period),period=gsub("all","All Period",period)) %>%
  arrange(padj) %>%
  dplyr::filter(pathway %in% tmp_gsea$pathway,padj>0.05) %>%
  mutate(log10padj= -log10(padj))

tmp_gsea_plot<-rbind(tmp_gsea,tmp_gsea_other) %>%
  mutate(sig=if_else(padj<0.05,"Significant","NonSignificant")) %>%
  dplyr::filter(pathway %in% tmp_gsea$pathway[which(tmp_gsea$padj<0.05)])

#tmp_gsea<-read_excel("../../../../../文章/多组学/4.supplymentary_cfRNA_GSEA.xlsx") %>%
#  mutate(pathway=gsub("\\(KEGG\\)","",pathway))

pathway_level=as.factor(unique(tmp_gsea_plot$pathway))


p<-ggplot(tmp_gsea_plot,aes(x=period ,y=pathway,size=log10padj,fill=NES,shape=sig))+
  geom_point(aes(size = ifelse(sig == "NonSignificant", 0, log10padj),stroke = ifelse(sig == "NonSignificant", 0, 0.04)))+
  #geom_point(size=2)+
  scale_size_continuous(range = c(3,7))+
  scale_y_discrete(limits=rev(pathway_level))+
  scale_fill_gradient(low="#62B197",high="#E18E6D")+
  scale_shape_manual(values = c("Significant" = 21,"NonSignificant" = 24), # 形状24对应三角形。确保名称与aes内的字符串一致
                     name = "Significant")+
  ylab("Pathways")+
  xlab("Trimester")+
  labs(size=TeX('$-log_{10}(p.adj)$',bold = T))+
  theme(axis.text.x = element_text(color="black",face="bold",size=11),
        axis.text.y = element_text(color="black",size=11),
        axis.title.x = element_text(color="black",face="bold",size=14),
        axis.title.y = element_text(color="black",face="bold",size=14),
        axis.title = element_text(face="bold"),
        plot.title = element_text(color="black",face="bold"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(color="grey92",linetype=1,linewidth = 0.2),
        panel.grid.major.x=element_line(color="grey92",linetype=1,linewidth = 0.2),
        legend.key = element_rect(fill = "white",colour = NA),
        legend.text = element_text(color="black",size=14),
        legend.title = element_text(color="black",face="bold",size=14))


##ori plot
ggplot(tmp_gsea,aes(x=period ,y=pathway,size=log10padj,fill=NES))+
  geom_point(shape=21,stroke=0.04)+
  scale_size_continuous(range = c(3,7))+
  scale_y_discrete(limits=rev(pathway_level))+
  scale_fill_gradient(low="#62B197",high="#E18E6D")+
  ylab("Pathways")+
  xlab("Trimester")+
  labs(size=TeX('$-log_{10}(p.adj)$',bold = T))+
  theme(axis.text.x = element_text(color="black",face="bold",size=11),
        axis.text.y = element_text(color="black",size=11),
        axis.title.x = element_text(color="black",face="bold",size=14),
        axis.title.y = element_text(color="black",face="bold",size=14),
        axis.title = element_text(face="bold"),
        plot.title = element_text(color="black",face="bold"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(color="grey92",linetype=1,linewidth = 0.2),
        panel.grid.major.x=element_line(color="grey92",linetype=1,linewidth = 0.2),
        legend.key = element_rect(fill = "white",colour = NA),
        legend.text = element_text(color="black",size=14),
        legend.title = element_text(color="black",face="bold",size=14))



##### GO #####
gsea_deg_go<-data.frame()
for (i in unique(deg_res_349_nocoef_addbatch$period)) {
  temp<-deg_res_349_nocoef_addbatch %>%
    dplyr::filter(period==i) %>%
    arrange(desc(logFC))
  gene_fc<-temp$logFC
  names(gene_fc)<-temp$SYMBOL
  
  #symbol_ensg_flt<-symbol_ensg[match(names(gene_fc),symbol_ensg$V2),]
  #all(symbol_ensg_flt$V2==names(gene_fc))
  #names(gene_fc)<-symbol_ensg_flt$V1
  
  set.seed(3490)
  gsea_bp<-gseGO(geneList = gene_fc,ont = "BP",OrgDb = org.Hs.eg.db,keyType = "SYMBOL", 
                    pAdjustMethod = "BH",pvalueCutoff = 0.5,seed = 3490) %>%
    mutate(period=i,ont="BP")
  gsea_mf<-gseGO(geneList = gene_fc,ont = "MF",OrgDb = org.Hs.eg.db,keyType = "SYMBOL", 
                 pAdjustMethod = "BH",pvalueCutoff = 0.5,seed = 3490) %>%
    mutate(period=i,ont="MF")
  gsea_cc<-gseGO(geneList = gene_fc,ont = "CC",OrgDb = org.Hs.eg.db,keyType = "SYMBOL", 
                 pAdjustMethod = "BH",pvalueCutoff = 0.5,seed = 3490) %>%
    mutate(period=i,ont="CC")
  

  
  gsea_deg_go<-rbind(gsea_deg_go,gsea_bp@result,gsea_mf@result,gsea_cc@result)
}



##plot
gsea_deg_go<-read_excel("gsea_deg_go_251230.xlsx")

gsea_deg_go<-gsea_deg_go %>%
  mutate(period=gsub("1","1st",period),period=gsub("2","2nd",period),
         period=gsub("3","3rd",period),period=gsub("all","All Period",period)) %>%
  arrange(p.adjust) %>%
  dplyr::filter(p.adjust<0.05) %>%
  dplyr::mutate(log10padj= -log10(p.adjust))

#gsea_deg_go<-read_excel("../../../../../文章/多组学/4.supplymentary_cfRNA_GSEA.xlsx") %>%
#  mutate(pathway=gsub("\\(KEGG\\)","",pathway))

data_plot<-gsea_deg_go %>%
  dplyr::filter(p.adjust<0.001)
data_plot<-data_plot[1:50,]
pathway_level=as.factor(unique(data_plot$Description))

data_plot <- data_plot %>%
  group_by(ont) %>%                       # 按 ont 分组
  mutate(
    Description = factor(Description,     # 重新设定 levels
                         levels = rev(unique(Description[order(match(Description, pathway_level))])))
  ) %>% ungroup()


p<-ggplot(data_plot,aes(x=period ,y=Description,size=log10padj,fill=NES))+
  geom_point(shape=21,color=rgb(0,0,0,0))+
  scale_size_continuous(range = c(3,7))+
  #scale_y_discrete(limits=rev(pathway_level))+
  scale_fill_gradient(low="#62B197",high="#E18E6D")+
  ylab("Pathways")+
  xlab("Trimester")+
  labs(size=TeX('$-log_{10}(p.adj)$',bold = T))+
  facet_grid(ont ~ ., scales = "free_y", space = "free_y")+
  theme(axis.text.x = element_text(color="black",face="bold",size=11),
        axis.text.y = element_text(color="black",size=11),
        axis.title.x = element_text(color="black",face="bold",size=14),
        axis.title.y = element_text(color="black",face="bold",size=14),
        axis.title = element_text(face="bold"),
        plot.title = element_text(color="black",face="bold"),
        panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(color="grey92",linetype=1,linewidth = 0.2),
        panel.grid.major.x=element_line(color="grey92",linetype=1,linewidth = 0.2),
        legend.key = element_rect(fill = "white",colour = NA),
        legend.text = element_text(color="black",size=14),
        legend.title = element_text(color="black",face="bold",size=14))
