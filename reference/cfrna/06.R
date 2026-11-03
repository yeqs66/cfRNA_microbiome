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

load("__PATH_001__")
#pathways <- pathways$transciptome_ENTREZID[grep("(KEGG)",names(pathways$transciptome_ENTREZID))]
pathways <- do.call(rbind, lapply(seq_along(pathways$transciptome_ENTREZID), function(i) {
  data.frame(name=rep(names(pathways$transciptome_ENTREZID[i]),length(pathways$transciptome_ENTREZID[[i]])), value=pathways$transciptome_ENTREZID[[i]], stringsAsFactors=FALSE)
})) %>%
  dplyr::filter(grepl("(KEGG)",name))
symbol_ensg<-read.table("__PATH_002__")

##### vs deg #####
#for (n in 1:100) {
deg_mRNA<-read.csv("limma_349_batch_only.csv",header=T)
  gsea_deg_reactomekegg<-data.frame()
  for (i in unique(deg_mRNA$period)) {
    #temp<-read_excel("__PATH_003__",sheet = i) %>%
    temp<-deg_mRNA %>%
      dplyr::filter(period==i) %>%
      arrange(desc(logFC))
    gene_fc<-temp$logFC
    names(gene_fc)<-temp$SYMBOL
    
    symbol_ensg_flt<-symbol_ensg[match(names(gene_fc),symbol_ensg$V2),]
    all(symbol_ensg_flt$V2==names(gene_fc))
    names(gene_fc)<-symbol_ensg_flt$V1
    
    set.seed(3490)
    gsea_other<-GSEA(geneList = gene_fc,TERM2GENE = pathways,pvalueCutoff = 1,pAdjustMethod = "BH")
    gsea_other<-as.data.frame(gsea_other@result) %>%
      mutate(period=i)
    #gsea_other<-fgsea(pathways, stats=gene_fc) %>%
    #  mutate(period=i)
    
    #gsea_other<-GSEA(geneList = gene_fc,TERM2GENE = pathways,pvalueCutoff = 1,pAdjustMethod = "BH")
    
    #gsea_other<-as.data.frame(gsea_other@result) %>%
    #  mutate(type="other",period=i)
    
    gsea_deg_reactomekegg<-rbind(gsea_deg_reactomekegg,gsea_other)
  }
  
#  if (length(which(gsea_deg_reactome$padj[grep("porph",gsea_deg_reactomekegg$pathway)]<0.05))>1) {
#    print(n)
#    break
#  } else {
#    next
#  }
#}


#gsea_all<-read_excel("../../GDM VS control/病例对照200对/新质控测试/基检数质控+样本质控/limma流程/new_matrix_limma/gsea_limma_merge_nocoef_all.xlsx") %>%
#  dplyr::filter(p.adjust<0.01) %>%
#  dplyr::mutate(log10padj= -log10(p.adjust))

gsea_deg_reactomekegg<-gsea_deg_reactomekegg %>%
  mutate(period=gsub("1","1st",period),period=gsub("2","2nd",period),
         period=gsub("3","3rd",period),period=gsub("all","All Period",period)) %>%
  arrange(padj) %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::mutate(log10padj= -log10(padj)) %>%
  dplyr::rename(Description=pathway) %>%
  mutate(Description=gsub("\\(KEGG\\)","",Description))

#gsea_deg_reactomekegg<-read_excel("../../../../../文章/多组学/4.supplymentary_cfRNA_GSEA.xlsx") %>%
#  mutate(pathway=gsub("\\(KEGG\\)","",pathway))

pathway_level=as.factor(unique(gsea_deg_reactomekegg$Description))

cairo_pdf("gsea_deg_kegg_260205.pdf",width = 10,height = 11,family = "Arial")
ggplot(gsea_deg_reactomekegg,aes(x=period ,y=Description,size=log10padj,fill=NES))+
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

dev.off()
#ggsave(filename = "正式/gsea_deg_vs_1.2.pdf",p,width = 10,height = 14)
write.xlsx(gsea_deg_reactomekegg,file = "__PATH_004__")


##### vs deg GO #####
gsea_deg_go<-data.frame()
for (i in c(10:13)) {
  temp<-deg_mRNA %>%
    arrange(desc(logFC))
  gene_fc<-temp$logFC
  names(gene_fc)<-temp$SYMBOL
  
  set.seed(3490)
  gsea_bp<-gseGO(gene_fc,ont = "BP", OrgDb = "org.Hs.eg.db",keyType ="SYMBOL",pvalueCutoff=1,seed = 10)
  gsea_mf<-gseGO(gene_fc,ont = "MF", OrgDb = "org.Hs.eg.db",keyType ="SYMBOL",pvalueCutoff=1,seed = 10)
  gsea_cc<-gseGO(gene_fc,ont = "CC", OrgDb = "org.Hs.eg.db",keyType ="SYMBOL",pvalueCutoff=1,seed = 10)

  gsea_bp<-as.data.frame(gsea_bp@result) %>%
    mutate(type="BP",period=i)
  gsea_mf<-as.data.frame(gsea_mf@result) %>%
    mutate(type="MF",period=i)
  gsea_cc<-as.data.frame(gsea_cc@result) %>%
    mutate(type="CC",period=i)
  
  gsea_deg_go<-rbind(gsea_deg_go,gsea_bp,gsea_mf,gsea_cc)
}


##plot
gsea_deg_go<-read_excel("gsea_deg_go_251230.xlsx")

gsea_deg_go<-gsea_deg_go %>%
  mutate(period=gsub("1","1st",period),period=gsub("2","2nd",period),
         period=gsub("3","3rd",period),period=gsub("all","All Period",period),
         Description=str_to_sentence(Description)) %>%
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

cairo_pdf("gsea_deg_go_260323.pdf",width = 13,height = 10,family = "Arial")

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
p

dev.off()

write.xlsx(gsea_deg_go,file = "__PATH_005__")


