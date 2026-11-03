library(ggplot2)
library(tidyverse)
library(circlize)
library(clusterProfiler)
library(fgsea)
library(ComplexHeatmap)

rm(list = ls())
setwd("__PATH_001__")

##### gsea #####
load("kogmt.rda")
load("klist_to_module_20240325_1524_updated.rda")
gmt <- merge(kogmt,klist_to_module,by=1) %>%
  dplyr::select(mod,V2)

# gmt for pathway
pathgmt <- read.table("ko00001.forR.KO_only.gmt", header=FALSE)

#gmt_merge<-pathgmt %>%
#  left_join(gmt) %>%
#  dplyr::select(-V2) %>%
#  unique()


gsea<-data.frame()
for (i in c("early","middle","late","overall")) {
  tmp<-read.csv(paste("new_deg/TJGDM1962s_ko_tri_95.zinb_",i,"_rpk.csv",sep = "")) %>%
    dplyr::filter(!log2MedianRatio %in% c(Inf,-Inf,NA)) %>%
    #na.omit() %>%
    arrange(desc(log2MedianRatio))
  gene_fc<-tmp$log2MedianRatio
  names(gene_fc)<-tmp$feature
  
  set.seed(123)
  gsea_tmp<-GSEA(gene_fc,pvalueCutoff = 1,pAdjustMethod = "BH",TERM2GENE = pathgmt,seed = T,minGSSize = 8)
  gsea_tmp<-gsea_tmp@result %>%
    mutate(period=i) %>%
    dplyr::select(-qvalue)
  
  gsea<-rbind(gsea,gsea_tmp)
  assign(paste("df_",i,sep = ""),gsea_tmp)
}
gsea<-gsea %>%
  arrange(pvalue)

pathway_p01<-gsea$ID[which(gsea$pvalue<0.1)] ##select sig pathways

##### plot #####
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
mypal <- colorRampPalette( c( "#62B197", "#E18E6D" ) )(1000)

colnames(df_early) <- paste("E_",colnames(df_early),sep="")
colnames(df_middle) <- paste("M_",colnames(df_middle),sep="")
colnames(df_late) <- paste("L_",colnames(df_late),sep="")
colnames(df_overall) <- paste("A_",colnames(df_overall),sep="")

bg <- aggregate(pathgmt$V2, by = list(pathgmt$V1), length)
colnames(bg) <- c("Path","count")
df <- merge(df_early,df_middle,by=0, all=T) %>%
  merge(.,df_late,by.x=1,by.y=0,all=T) %>%
  merge(.,df_overall,by.x=1,by.y=0,all=T) %>%
  merge(.,bg,by=1)
df[is.na(df)] <- 0
df <- df %>%
  dplyr::mutate(E_color=map2color(E_NES,mypal),
                M_color=map2color(M_NES,mypal),
                L_color=map2color(L_NES,mypal),
                A_color=map2color(A_NES,mypal),
                #E_ratio=E_setSize/count,
                #M_ratio=M_setSize/count,
                #L_ratio=L_setSize/count,
                #A_ratio=A_setSize/count,
                E_log10padj= -log10(E_p.adjust),
                M_log10padj= -log10(M_p.adjust),
                L_log10padj= -log10(L_p.adjust),
                A_log10padj= -log10(A_p.adjust),
                E_ratio=ifelse(E_pvalue<0.1,gsub(".*list=","",E_leading_edge),0),
                E_ratio=as.numeric(gsub("%.*","",E_ratio)) / 100,
                M_ratio=ifelse(M_pvalue<0.1,gsub(".*list=","",M_leading_edge),0),
                M_ratio=as.numeric(gsub("%.*","",M_ratio)) / 100,
                L_ratio=ifelse(L_pvalue<0.1,gsub(".*list=","",L_leading_edge),0),
                L_ratio=as.numeric(gsub("%.*","",L_ratio)) / 100,
                A_ratio=ifelse(A_pvalue<0.1,gsub(".*list=","",A_leading_edge),0),
                A_ratio=as.numeric(gsub("%.*","",A_ratio)) / 100
  ) %>%
  arrange(desc(Row.names)) %>%
  dplyr::filter(Row.names %in% pathway_p01)



#save(df,file="__PATH_002__")

library(circlize)
cairo_pdf("GSEA_test.pdf", width=6, height=7,family="Arial")
circos.par(start.degree = 270, gap.degree = 180,
           track.height=0.18,circle.margin=c(0.1,0.1))
circos.initialize(letters[1], xlim = c(0, 30))
circos.track(ylim = c(0, 0.6), panel.fun = function(x, y) {
  circos.axis(sector.index = "a", h = "top", major.at = 1:30-0.5,
              labels = gsub("ko","KO",gsub("]","",gsub(".*PATH:","",df$Row.names))),
              minor.ticks = 0,
              labels.cex = 0.6,direction = "outside",
              labels.niceFacing = F,labels.facing="reverse.clockwise")
  circos.barplot(df$A_ratio, 1:30-0.5,col=df$A_color, border = "black")

  circos.yaxis(side = "left", sector.index = "a",
               labels.cex=0.5,at=c(0,0.3,0.6))
  circos.text(30, 0, "All Period", facing = "reverse.clockwise",
              adj = c(0.95, -1), cex = 0.7)
}, bg.border = NA)
circos.track(ylim = c(0, 0.6), panel.fun = function(x, y) {
  value = matrix(df$L_ratio,ncol=1)
  circos.barplot(df$L_ratio, 1:30-0.5,col=df$L_color, border = "black")

  circos.yaxis(side = "left", sector.index = "a",
               labels.cex=0.5,at=c(0,0.3,0.6))
  circos.text(30, 0, "Late", facing = "reverse.clockwise",
              adj = c(1.3, -1), cex = 0.7)
}, bg.border = NA)
circos.track(ylim = c(0, 0.6), panel.fun = function(x, y) {
  value = matrix(df$M_ratio,ncol=1)
  circos.barplot(df$M_ratio, 1:30-0.5,col=df$M_color, border = "black")

  circos.yaxis(side = "left", sector.index = "a",
               labels.cex=0.5,at=c(0,0.3,0.6))
  circos.text(30, 0, "Middle", facing = "reverse.clockwise",
              adj = c(1.1, -1), cex = 0.7)
}, bg.border = NA)
circos.track(ylim = c(0, 0.6), panel.fun = function(x, y) {
  circos.barplot(df$E_ratio, 1:30-0.5,col = df$E_color, border = "black")
  #    circos.axis(sector.index = "a", h = "bottom", major.at = 1:30-0.5,
  #    labels = gsub(".*PATH:","",df$Row.names), minor.ticks = 0,
  #    labels.cex = 0.5,direction = "inside",
  #    labels.niceFacing = F,labels.facing="reverse.clockwise")

  circos.yaxis(side = "left", sector.index = "a",
               labels.cex=0.5,at=c(0,0.3,0.6))
  circos.text(30, 0, "Early", facing = "reverse.clockwise",
              adj = c(1.3, -1), cex = 0.7)
}, bg.border = NA)

##add legend of NES
generate_color <- function(nes) {
  ifelse(nes > 0, "#62B197", "#E18E6D")  # 红/蓝
}

df$E_color <- generate_color(df$E_NES)
df$M_color <- generate_color(df$M_NES)
df$L_color <- generate_color(df$L_NES)
df$A_color <- generate_color(df$A_NES)


circos.clear()

##add legend
col_fun <- colorRamp2(c(-2, 2), c("#62B197", "#E18E6D"))
lgd <- Legend(col_fun = col_fun, 
              title = "NES",
              direction = "vertical",
              legend_width = unit(0.5, "cm"),
              legend_height = unit(0.5,"cm"),
              labels_gp = gpar(fontsize = 6),
              title_gp = gpar(fontsize = 7, fontface = "bold"),
              grid_height = unit(0.2,"cm"))

# 放置在底部中央
draw(lgd, x = unit(0.54, "npc"), y = unit(0.45, "npc"), just = c("right", "bottom"))

#legend("topright", 
#       legend = "p.adjust < 0.1", 
#       fill = NA, 
#       border = "black",
#       bty = "n",
#       cex = 0.7)


dev.off()


fortable <- data.frame("KO"=df$Row.names,"name"=df$Row.names)
fortable <- fortable %>%
  dplyr::mutate(KO = gsub("ko","KO",
                          gsub("]","",
                               gsub(".*PATH:","",KO))),
                name = gsub("_", " ",
                            gsub("^[0-9]*_","",
                                 gsub("_\\[PATH:.*","",name)))) %>%
  arrange(KO)

write.csv(fortable, 'GSEA_KO_plot_table_260226.csv')
gsea_ko<-gsea
save(gsea,file="gsea_deg_KO_260226.rda")

