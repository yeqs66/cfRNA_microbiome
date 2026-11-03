library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(openxlsx)
library(readxl)
library(clusterProfiler)
library(tidyverse)
library(KEGG.db)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(gridExtra)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)

rm(list=ls())
oridir<-getwd()
setwd("__PATH_001__")

load("../../../GDM VS control/临床信息/labelledclinical_20240223.rda")
load("../select_feature.rda")
load("../KEGG_db_240228_ENTREZID&SYMBOL.rda")  #使用createKEGGdb自建KEGG_db数据库
load("../../../../../../数据库/pathways_ENSG.rda")
#load("../class_custom.rda") #自定义KEGG class数据库
load("../keggclass_custom.rda") #自定义KEGG class数据库
deg_40<-read.csv("__PATH_002__",header = T) %>%
  dplyr::filter(adj.P.Val<0.1)
deg_40<-as.character(deg_40$SYMBOL)


module_gene<-read.csv("__PATH_003__",header = T) %>%
  dplyr::rename(SYMBOL=X)

##### set plot para #####
fontsize_row_1<-4
fontsize_row_2<-8

fontsize_col<-5

fontsize_legend_title<-9
fontsize_legend_label<-8

##### set plot para #####

##### prepare kegg_gmt #####
#kegg_pathway_hierarchy层级文件
kegg_pathway_hierarchy<-read_excel("../../../../../../../GSEA gmt文件/msigdb_v2023.2.Hs_files_to_download_locally/kegg pathway层级关系-msidgb改名.xlsx")
kegg_gmt<-kegg_db[,c(1,2,4)] %>%
  dplyr::rename(term_id=ID,gene=SYMBOL) %>%
  left_join(kegg_pathway_hierarchy[,-6],by="term_id")

pathways <- do.call(rbind, lapply(seq_along(pathways$transciptome_ENTREZID), function(i) {
       data.frame(term=rep(names(pathways$transciptome_ENTREZID[i]),length(pathways$transciptome_ENTREZID[[i]])), value=pathways$transciptome_ENTREZID[[i]], stringsAsFactors=FALSE)
     })) %>%
       dplyr::filter(grepl("(KEGG)",term)) %>%
  mutate(term=gsub("\\(KEGG\\) ","",term))


symbol_ensg<-read.table("../../../GDM VS control/病例对照200对/新质控测试/基检数质控+样本质控/limma流程/new_matrix_limma/m_miRNA_SYMBOL_ENSG.txt",header = T)

#tmp<-AnnotationDbi::select(org.Hs.eg.db,keys = tmp,columns = "ENSEMBL",keytype = "SYMBOL")

##### prepare kegg_gmt #####

##### module gene number #####

module_gene_num<-as.numeric()
for (i in unique(module_gene$moduleColors)) {
  module_gene_num[i]<-length(module_gene$SYMBOL[which(module_gene$moduleColors==i)])
}
module_gene_num<-data.frame(num=module_gene_num) %>%
  dplyr::arrange(desc(num)) %>%
  mutate(order=1:length(unique(module_gene$moduleColors)),module=rownames(.),order=as.character(order)) %>%
  dplyr::rename(moduleColors=module)
rownames(module_gene_num)<-module_gene_num$order
#module_gene_num<-module_gene_num[arrange_module,]

#module_gene_num<-module_gene_num[match(arrange_module,names(module_gene_num))]
module_gene<-module_gene %>% 
  left_join(module_gene_num,by="moduleColors")

##### module gene number #####



##### fisher's test for DEG #####
modules_fisher_deg<-as.numeric()
for (i in unique(module_gene$order)) {
  module_deg<-as.data.frame(matrix(ncol=2,nrow = 2,data = 0))
  module_deg[1,1]<-as.numeric(length(intersect(module_gene$SYMBOL[which(module_gene$order==i)],deg_40)))
  module_deg[1,2]<-as.numeric(length(deg_40) - as.numeric(module_deg[1,1]))
  module_deg[2,1]<-as.numeric(length(module_gene$SYMBOL[which(module_gene$order==i)])-module_deg[1,1])
  module_deg[2,2]<-as.numeric(length(module_gene$SYMBOL) - length(deg_40))
  test_result <- fisher.test(module_deg)
  modules_fisher_deg[i]<-test_result$p.value
}

modules_fisher_deg<-data.frame(module=names(modules_fisher_deg),pvalue=modules_fisher_deg,padj=p.adjust(modules_fisher_deg,method = "BH",length(modules_fisher_deg))) %>%
  mutate(pch=ifelse(padj < 0.001, "***", ifelse(padj<0.01,"**", ifelse(padj < 0.05, "*", ""))))
#加**

#modules_fisher_deg<-modules_fisher_deg[arrange_module,]

##### fisher's test for DEG #####



##### pick 10 top gene set in each cluster #####
module_enrich<-list()
for (i in unique(module_gene$order)) {
  tmp<-as.character(module_gene$SYMBOL[which(module_gene$order==i)])
  tmp<-gsub("[3,5]p$","",tmp)
  tmp<-gsub("hsa.miR.","MIR",tmp)
  tmp<-gsub("hsa.let","MIRLET",tmp)
  tmp<-gsub("MIR\\.","MIR",tmp)
  temp<-grep("[0-9]\\.[0-9]",tmp)
  tmp[temp]<-gsub("\\.","-",tmp[temp])
  tmp<-gsub("\\.","",tmp)
  #tmp<-bitr(geneID = tmp,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  tmp<-unique(tmp)
  tmp<-toupper(tmp)
  tmp<-data.frame(SYMBOL=tmp) %>%
    left_join(symbol_ensg,by="SYMBOL")
  
  module_enrich[[i]]<-tmp$ENSEMBL
}

formula_res_kegg <- compareCluster(
  module_enrich,
  fun = "enricher",
  TERM2GENE=pathways,
  #style="margin: 0px; padding: 0px; font-size: inherit; line-height: inherit; color: rgb(80, 161, 79); overflow-wrap: inherit !important; word-break: inherit !important;">"BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

cluster_kegg<-formula_res_kegg@compareClusterResult %>%
  mutate(setsize=gsub(".*\\/","",GeneRatio),
         GeneRatio_2=as.numeric(Count) / as.numeric(setsize),
         log10padj= -log10(p.adjust)) %>%
  dplyr::rename(term=Description)

cluster_kegg_flt<-cluster_kegg[which(cluster_kegg$p.adjust<0.001),]
#pathway_level<-unique(cluster_kegg_flt$term[order(cluster_kegg_flt$term,cluster_kegg_flt$GeneRatio_2,decreasing = T)])



cluster_kegg_flt_top15<-data.frame()
for (i in unique(cluster_kegg_flt$term)) {
  tmp<-cluster_kegg_flt[which(cluster_kegg_flt$term==i),] %>%
    arrange(p.adjust)
  if (nrow(tmp)>14) {
    tmp<-tmp[1:15,]
  }
  cluster_kegg_flt_top15<-rbind(cluster_kegg_flt_top15,tmp)
}

cluster_kegg_flt_top15<-cluster_kegg_flt_top15 %>%
  left_join(kegg_pathway_hierarchy[,c(3,6)],by="term") %>%
  arrange(class) %>%
  left_join(keggclass_custom,by="class") %>% ##自定义class
  #left_join(kegg_pathway_hierarchy[,c(4,6)]) %>%
  #dplyr::rename(class_reorder=class) %>%
  mutate(log10padj=as.numeric(log10padj))


##### para_heatmap_gene #####
heatmap_wgcna_enrichKEGG<-cluster_kegg_flt_top15[,c(1,16,3)] %>% 
  pivot_wider(names_from = term,values_from = log10padj)

#heatmap_wgcna_enrichKEGG<-cluster_kegg_flt_top15[,c(1,13,3)] %>% 
#  pivot_wider(names_from = term,values_from = log10padj,values_fn = list)
#tmp<-data.frame(matrix(nrow = 9, ncol = 0))
#for (i in names(heatmap_wgcna_enrichKEGG)[-1]) {
#  temp<-gsub("NULL",0,heatmap_wgcna_enrichKEGG[[i]])
#  temp<-as.numeric(temp)
#  names(temp)<-heatmap_wgcna_enrichKEGG$Cluster
#  print(is.na(temp))
#  print(i)
#  tmp<-cbind(tmp,temp)
#}
#heatmap_wgcna_enrichKEGG<-tmp

heatmap_wgcna_enrichKEGG<-as.data.frame(heatmap_wgcna_enrichKEGG)
rownames(heatmap_wgcna_enrichKEGG)<-heatmap_wgcna_enrichKEGG$Cluster
heatmap_wgcna_enrichKEGG<-heatmap_wgcna_enrichKEGG %>%
  #dplyr::select(-Cluster) %>%
  dplyr::select(cluster_kegg_flt_top15$term)
heatmap_wgcna_enrichKEGG[is.na(heatmap_wgcna_enrichKEGG)]<-0


#coul <- colorRampPalette(brewer.pal(11, "PuBu"))(50)

anno_col_wgcna_heatmap<-list()
for (i in 1:length(unique(cluster_kegg_flt_top15$class_reorder))) {
  #  anno_col_wgcna_heatmap[["pathway_class"]][unique(cluster_kegg_flt_top15$class)[i]]<-paletteer::paletteer_d("trekcolors::lcars_series",length(unique(cluster_kegg_flt_top15$class)))[i]
  anno_col_wgcna_heatmap[["pathway_class"]][unique(cluster_kegg_flt_top15$class_reorder)[i]]<-colorRampPalette(brewer.pal(8, "Set3"))(length(unique(cluster_kegg_flt_top15$class_reorder)))[i]
}

anno_pathway<-data.frame(term=colnames(heatmap_wgcna_enrichKEGG)) %>%
  left_join(unique(cluster_kegg_flt_top15[,c(3,15)]),by="term")

anno_pathway<-HeatmapAnnotation(pathway_class=factor(anno_pathway$class_reorder,levels = unique(anno_pathway$class_reorder)),
                                #title_name=anno_text(c(rep("",round(length(colnames(heatmap_wgcna_enrichKEGG))/2-0.1)),"Pathway Class",rep("",round(length(colnames(heatmap_wgcna_enrichKEGG))/2-0.1))),
                                #                     rot=90,
                                #                     location=1),
                                col = anno_col_wgcna_heatmap["pathway_class"],
                                #annotation_label=NULL,
                                show_annotation_name = F,
                                show_legend = T,
                                which = "row",
                                #annotation_name_gp = gpar(fontsize=4),
                                annotation_legend_param = list(
                                  pathway_class=list(
                                    title = "Pathway Class",
                                    title_gp = gpar(fontsize=fontsize_legend_title,fontface="bold"),
                                    labels_gp = gpar(fontsize = fontsize_legend_label,fontface="plain"),
                                    grid_height=unit(2,units = "mm"),
                                    grid_width=unit(2,units = "mm")
                                  )
                                )
)


#anno_pathway<-Legend(labels = factor(anno_pathway$class_reorder,levels = unique(anno_pathway$class_reorder)),
#                                legend_gp = gpar(anno_col_wgcna_heatmap["pathway_class"]),
#                                title ="Pathway Class",
#                                direction  = "horizontal"
#)

heatmap_palette_pathway<-colorRamp2(c(0,5,10,60),c("white","#6BAED6","#4292C6","#08519C"))

heatmap_wgcna_enrichKEGG<-as.data.frame(t(heatmap_wgcna_enrichKEGG)) 

tmp<-as.data.frame(matrix(ncol=length(setdiff(as.character(module_gene_num$order),names(heatmap_wgcna_enrichKEGG))),nrow = nrow(heatmap_wgcna_enrichKEGG),data = 0)) %>%
  dplyr::rename_all(~ setdiff(as.character(module_gene_num$order),names(heatmap_wgcna_enrichKEGG)))

heatmap_wgcna_enrichKEGG<-cbind(heatmap_wgcna_enrichKEGG,tmp)

module_gene_num<-module_gene_num[as.character(names(heatmap_wgcna_enrichKEGG)),]
modules_fisher_deg<-modules_fisher_deg[as.character(names(heatmap_wgcna_enrichKEGG)),]

#heatmap_wgcna_enrichKEGG<-heatmap_wgcna_enrichKEGG %>%
#  dplyr::rename(!!!setNames(module_gene_num$moduleColors, module_gene_num$order))


#is_sig_0.01 = modules_fisher_deg$padj < 0.01
#is_sig_0.05 = 0.01 < modules_fisher_deg$padj < 0.05
#pch = 

## set color
col_gene_deg<-colorRamp2(c(1,0.5,0),c("#EFEDF5","#DADAEB","#9E9AC8"))

## anno for heatmap
anno_gene_num<-HeatmapAnnotation(
  module_name_1=anno_text(module_gene_num$order,gp = gpar(fontsize=fontsize_col),rot = 0),
  #gene_num=anno_barplot(as.numeric(module_gene_num$num),
  #                      gp = gpar(fill = "#5F559BFF",
  #                                col = NA),
  #                      height = unit(2,"cm")),
  "Gene num"=anno_lines(as.numeric(module_gene_num$num),
                      gp=gpar(col="#74C476",alpha=0.5),
                      add_points = T,
                      pch=21,
                      size = unit(1,"mm"),
                      pt_gp = gpar(fill="#74C476",alpha=0.8,pch=21),
                      axis_param = list(gp=gpar(fontsize=fontsize_row_2))),
  "DEG"=anno_simple(modules_fisher_deg$padj,
                       col = col_gene_deg,
                       pch = modules_fisher_deg$pch,
                       pt_size = unit(0.3,"npc"),
                       border = T),
  #gp=gpar(fontsize=4),
  annotation_name_side="right",
  #col=list(gene_deg=col_gene_deg),
  gap = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize=fontsize_row_2)
)

## padj's sig sign for gene_deg 
lgd_deg_pvalue = Legend(title = "Enrichment(FDR,Fisher test)", col_fun = col_gene_deg, at = c(0, 0.5, 1), 
                        labels = c("0", "0.5",  "1"),title_gp = gpar(fontsize=fontsize_legend_label,fontface="bold"),
                        direction = "horizontal",labels_gp = gpar(fontsize=fontsize_legend_title,fontface="plain"),
                        #title_gap = unit(0.5,units = "mm"),
)

lgd_deg_sig_0.001 = Legend(pch = "*", type = "points", labels = "FDR < 0.05",size = unit(5,"mm"), background  = "white",labels_gp = gpar(fontsize=fontsize_legend_label,fontface="plain"),title_gap = unit(0.5,units = "mm"),)
lgd_deg_sig_0.01 = Legend(pch = "**", type = "points", labels = "FDR < 0.01",size = unit(5,"mm"), background  = "white",labels_gp = gpar(fontsize=fontsize_legend_label,fontface="plain"),title_gap = unit(0.5,units = "mm"),)
lgd_deg_sig_0.05 = Legend(pch = "***", type = "points", labels = "FDR < 0.001",size = unit(5,"mm"), background  = "white",labels_gp = gpar(fontsize=fontsize_legend_label,fontface="plain"),title_gap = unit(0.5,units = "mm"),)

#lgd_gene_deg = Legend(title = "Enrichment(FDR,Fisher test)",
#                      col_fun=colorRamp2(c(1,0.5,0),c("#EFEDF5","#DADAEB","#9E9AC8")),
#                      at=(c(0,0.5,1)),
#                      labels = c(0,0.5,1),direction = "horizontal")

ht_gene<-ComplexHeatmap::Heatmap(heatmap_wgcna_enrichKEGG,
                                 col = heatmap_palette_pathway,
                                 cluster_rows = F,cluster_columns = T,
                                 column_names_gp = gpar(fontsize=fontsize_col),
                                 show_heatmap_legend = T,
                                 heatmap_legend_param = list(title = "-log10(FDR)",title_gp = gpar(fontsize=fontsize_legend_title,fontface="bold"),
                                                             labels_gp = gpar(fontsize = fontsize_legend_label,fontface="plain")),
                                 row_names_side = "left",
                                 show_row_names = T, #显示rownames
                                 row_names_gp = gpar(fontsize=fontsize_row_1),
                                 column_title_rot = 90,
                                 column_names_rot = 0,
                                 height = unit(12, "cm"),
                                 width = unit(4,units = "cm"),
                                 #heatmap_height = unit(1, "npc"),
                                 #heatmap_width = unit(1, "npc"),
                                 right_annotation = anno_pathway,
                                 bottom_annotation = anno_gene_num,
                                 column_names_side = "top",
                                 border = T,
                                 name = "ht_gene",
                                 #column_split = T
                                 #row_title_gp = gpar(
                                 #  margin=unit(c(2,0,0,0),"mm")
                                 #)
) ## 需要添加cluster间的间隔线

##### para_heatmap_gene #####


##### para_heatmap_clinical #####
#load("../../代谢组/wgcna/wgcna merge-cfRNA &  metabolite/test/coldata_wgcna_0709.rda")
load("__PATH_004__")
feature_clinical<-data.frame(term=select_feature$clinical,class=gsub(".*_","",select_feature$clinical)) %>%
  mutate(term_new=term,
         term_new=gsub("UA\\.x","UA",term_new),
         term_new=gsub("Age_mo_Anthropometric","Age",term_new),
         term_new=gsub("height_mo_Anthropometric","Height of mother",term_new),
         term_new=gsub("weight_mo_first_check_Anthropometric","Prepregnancy weight of mother",term_new),
         term_new=gsub("BMI_mo_Anthropometric","Prepregnancy BMI",term_new),
         term_new=gsub("SBP_first_check_BloodPressuse","SBP (1st trimester)",term_new),
         term_new=gsub("DBP_first_check_BloodPressuse","DBP (1st trimester)",term_new),
         term_new=gsub("SBP_42d_BloodPressuse","SBP (42 days)",term_new),
         term_new=gsub("DBP_42d_BloodPressuse","DBP (42 days)",term_new),
         term_new=gsub("OGTT_0h_Glycaemic","OGTT 0h",term_new),
         term_new=gsub("OGTT_1h_Glycaemic","OGTT 1h",term_new),
         term_new=gsub("OGTT_2h_Glycaemic","OGTT 2h",term_new),
         term_new=gsub("OGTT_HbA1c_Glycaemic","OGTT HbA1c",term_new),
         term_new=gsub("birth_BMI_OffspringRelated","BMI of child at birth",term_new),
         term_new=gsub("delivery_mode_kid_OffspringRelated","Delivery mode",term_new),
         #term_new=gsub("^([a-z])",toupper,term_new,perl=T),
         #term_new=str_to_title(term_new,),
         term_new=gsub("_[^_]*$","",term_new),
         term_new=gsub("_"," ",term_new),
         term_new=gsub("TyG$","TyG index",term_new),
         term_new=gsub("TyG-BMI","TyG-BMI index",term_new),
         term_new=gsub("TG/HDL-C","TG/HDL-C ratio",term_new),
         class=gsub("kidney","Kidney",class))

wgcna_rna_cor<-read.csv("__PATH_005__",row.names = 1,check.names = F)
wgcna_rna_pvalue<-read.csv("__PATH_006__",row.names = 1,check.names = F)
#wgcna_rna_cor<-wgcna_rna_cor[-which(rownames(wgcna_rna_cor)=="seqff"),]
#rownames(wgcna_rna_cor)<-colnames(coldata_wgcna)

wgcna_rna_cor<-wgcna_rna_cor[which(rownames(wgcna_rna_cor) %in% intersect(rownames(wgcna_rna_cor),colnames(coldata_wgcna))),]
wgcna_rna_pvalue<-wgcna_rna_pvalue[which(rownames(wgcna_rna_pvalue) %in% intersect(rownames(wgcna_rna_pvalue),colnames(coldata_wgcna))),]


wgcna_rna_cor_flt<-na.omit(wgcna_rna_cor)
wgcna_rna_cor_flt<-wgcna_rna_cor_flt[apply(wgcna_rna_cor_flt,1,function(x) max(x) > 0.10| min(x) < -0.10),]
wgcna_rna_cor_flt<-wgcna_rna_cor_flt[-which(rownames(wgcna_rna_cor_flt)=="seqff"),] %>%
  dplyr::rename_with(~ gsub("ME", "", .x)) %>%
  dplyr::rename(!!!setNames(module_gene_num$module,module_gene_num$order))

wgcna_rna_cor_flt_clinical<-wgcna_rna_cor_flt[which(rownames(wgcna_rna_cor_flt) %in% feature_clinical$term),colnames(heatmap_wgcna_enrichKEGG)]

wgcna_rna_pvalue_flt<-wgcna_rna_pvalue %>%
  dplyr::rename_with(~ gsub("ME", "", .x)) %>%
  dplyr::rename(!!!setNames(module_gene_num$module,module_gene_num$order))
wgcna_rna_pvalue_flt_clinical<-wgcna_rna_pvalue_flt[rownames(wgcna_rna_cor_flt_clinical),colnames(wgcna_rna_cor_flt_clinical)]

## set clinical class color
for (i in 1:length(unique(feature_clinical$class))) {
  anno_col_wgcna_heatmap[["clinical_class"]][unique(feature_clinical$class)[i]]<-ggsci::pal_frontiers()(length(unique(feature_clinical$class)))[i]
  
}

anno_clinical<-HeatmapAnnotation(clinical_class=factor(feature_clinical$class[which(feature_clinical$term %in% rownames(wgcna_rna_cor_flt_clinical))],levels = unique(feature_clinical$class[which(feature_clinical$term %in% rownames(wgcna_rna_cor_flt_clinical))])),
                                 #title_name=anno_text(c(rep("",round(length(rownames(wgcna_rna_cor_flt_clinical))/2-0.1)),"Clinical Class",rep("",round(length(rownames(wgcna_rna_cor_flt_clinical))/2-0.1))),
                                 #                     rot=90,
                                 #                    location=1),
                                 col = anno_col_wgcna_heatmap["clinical_class"],
                                 #annotation_label=NULL,
                                 show_annotation_name = F,
                                 show_legend = T,
                                 #annotation_name_gp = gpar(fontsize=4),
                                 which = "row",
                                 annotation_legend_param = list(
                                   clinical_class=list(
                                     title = "Clinical Class",
                                     title_gp = gpar(fontsize=fontsize_legend_title,fontface="bold"),
                                     labels_gp=gpar(fontsize=fontsize_legend_label,fontface="plain"),
                                     grid_height=unit(2,units = "mm"),
                                     grid_width=unit(2,units = "mm")
                                     #legend_height=unit(2,units = "cm")
                                   )
                                 ))

anno_col_clinical<-HeatmapAnnotation(
  module_name_2=anno_text(module_gene_num$order,gp = gpar(fontsize=fontsize_col),rot = 0)
)

heatmap_palette_clinical<-colorRamp2(c(0.5,0,-0.5),c("#6D3580","white","#FFE26F"))
#heatmap_palette_clinical<-colorRamp2(c(0.5,0,-0.5),c("FFF7AC",""))

#rownames(wgcna_rna_cor_flt_clinical)<-gsub("_[^_]*$","",rownames(wgcna_rna_cor_flt_clinical))
#rownames(wgcna_rna_cor_flt_clinical)<-gsub("_"," ",rownames(wgcna_rna_cor_flt_clinical))
#rownames(wgcna_rna_cor_flt_clinical)<-gsub("UA\\.x","UA",rownames(wgcna_rna_cor_flt_clinical))

wgcna_rna_cor_flt_clinical<-wgcna_rna_cor_flt_clinical %>%
  t() %>% as.data.frame() %>%
  dplyr::rename(!!!setNames(feature_clinical$term[which(feature_clinical$term %in% names(.))],
                            feature_clinical$term_new[which(feature_clinical$term %in% names(.))])) %>%
  t() %>% as.data.frame()

ht_clinical<-ComplexHeatmap::Heatmap(wgcna_rna_cor_flt_clinical,
                                     col = heatmap_palette_clinical,
                                     show_heatmap_legend = T,
                                     cluster_rows = F,cluster_columns = F,
                                     row_names_side = "left",
                                     show_column_names = F,
                                     show_row_names = T, #显示rownames
                                     row_names_gp = gpar(fontsize=fontsize_row_2),
                                     heatmap_legend_param = list(
                                       title = "Correlation",
                                       title_gp=gpar(fontsize=fontsize_legend_title,fontface="bold"),
                                       labels_gp=gpar(fontsize=fontsize_legend_label,fontface="plain")),
                                     column_names_gp = gpar(fontsize=fontsize_col),
                                     column_title_rot = 90,
                                     column_names_rot = 0,
                                     right_annotation = anno_clinical,
                                     top_annotation = anno_col_clinical,
                                     column_names_side = "top",
                                     height = unit(4/15*nrow(wgcna_rna_cor_flt_clinical),units = "cm"),
                                     width = unit(4,units = "cm"),
                                     #heatmap_height = unit(1, "npc"),
                                     #heatmap_width = unit(1, "npc"),
                                     border = T,
                                     name = "ht_clinical",
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       if (wgcna_rna_pvalue_flt_clinical[i,j]<0.05 & wgcna_rna_pvalue_flt_clinical[i,j]>0.01) {
                                         grid.rect(x = x, y = y,width = width,height = height,gp=gpar(lwd=1,fill="transparent"))
                                         grid.text("*", x , y , gp=gpar(fontsize=3))
                                       }
                                       if (wgcna_rna_pvalue_flt_clinical[i,j]<0.01 & wgcna_rna_pvalue_flt_clinical[i,j]>0.001) {
                                         grid.rect(x = x, y = y,width = width,height = height,gp=gpar(lwd=1,fill="transparent"))
                                         grid.text("*", x , y , gp=gpar(fontsize=3))
                                       }
                                       if (wgcna_rna_pvalue_flt_clinical[i,j]<0.001) {
                                         grid.rect(x = x, y = y,width = width,height = height,gp=gpar(lwd=1,fill="transparent"))
                                         grid.text("***", x , y , gp=gpar(fontsize=3))
                                       }
                                     }
                                     #column_split = T
                                     #row_title_gp = gpar(
                                     #  margin=unit(c(2,0,0,0),"mm")
                                     #)
)



##### para_heatmap_clinical #####


##### plot #####
ht_list<-ht_gene %v% ht_clinical

cairo_pdf("__PATH_007__", width = 10, height = 10, 
          family = "Arial",  # 使用标准字体
          onefile = FALSE) 
draw(ht_list,annotation_legend_list=list(lgd_deg_pvalue ,lgd_deg_sig_0.001, lgd_deg_sig_0.01, lgd_deg_sig_0.05))
dev.off()


###分开画
cairo_pdf("__PATH_008__",height = 10,width = 10,family = "Arial",onefile = F)
draw(ht_gene,merge_legend = FALSE)
dev.off()

cairo_pdf("__PATH_009__",height = 10,width = 10,family = "Arial",onefile = F)
draw(ht_clinical,merge_legend = F)
dev.off()

#draw(ht_list)
#output
setwd(oridir)
