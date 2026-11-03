library(MOFA2)
library(tidyverse)
library(readxl)
library(reticulate)
library(WGCNA)
library(latex2exp)
rm(list = ls())
setwd("__PATH_001__")

##### MOFA pipeline #####
load("__PATH_002__")

load("__PATH_003__")
load("__PATH_004__")
sample_meta<-read_excel("__PATH_005__",sheet = 10) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_meta=sampleid)
sample_cfRNA<-read_excel("__PATH_005__",sheet = 8) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_cfRNA=sampleid) %>%
  left_join(sample_meta[,c("sign","sampleid_meta")])
##ko rpm
rpm_meta<-read.csv("__PATH_006__",header = T,row.names = 1)
sample_cfRNA<-sample_cfRNA %>%
  dplyr::filter(sampleid_meta %in% names(rpm_meta))
rpm_meta<-rpm_meta %>%
  dplyr::select(sample_cfRNA$sampleid_meta) %>%
  dplyr::rename(!!!setNames(sample_cfRNA$sampleid_meta,sample_cfRNA$sign)) %>%
  t() %>% as.data.frame()
rpm_meta<-na.omit(rpm_meta[rownames(coldata_wgcna),])
missing_sample <- setdiff(rownames(coldata_wgcna), rownames(rpm_meta))
rpm_meta[missing_sample, ] <- 0
rpm_meta<-rpm_meta[rownames(coldata_wgcna),]
rpm_meta<-as.data.frame(log2(rpm_meta+1))
all(rownames(rpm_meta)==rownames(coldata_wgcna))


##vst_cfRNA
load("__PATH_007__")
#vst_cfRNA<-as.data.frame(t(vst_cfRNA))
##long data frame
tmp_meta<-rpm_meta %>%
  rownames_to_column() %>%
  dplyr::rename(sample=rowname) %>%
  pivot_longer(cols = !sample,names_to = "feature",values_to = "value") %>%
  mutate(view="KO",
         value=as.numeric(value)) %>%
  dplyr::select(sample,feature,view,value)
tmp_rna<-vst_cfRNA %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(sample=rowname) %>%
  pivot_longer(cols = !sample,names_to = "feature",values_to = "value") %>%
  mutate(view="rna",
         value=as.numeric(value)) %>%
  dplyr::select(sample,feature,view,value)
MOFAobject<-rbind(tmp_rna,tmp_meta)


MOFAobject <- create_mofa(MOFAobject)
print(MOFAobject)
plot_data_overview(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
head(model_opts)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose <- T
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
MOFAobject.trained <- run_mofa(MOFAobject, outfile = "__PATH_008__", use_basilisk=TRUE)

##### downstream analysis #####
model<-load_model("MOFA2/model_tpm_ko_topVarience50_260120.hdf5")
#load("__PATH_003__")
load("mofa_input_ko_m_mi_lncRNA_topVarience50_260120.rda")
#load("__PATH_003__")
coldata_mofa<-coldata %>%
  mutate(sample=rownames(.))
#Nsamples = sum(get_dimensions(model)[["N"]])
#sample_metadata <- data.frame(sample = samples_names(model)[[1]],replace = TRUE) %>%
#  left_join(coldata_mofa[,c("sample","disease","period")])
samples_metadata(model) <- coldata_mofa %>%
  dplyr::select(sample,disease,period,ALT_Liver,AST_Liver,Preterm_OffspringRelated) %>%
  mutate(group=period) %>%   ##group by period
  mutate(disease=gsub("0","control",disease),
         disease=gsub("1","GDM",disease))

head(get_variance_explained(model)$r2_total[[1]])
get_variance_explained(model)$r2_per_factor[[1]]



##plot varience
cairo_pdf(filename = "__PATH_009__",
          width = 6,height = 3,family = "Arial")

p<-plot_variance_explained(model, x="view", y="factor",factors = c(1:5))
p+ facet_grid(~factor(group,
                      levels = c("1","2","3"),
                      labels = c("1st","2nd","3rd"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()



#load("__PATH_003__")
factors <- as.data.frame(get_factors(model, factors = "all",groups = 1,as.data.frame = F, scale = T)) 

weights <- get_weights(model, views = "all", factors = "Factor1",as.data.frame = T) %>%
  #dplyr::filter(factor=="Factor2") %>%
  mutate(value_abs=abs(value)) %>%
  arrange(factor,desc(value_abs)) %>%
  dplyr::filter(view=="mRNA")

plot_top_weights(model,
                 view = "KO",
                 factor = 1,
                 nfeatures = 20,
                 scale = T,
                 abs = F
)

lapply(weights,dim)

#factors<-as.data.frame(factors) %>%
#  rownames_to_column() %>%
#  dplyr::rename(sample=rowname) %>%
#  left_join(coldata_wgcna_mofa[which(coldata_wgcna_mofa$period==1),]) %>%
#  dplyr::select(-c(sample,period))

coldata_wgcna_mofa<-coldata[,c(1,4,14:61,65:85,87:129,132,134,138:191,195:199,201:205)]
coldata_wgcna_mofa<-coldata_wgcna_mofa[rownames(factors),c(3:45)] 

#factors<-factors %>%
#  rownames_to_column("sample") %>%
#  left_join(coldata_mofa[,c("sample","disease")],by="sample")

## factor cor traits
cor_mofa<-cor(x=coldata_wgcna_mofa,y=factors,use="p",method = "pearson")
cor_mofa_p<-as.data.frame(corPvalueStudent(cor_mofa,nSamples = nrow(factors)))
cor_mofa_padj<-data.frame(row.names = rownames(cor_mofa_p))

for (i in colnames(cor_mofa_p)) {
  tmp<-p.adjust(cor_mofa_p[[i]])
  cor_mofa_padj[[i]]<-tmp
}
#cor_mofa_padj<-na.omit(cor_mofa_padj)
tmp<-na.omit(cor_mofa_p[,1:5])
cor_mofa_p_sig<-data.frame(row.names = rownames(tmp))
for (row in rownames(tmp)) {
  if (any(tmp[row,] < 0.05)) {
    cor_mofa_p_sig <- rbind(cor_mofa_p_sig,tmp[row,])
  }
}

#tmp<-na.omit(cor_mofa_padj[,1:5])
#cor_mofa_padj_sig<-data.frame(row.names = rownames(tmp))
#for (row in rownames(tmp)) {
#  if (any(tmp[row,] < 0.05)) {
#    cor_mofa_padj_sig <- rbind(cor_mofa_padj_sig,tmp[row,])
#  }
#}

cor_mofa_padj_sig<-cor_mofa_padj[rownames(cor_mofa_p_sig),colnames(cor_mofa_p_sig)]

cor_mofa_sig<-cor_mofa[rownames(cor_mofa_p_sig),colnames(cor_mofa_p_sig)]

temp<-as.data.frame(cor_mofa) %>%
  dplyr::rename_with(~gsub("X..","",.x))
write.xlsx(temp,file="mofa_cor_pearson_trait_late.xlsx",rowNames=T)
temp<-cor_mofa_p %>%
  dplyr::rename_with(~gsub("X..","",.x))
write.xlsx(temp,file="mofa_cor_pearson_pvalue_trait_late.xlsx",rowNames=T)
temp<-cor_mofa_padj %>%
  dplyr::rename_with(~gsub("X..","",.x))
write.xlsx(temp,file="mofa_cor_pearson_pvalue_fdr_trait_late.xlsx",rowNames=T)

##### heatmap plot #####
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(paletteer) 

### set plot para
fontsize_row_1<-4
fontsize_row_2<-8

fontsize_col<-8

fontsize_legend_title<-9
fontsize_legend_label<-8

##prepare data
#data_cor_trait<- -log10(na.omit(cor_mofa_sig))

data_cor_trait<-as.data.frame(cor_mofa_sig) %>%
  dplyr::rename_with(~ gsub("X..", "", .x)) %>%
  dplyr::select(Factor1,Factor2,Factor3,Factor4,Factor5)




data_cor_trait<-data_cor_trait[rownames(cor_mofa_sig),]
#plot_data<-data_cor_trait[which(rownames(wgcna_rna_flt) %in% feature_clinical$term),colnames(heatmap_wgcna_enrichKEGG)]

feature_clinical<-data.frame(term=rownames(data_cor_trait),
                             class=gsub(".*_","",rownames(data_cor_trait))) %>%
  mutate(term=str_replace(term, "_[^_]+$", ""))

data_cor_trait<-data_cor_trait %>%
  rownames_to_column(var = "term") %>%
  #mutate(term=str_replace(term, "_[^_]+$", "")) %>%
  mutate(term=sub("_[^_]*$", "",term)) %>%
  column_to_rownames(var = "term")

## set clinical class color
anno_col_mofa_heatmap<-list()
for (i in 1:length(unique(feature_clinical$class))) {
  anno_col_mofa_heatmap[["clinical_class"]][unique(feature_clinical$class)[i]]<-ggsci::pal_npg(alpha = 0.6)(length(unique(feature_clinical$class)))[i]
  
}

anno_clinical<-HeatmapAnnotation(clinical_class=factor(feature_clinical$class[which(feature_clinical$term %in% rownames(data_cor_trait))],levels = unique(feature_clinical$class[which(feature_clinical$term %in% rownames(data_cor_trait))])),
                                 #title_name=anno_text(c(rep("",round(length(rownames(data_cor_trait))/2-0.1)),"Clinical Class",rep("",round(length(rownames(data_cor_trait))/2-0.1))),
                                 #                     rot=90,
                                 #                    location=1),
                                 col = anno_col_mofa_heatmap["clinical_class"],
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
  factor_name=anno_text(gsub("Factor","",colnames(data_cor_trait)),gp = gpar(fontsize=fontsize_row_2),rot = 0)
)

heatmap_palette_clinical<-colorRamp2(c(0.3,0,-0.3),c("#ECA799","white","#C3E3DA"))
#heatmap_palette_clinical<-colorRamp2(c(0.5,0,-0.5),c("FFF7AC",""))

ht_clinical<-ComplexHeatmap::Heatmap(data_cor_trait,
                                     col = heatmap_palette_clinical,
                                     show_heatmap_legend = T,
                                     cluster_rows = F,cluster_columns = F,
                                     row_names_side = "left",
                                     show_column_names = F,
                                     show_row_names = T, #显示rownames
                                     row_names_gp = gpar(fontsize=fontsize_row_2),
                                     heatmap_legend_param = list(
                                       title = "-log10 padj",
                                       title_gp=gpar(fontsize=fontsize_legend_title,fontface="bold"),
                                       labels_gp=gpar(fontsize=fontsize_legend_label,fontface="plain")),
                                     column_names_gp = gpar(fontsize=fontsize_col),
                                     column_title_rot = 90,
                                     column_names_rot = 90,
                                     right_annotation = anno_clinical,
                                     top_annotation = anno_col_clinical,
                                     column_names_side = "top",
                                     height = unit(4/15*nrow(data_cor_trait),units = "cm"),
                                     width = unit(4,units = "cm"),
                                     #heatmap_height = unit(1, "npc"),
                                     #heatmap_width = unit(1, "npc"),
                                     border = T,
                                     name = "ht_clinical",
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       if (cor_mofa_padj_sig[i,j]<0.05 & cor_mofa_padj_sig[i,j]>0.01) {
                                         grid.rect(x = x, y = y,width = width,height = height,gp=gpar(lwd=1,fill="transparent"))
                                         grid.text("*", x , y , gp=gpar(fontsize=6))
                                       }
                                       if (cor_mofa_padj_sig[i,j] <0.01 & cor_mofa_padj_sig[i,j]>0.001) {
                                         grid.rect(x = x, y = y,width = width,height = height,gp=gpar(lwd=1,fill="transparent"))
                                         grid.text("**", x , y , gp=gpar(fontsize=6))
                                       }
                                       if (cor_mofa_p_sig[i,j]<0.001) {
                                         grid.rect(x = x, y = y,width = width,height = height,gp=gpar(lwd=1,fill="transparent"))
                                         grid.text("***", x , y , gp=gpar(fontsize=6))
                                       }
                                     }
                                     #column_split = T
                                     #row_title_gp = gpar(
                                     #  margin=unit(c(2,0,0,0),"mm")
                                     #)
)

cairo_pdf(filename = "__PATH_010__",
          width = 5,height = 4,family = "Arial")
draw(ht_clinical)
dev.off()



##### heatmap to distingrate #####
plot_data_heatmap(model, 
                  view = "mRNA", 
                  factor = 1, 
                  features = 200,
                  groups = 1,
                  show_rownames = F, show_colnames = F, 
                  cluster_rows = T, cluster_cols = T,
                  annotation_samples = "disease"
)


##### gsea for each factor #####
library(clusterProfiler)
load("__PATH_011__")
pathways <- do.call(rbind, lapply(seq_along(pathways$transcriptome), function(i) {
  data.frame(term=rep(names(pathways$transcriptome[i]),length(pathways$transcriptome[[i]])), value=pathways$transcriptome[[i]], stringsAsFactors=FALSE)
})) %>%
  dplyr::filter(grepl("(KEGG)",term)) %>%
  mutate(term=gsub("\\(KEGG\\) ","",term))

gsea_mRNA<-data.frame()
gsea_KO<-data.frame()
for (i in 1:nrow(model@cache$variance_explained$r2_per_factor[[1]])) {   ##factor num
  weights <- get_weights(model, views = "mRNA", factors = "all",as.data.frame = T) %>%
    dplyr::filter(factor==paste0("Factor",i)) %>%
    arrange(desc(value))
  
  tmp<-weights$value
  names(tmp)<-weights$feature
  
  set.seed(3490)
  tmp_gsea_mRNA<-GSEA(tmp,TERM2GENE = pathways,pvalueCutoff = 1,pAdjustMethod = "BH",minGSSize = 8,seed = T)
  
  weights <- get_weights(model, views = "KO", factors = "all",as.data.frame = T) %>%
    dplyr::filter(factor==paste0("Factor",i)) %>%
    arrange(desc(value))
  
  tmp<-weights$value
  names(tmp)<-weights$feature
  tmp_gsea_KO<-GSEA(tmp,TERM2GENE = pathgmt,pvalueCutoff = 1,pAdjustMethod = "BH",seed = T,minGSSize = 8)
  
  tmp_gsea_mRNA<-as.data.frame(tmp_gsea_mRNA@result) %>%
    mutate(factor=i)
  tmp_gsea_KO<-as.data.frame(tmp_gsea_KO@result) %>%
    mutate(factor=i)
  
  gsea_mRNA<-rbind(gsea_mRNA,tmp_gsea_mRNA)
  gsea_KO<-rbind(gsea_KO,tmp_gsea_KO)
}


### plot
##KO
data<-gsea_KO %>%
  dplyr::filter(p.adjust<0.05) %>%
  mutate(factor=as.character(factor),
         log10padj= -log10(p.adjust),
         ID = gsub("ko","KO",
                   gsub("]","",
                        gsub(".*PATH:","",ID))),
         Description = gsub("_", " ",
                     gsub("^[0-9]*_","",
                          gsub("_\\[PATH:.*","",Description)))) %>%
  arrange(p.adjust)
pathway_level=as.factor(unique(data$Description))

p<-ggplot(data,aes(x=factor ,y=Description,size=log10padj,fill=NES))+
  geom_point(shape=21,stroke=0.04)+
  #scale_size_continuous(range = c(3,7))+
  scale_y_discrete(limits=rev(pathway_level))+
  scale_fill_gradient(low="#62B197",high="#E18E6D")+
  ylab("Pathways")+
  xlab("Factor")+
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

##cfRNA
data<-gsea_mRNA %>%
  dplyr::filter(p.adjust<0.005) %>%
  mutate(factor=as.character(factor),
         log10padj= -log10(p.adjust)) %>%
  arrange(p.adjust)
data<-data[1:50,]
pathway_level=as.factor(unique(data$Description))

p<-ggplot(data,aes(x=factor ,y=Description,size=log10padj,fill=NES))+
  geom_point(shape=21,stroke=0.04)+
  #scale_size_continuous(range = c(3,7))+
  scale_y_discrete(limits=rev(pathway_level))+
  scale_fill_gradient(low="#62B197",high="#E18E6D")+
  ylab("Pathways")+
  xlab("Factor")+
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


gsea_mRNA<-gsea_mRNA %>%
  mutate(view="cf-mRNA")
gsea_KO<-gsea_KO %>%
  mutate(view="microbiome KO")

gsea_mofa<-rbind(gsea_mRNA,gsea_KO)

write.xlsx(gsea_mofa,file = "GSEA_MOFA_260311.xlsx")

##### gsea by tutorial #####
pathways_df<-pathways %>% 
    mutate(tmp = 1) %>%                    # 先打标记
    pivot_wider(
        id_cols     = term,                    # 行名
        names_from  = value,                    # 列名
        values_from = tmp,
        values_fill = 0,                       # 缺失填 0
        values_fn   = length                   # 防重复：计数→1
    ) %>% column_to_rownames("term") %>% as.matrix()

res_kegg_mRNA_tutorial <- run_enrichment(model, 
                               feature.sets = pathways_df, 
                               view = "mRNA",
                               sign = "positive"
)


##### top 10 features enrichment #####
library(clusterProfiler)
library(org.Hs.eg.db)
load("__PATH_011__")
pathways <- do.call(rbind, lapply(seq_along(pathways$transciptome_ENTREZID), function(i) {
  data.frame(term=rep(names(pathways$transciptome_ENTREZID[i]),length(pathways$transciptome_ENTREZID[[i]])), value=pathways$transciptome_ENTREZID[[i]], stringsAsFactors=FALSE)
})) %>%
  dplyr::filter(grepl("(KEGG)",term)) %>%
  mutate(term=gsub("\\(KEGG\\) ","",term))

symbol_ensg<-read.table("__PATH_012__") %>%
  dplyr::rename(feature=V2,ENSG=V1)

pathgmt<-read.table("microbiome ko gsea复现/ko00001.forR.KO_only.gmt")

#res_gobp_mRNA<-data.frame()
res_kegg_mRNA<-data.frame()
res_kegg_KO<-data.frame()
weight_mi_lncRNA_flt<-data.frame()
for (i in 1:nrow(model@cache$variance_explained$r2_per_factor[[1]])) {   ##factor num
  weights <- get_weights(model, views = "all", factors = "all",as.data.frame = T,abs = T) %>%
    dplyr::filter(factor==paste0("Factor",i)) %>%
    arrange(desc(value))
  
  tmp<-weights[1:floor(nrow(weights)/10),] %>%
    dplyr::filter(view=="mRNA") %>%
    left_join(symbol_ensg,by="feature") 
  tmp<-tmp$ENSG
  
  set.seed(3490)
  #tmp_mRNA_gobp<-enrichGO(tmp,OrgDb = org.Hs.eg.db,keyType = "ENSEMBL",ont = "BP",pvalueCutoff = 1,pAdjustMethod = "BH")
  
  tmp_mRNA_kegg<-enricher(tmp,TERM2GENE = pathways,pvalueCutoff = 1,pAdjustMethod = "BH",minGSSize = 8)
  
  tmp<-weights[1:floor(nrow(weights)/10),] %>%
    dplyr::filter(view=="KO")
  tmp_KO<-enricher(tmp$feature,TERM2GENE = pathgmt,pvalueCutoff = 1,pAdjustMethod = "BH",minGSSize = 8)
  
  if (!is.null(tmp_KO)) {
    tmp_KO<-as.data.frame(tmp_KO@result) %>%
      mutate(factor=i)
    res_kegg_KO<-rbind(res_kegg_KO,tmp_KO)
  }
  
  #if (!is.null(tmp_mRNA_gobp)) {
  #  tmp_mRNA_gobp<-as.data.frame(tmp_mRNA_gobp@result) %>%
  #    mutate(factor=i)
  #  res_gobp_mRNA<-rbind(res_gobp_mRNA,tmp_mRNA_gobp)
  #}

  if (!is.null(tmp_mRNA_kegg)) {
    tmp_mRNA_kegg<-as.data.frame(tmp_mRNA_kegg@result) %>%
      mutate(factor=i)
    res_kegg_mRNA<-rbind(res_kegg_mRNA,tmp_mRNA_kegg)
  }
  
  tmp<-weights[1:floor(nrow(weights)),] 
  tmp_miRNA<-tmp %>%
    dplyr::filter(view=="miRNA") %>%
    mutate(factor=i)
  tmp_lncRNA<-tmp %>%
    dplyr::filter(view=="lncRNA") %>%
    mutate(factor=i)
  weight_mi_lncRNA_flt<-rbind(weight_mi_lncRNA_flt,tmp_miRNA,tmp_lncRNA)
}

res_kegg_KO<-res_kegg_KO %>%
  mutate(ID = gsub("ko","KO",
                   gsub("]","",
                        gsub(".*PATH:","",ID))),
         Description = gsub("_", " ",
                            gsub("^[0-9]*_","",
                                 gsub("_\\[PATH:.*","",Description))),
         view="microbiome KO")

res_kegg_mRNA<-res_kegg_mRNA %>%
  mutate(view="cf-mRNA")
  

tmp<-res_kegg_KO %>%
  dplyr::filter(p.adjust<0.01,factor=="3") %>%
  separate(GeneRatio, into = c("a", "b"), sep = "/", convert = TRUE) %>% # 分割为 "a" 和 "b"，并转为数值型
  mutate(GeneRatio = a / b,
         log10padj= -log10(p.adjust)) %>%
  dplyr::select(-a,-b)

cairo_pdf(filename = "__PATH_013__",
          width = 7,height = 5,family = "Arial")
ggplot(tmp,aes(x=GeneRatio,y=Description,fill=log10padj)) +
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_gradient(low="#795B3466", high="#795B34",name=TeX('$-log_{10}(p.adjust)$'))+
  ylab("KEGG Pathway")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank()          # 移除轴线（避免重复）
  )
dev.off()

ggsave(p,filename = "__PATH_013__",width = 9,height = 6)

tmp<-res_kegg_mRNA %>%
  dplyr::filter(p.adjust<0.05,factor=="2") %>%
  separate(GeneRatio, into = c("a", "b"), sep = "/", convert = TRUE) %>% # 分割为 "a" 和 "b"，并转为数值型
  mutate(GeneRatio = a / b) %>%
  dplyr::select(-a,-b)


##### target mRNA analysis of miRNA and lncRNA in factor #####
library("multiMiR")
library(org.Hs.eg.db)
library(biomaRt)
library(patchwork)

res_target_miRNA_kegg<-data.frame()
res_pred_mrna<-data.frame()
for (i in 1:2) {  
  intested_mirna<-weight_mi_lncRNA_flt %>%
    dplyr::filter(factor==i,view=="miRNA")
  intested_mirna<-intested_mirna[1:20,]
  pred_mrna<-get_multimir(org = "hsa",
                          mirna = intested_mirna$feature,
                          summary = T,
                          predicted.cutoff = 15,
                          predicted.cutoff.type = "p",
                          use.tibble = T)
  pred_mrna_res<-as.data.frame(pred_mrna@data)
  pred_mrna_res <- pred_mrna_res %>%
    group_by(mature_mirna_id, target_ensembl) %>%
    summarise(n_databases = n_distinct(database),
              databases = paste(unique(database), collapse = ";"),
              .groups = "drop") %>%
    filter(n_databases >= 2) %>%
    mutate(factor=i)
  
  res_pred_mrna<-rbind(res_pred_mrna,pred_mrna_res)
  
  tmp<-enricher(gene = pred_mrna_res$target_ensembl,pvalueCutoff = 0.05,pAdjustMethod = "BH",
                TERM2GENE = pathways)
  tmp<-tmp@result %>%
    mutate(factor=i)
  
  res_target_miRNA_kegg<-rbind(res_target_miRNA_kegg,tmp)
}

##plot
tmp<-res_target_miRNA_kegg %>%
  dplyr::filter(factor %in% c(1,2),
                p.adjust<0.05) %>%
  mutate(factor=gsub("1","Factor1",factor),
         factor=gsub("2","Factor2",factor),
         log10padj= -log10(p.adjust),
         tmp1=gsub(".*\\/","",GeneRatio),
         tmp2=gsub("\\/.*","",GeneRatio),
         GeneRatio=as.numeric(tmp2)/as.numeric(tmp1)) %>%
  arrange(log10padj)

p<-ggplot(tmp,aes(x=Description ,y=p.adjust,fill=GeneRatio))+
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_gradient(low="#C17F9E4D", high="#C17F9E")+
  coord_flip() +
  facet_wrap(~factor) +
  xlab("KEGG Pathway") +
  ylab("-log10(p.adjust)") +
  theme_classic()

ggsave(p,filename = "__PATH_014__",width = 9,height = 6)

###背靠背棒棒糖
p1<-ggplot(tmp,aes(x=Description ,y=p.adjust,fill=GeneRatio))+
  geom_bar(stat="identity",position = "dodge")+
  scale_fill_gradient(low="#C17F9E4D", high="#C17F9E")+
  coord_flip() +
  facet_wrap(~factor) +
  xlab("KEGG Pathway") +
  ylab("-log10(p.adjust)") +
  theme_classic()

ggsave(p,filename = "__PATH_014__",width = 9,height = 6)
save(res_target_miRNA_kegg,res_pred_mrna,file="res_pred_mRNA_noperiod_260311.xlsx")

###背靠背棒棒糖
pathway_level<-factor(unique(tmp$Description),levels = unique(tmp$Description))

p1<-ggplot(tmp[which(tmp$factor=="Factor1"),], aes(x = log10padj, y = Description)) +
  geom_col(fill ="#a52a2a", width = 0.1) +
  geom_point(aes(size = GeneRatio), color ="#a52a2a") +
  scale_size_continuous(range = c(0, 5)) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 2)) + #x轴刻度自定义
  scale_y_discrete(limits=pathway_level)+
  labs(title ='Factor 1') + #添加标题
  xlab("-log10(p.adj)")+
  theme_bw() +
  theme(
    plot.title = element_text(face ='bold', color ='#a42929', size = 10),##修改标题字体、字号、颜色
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.2,size = 0.3),
    axis.text.y = element_text(size = 8,colour ="black",hjust = 0.5,face = "bold"),
    axis.text.x = element_text(size = 6,face ="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )



p2<-ggplot(tmp[which(tmp$factor=="Factor2"),], aes(x = log10padj, y = Description)) +
  geom_col(fill ="#dc968d", width = 0.1) +
  geom_point(aes(size = GeneRatio), color ="#dc968d") +
  scale_size_continuous(range = c(0, 5)) +
  scale_x_reverse(limits = c(4, 0), breaks = seq(0, 4, 2)) + #x轴刻度自定义
  scale_y_discrete(limits=pathway_level,position = "right")+
  labs(title ='Factor 2') + #添加标题
  xlab("-log10(p.adj)")+
  theme_bw() +
  theme(
    plot.title = element_text(face ='bold', color ='#dc968d', size = 10),##修改标题字体、字号、颜色
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.2,size = 0.3),
    #axis.text.y = element_text(size = 7,colour ="black",hjust = 0.5),
    axis.text.x = element_text(size = 6,face ="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p <- (p2 | p1) + plot_layout(guides ='collect') &
  theme( legend.justification = c("right","bottom"),
         legend.text = element_text(face ="bold",size=8)
  )

cairo_pdf(filename = "MOFA_pred_miRNA_mRNA_kegg_enrich_factor1_2_260324.pdf",width = 8,height = 6,family = "Arial")
p
dev.off()


##### test #####
model<-load_model("__PATH_015__")
#load("__PATH_003__")
load("mofa_input_ko_m_mi_lncRNA_topVarience50_260120.rda")
load("__PATH_003__")
coldata_mofa<-coldata %>%
  mutate(sample=rownames(.))
#Nsamples = sum(get_dimensions(model)[["N"]])
#sample_metadata <- data.frame(sample = samples_names(model)[[1]],replace = TRUE) %>%
#  left_join(coldata_mofa[,c("sample","disease","period")])
samples_metadata(model) <- coldata_mofa %>%
  dplyr::select(sample,disease,period) %>%
  mutate(group=period)

head(get_variance_explained(model)$r2_total[[1]])
get_variance_explained(model)$r2_per_factor[[1]]


coldata_wgcna_mofa<-coldata[,c(1,4,14:61,65:85,87:129,132,134,138:191,195:199,201:205)]

factor <- data.frame()
for (i in 1:3) {
  tmp<-as.data.frame(get_factors(model, factors = "all",groups = i,
                                 as.data.frame = F, scale = T)) %>%
    dplyr::rename_with(~gsub("^X..F","F",.x))
  factor<-rbind(factor,tmp)
}

coldata_wgcna_mofa<-coldata_wgcna_mofa[rownames(factor),] %>%
  rownames_to_column() %>%
  dplyr::rename(sample=rowname)

factor<-factor %>%
  rownames_to_column() %>%
  dplyr::rename(sample=rowname) %>%
  left_join(coldata_wgcna_mofa) %>%
  mutate(disease=gsub("0","control",disease),
         disease=gsub("1","GDM",disease)) %>%
  pivot_longer(cols = names(.)[2:11],names_to = "factor",values_to = "factor_value") %>%
  mutate(factor=factor(factor,levels=c("Factor1","Factor2","Factor3","Factor4","Factor5",
                                          "Factor6","Factor7","Factor8","Factor9","Factor10")))

ggplot(factor,aes(x=disease,y=factor_value)) +
  geom_boxplot() +
  ggsignif::stat_signif(comparisons = list(c("GDM","control")),
                        test = "wilcox.test",
                        map_signif_level = F,   # 自动 */**/***
                        y_position  = 0.5,
                        textsize = 3,
                        color="grey50")+
  facet_wrap(~factor)
  



