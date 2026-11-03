library(MOFA2)
library(tidyverse)
library(readxl)
library(reticulate)

rm(list = ls())
setwd("__PATH_001__")
use_condaenv("__PATH_002__",required = T)

#load("mofa2_tmp.rda")
#MOFAobject.trained <- run_mofa(MOFAobject, outfile = "tmp_model.hdf5", use_basilisk=F)

#####
load("input/coldata_wgcna_260105.rda")
#load("__PATH_003__")
sample_meta<-read_excel("input/tianjinsamples_240403.xlsx",sheet = 10) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_meta=sampleid)
sample_cfRNA<-read_excel("input/tianjinsamples_240403.xlsx",sheet = 8) %>%
  mutate(visit=gsub("孕早期","early",visit),visit=gsub("孕中期","middle",visit),visit=gsub("孕晚期","late",visit),
         sign=paste(sjid_mo,visit,sep = "")) %>%
  dplyr::rename(sampleid_cfRNA=sampleid) %>%
  left_join(sample_meta[,c("sign","sampleid_meta")])


##ko rpm
rpm_meta<-read.csv("input/TJGDM1962s_ko_tri_95_rpk.csv",header = T,row.names = 1)
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


##strain rpm
rpm_meta<-read.csv("input/TJGDM1962s_microbiome_ab_tri_95.csv",header = T,row.names = 1)
rpm_meta<-rpm_meta[apply(rpm_meta,1,function(x) sum(x>0) > floor(ncol(rpm_meta)/10)),]
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
#load("input/wgcna_vst_m_mi_lnc_260105.rda")

## cfRNA TPM
tpm_m_lncrna<-read.table("__PATH_004__",row.names = 1,header = T)
tpm_mirna<-read.table("__PATH_005__",header = T,row.names = 1)

all(names(tpm_m_lncrna)==names(tpm_mirna))
tpm_m_mi_lncrna<-rbind(tpm_m_lncrna,tpm_mirna)
tpm_m_mi_lncrna_flt<-tpm_m_mi_lncrna[colnames(vst_cfRNA),]
tmp<-sample_cfRNA %>%
  mutate(sampleid=gsub("^","S",sampleid_cfRNA)) %>%
  dplyr::filter(sampleid %in% names(count_mirna))
tpm_m_mi_lncrna_flt<-tpm_m_mi_lncrna_flt %>%
  dplyr::rename(!!!setNames(tmp$sampleid,tmp$sign))
tpm_m_mi_lncrna_flt<-tpm_m_mi_lncrna_flt[,rownames(vst_cfRNA)]
save(tpm_m_mi_lncrna_flt,file="__PATH_006__")

load("input/tpm_multiRNA_flt_wgcna_260112.rda")

## cfRNA first merge count matrix then vst
count_mrna<-read.table("__PATH_007__",row.names = 1,header = T)
count_lncrna<-read.table("__PATH_008__",row.names = 1,header = T)
count_mirna<-read.table("__PATH_009__",row.names = 1,header = T)

all(names(count_mrna)==names(count_mirna))
all(names(count_mrna)==names(count_lncrna))
count_m_mi_lncrna<-rbind(count_mrna,count_lncrna,count_mirna)





#vst_cfRNA<-as.data.frame(t(vst_cfRNA))
##long data frame
tmp_meta<-rpm_meta %>%
  rownames_to_column() %>%
  dplyr::rename(sample=rowname) %>%
  pivot_longer(cols = !sample,names_to = "feature",values_to = "value") %>%
  mutate(view="KO",
         value=as.numeric(value)) %>%
  dplyr::select(sample,feature,view,value)
#tmp_rna<-vst_cfRNA %>%
#  as.data.frame() %>%
#  rownames_to_column() %>%
#  dplyr::rename(sample=rowname) %>%
#  pivot_longer(cols = !sample,names_to = "feature",values_to = "value") %>%
#  mutate(view="rna",
#         value=as.numeric(value)) %>%
#  dplyr::select(sample,feature,view,value)

tmp_rna<-tpm_m_mi_lncrna_flt %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(sample=rowname) %>%
  pivot_longer(cols = !sample,names_to = "feature",values_to = "value") %>%
  mutate(view="rna",
         value=as.numeric(value)) %>%
  dplyr::select(sample,feature,view,value)

MOFAobject<-rbind(tmp_rna,tmp_meta)

##add groups
MOFAobject<-list(KO=t(rpm_meta),
                 mRNA=t(vst_mrna),
                 miRNA=t(vst_mirna),
                 lncRNA=t(vst_lncrna)
)

MOFAobject<-list(KO=t(rpm_meta),
                 mRNA=(tpm_mRNA_combat_log),
                 miRNA=(tpm_miRNA_combat_log),
                 lncRNA=(tpm_lncRNA_combat_log)
)



tmp<-coldata_wgcna %>%
  mutate(sample=rownames(.))
tmp<-tmp[colnames(MOFAobject$mRNA),]
MOFAobject <- create_mofa(MOFAobject, groups=as.character(tmp$period))


##no group
MOFAobject <- create_mofa(MOFAobject)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views<-T
train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose <- T

train_opts$convergence_mode <- "slow"

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
MOFAobject.trained <- run_mofa(MOFAobject, outfile = "model_tpm_scale.hdf5", use_basilisk=F)

head(get_variance_explained(model)$r2_total[[1]])
get_variance_explained(model)$r2_per_factor[[1]]



##### downstream analysis #####
coldata_wgcna<-coldata_wgcna %>%
  mutate(sample=rownames(.))

Nsamples = sum(get_dimensions(model)[["N"]])

sample_metadata <- data.frame(sample = samples_names(model)[[1]],replace = TRUE) %>%
  left_join(coldata_wgcna[,c("sample","disease","period")])

samples_metadata(model) <- sample_metadata


## plot
p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "disease",
                 dot_size = 3,        # change dot size
                 dodge = TRUE,           # dodge points with different colors
                 legend = FALSE,          # remove legend
                 add_violin = TRUE,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("1"="black", "0"="red")) +
  scale_fill_manual(values=c("1"="black", "0"="red"))


##scatter plot
plot_factors(model, 
             factors = c(1,2,3),
             color_by = "disease"
)

