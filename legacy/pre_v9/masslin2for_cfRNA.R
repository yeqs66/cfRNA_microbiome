rm(list = ls())
library(tidyverse)
library(ggplot2)
library(vegan)
library(sjmisc)
library(ggpubr)
library(rstatix)
library(sjmisc)
library(MetBrewer)
library(Maaslin2)


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/")
load("metabolics_matrix_nottranform.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/trajGWAS/microbiome_sampleinfo.rda")
load("final_files/data347final_confirmed.rda")
TJGDM1962s_microbiome_ab_tri_95 <- read.csv("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/trajGWAS/TJGDM1962s_microbiome_ab_tri_95.csv")




abundance <- TJGDM1962s_microbiome_ab_tri_95 %>%
  column_to_rownames("X") %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  inner_join(microbiome_sampleinfo, by = "sample") %>%
  left_join(data347newnewnew20240627,by ="ID")


inputmatrix <- abundance %>%
  dplyr::select(sample,starts_with("s__")) %>%
  column_to_rownames("sample")



inputmetainfo <- abundance %>%
  column_to_rownames("sample")


fit_adjust = Maaslin2(input_data  = inputmatrix, 
                           input_metadata = inputmetainfo, 
                           output         = "masslin2result/compare", 
                           fixed_effects  = c("disease","trimester","smoking_mo","drinking_mo_18topreg","ethnicity_mo","BMI_mo"),
                           random_effects = c('ID'),
                           reference      = c("disease,Control",
                                              "trimester,early",
                                              "smoking_mo,No",
                                              "drinking_mo_18topreg,No",
                                              "ethnicity_mo,Han"),
                           min_abundance = 0.05,
                           min_prevalence = 0.05,
                           correction = "BH",
                           plot_scatter = FALSE)



fit_unadjust = Maaslin2(input_data  = inputmatrix, 
                           input_metadata = inputmetainfo, 
                           output         = "masslin2result/unadjustcompare", 
                           fixed_effects  = c("disease","trimester"),
                           random_effects = c('ID'),
                           reference      = c("disease,Control",
                                              "trimester,early"),
                           min_abundance = 0.05,
                           min_prevalence = 0.05,
                           correction = "BH",
                           plot_scatter = FALSE)