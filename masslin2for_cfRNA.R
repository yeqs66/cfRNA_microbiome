rm(list = ls())

# ============================================================
# Maaslin2 Association Analysis: Gut Microbiome and GDM
# ============================================================
# This script uses the Maaslin2 (Multivariable Association with
# Linear Models) framework to identify gut microbial species
# (at species level, "s__" prefix) that are significantly
# associated with gestational diabetes mellitus (GDM) across
# pregnancy trimesters.
#
# Two models are fitted:
#   1. Adjusted model  – includes multiple covariates
#   2. Unadjusted model – disease status and trimester only
#
# Maaslin2 applies generalised linear models with optional
# random effects, filters on minimum abundance/prevalence,
# and corrects for multiple testing (Benjamini-Hochberg by default).
# ============================================================

library(tidyverse)
library(ggplot2)
library(vegan)
library(sjmisc)
library(ggpubr)
library(rstatix)
library(MetBrewer)
library(Maaslin2)

# ============================================================
# NOTE: Update the paths below to match your local data location
# before running this script.
# ============================================================
setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/")
load("metabolics_matrix_nottranform.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/trajGWAS/microbiome_sampleinfo.rda")
load("final_files/data347final_confirmed.rda")
TJGDM1962s_microbiome_ab_tri_95 <- read.csv("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/trajGWAS/TJGDM1962s_microbiome_ab_tri_95.csv")


# ============================================================
# Data preparation
# ============================================================

# Combine species-level microbiome abundance table with sample
# metadata (trimester, disease status, covariates) by sample ID
abundance <- TJGDM1962s_microbiome_ab_tri_95 %>%
  column_to_rownames("X") %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  inner_join(microbiome_sampleinfo, by = "sample") %>%
  left_join(data347newnewnew20240627,by ="ID")

# Feature matrix: rows = samples, columns = species (s__ prefix)
inputmatrix <- abundance %>%
  dplyr::select(sample,starts_with("s__")) %>%
  column_to_rownames("sample")

# Metadata matrix: rows = samples, columns = all covariates
inputmetainfo <- abundance %>%
  column_to_rownames("sample")


# ============================================================
# Model 1: Adjusted association analysis
# Fixed effects: disease, trimester, and confounders
# Random effect: subject ID (repeated measures across trimesters)
# ============================================================
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
                           min_abundance = 0.05,   # exclude species with mean relative abundance < 5%
                           min_prevalence = 0.05,  # exclude species present in < 5% of samples
                           correction = "BH",
                           plot_scatter = FALSE)



# ============================================================
# Model 2: Unadjusted association analysis
# Fixed effects: disease status and trimester only
# ============================================================
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
