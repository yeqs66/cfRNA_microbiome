library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(lmerTest)
library(lme4)
library(lsmeans)
library(ggbreak)



coldata_rna<-read.csv("__PATH_001__",row.names = 1) %>%
  mutate(ID=as.character(ID),
          trimester=case_when(
            trimester=="1"~"early",
            trimester=="2"~"middle",
            trimester=="3"~"late"
          ),
        trimester=factor(trimester,levels=c("early","middle","late"))) %>%
  dplyr::select(sample,ID,trimester,disease,batch)
res_deconformer_preg<-read.table("细胞溯源/res_deconv_349.txt") %>%
  dplyr::rename_with(~paste0("preg_",.x)) %>%
  mutate(sample=rownames(.)) %>%
  left_join(coldata_rna[,c("sample","trimester","disease","batch","ID")])

res_deconformer_ori<-read.table("__PATH_002__")%>%
  dplyr::rename_with(~paste0("ori_",.x)) %>%
  mutate(sample=rownames(.))

m <- inner_join(res_deconformer_ori, res_deconformer_preg, by = "sample") %>%
  mutate(trimester=factor(trimester,levels=c("early","middle","late")), 
          disease=factor(disease,levels=c("GDM","control")), 
          batch=as.character(batch)) %>%
  column_to_rownames("sample")

##### 合并数据 #####
meta_cols <- c("sample", "trimester", "disease", "ID", "batch")
cell_a <- setdiff(names(res_deconformer_ori), meta_cols)
cell_p <- setdiff(names(res_deconformer_preg), meta_cols)
shared <- intersect(sub("^ori_", "", cell_a), sub("^preg_", "", cell_p))   # 共有成人类型
tropho <- setdiff(sub("^preg_", "", cell_p), shared)      

# ---- Step 1: 比例一致性（Spearman ρ + Lin's CCC）----
ccc <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  (2 * cov(x, y)) / (var(x) + var(y) + (mx - my)^2)
}

frac_corr <- map_dfr(shared, function(ct) {
  xa <- m[[paste0("ori_", ct)]]; xp <- m[[paste0("preg_", ct)]]
  sp <- cor.test(xa, xp, method = "spearman")
  tibble(cell_type = ct,
         spearman_r = sp$estimate, spearman_p = sp$p.value,
         ccc = ccc(xa, xp))
})


##### remove preg cell #####
renorm <- function(df, cols) {
  tot <- rowSums(df[cols], na.rm = TRUE)
  df[cols] <- sweep(as.matrix(df[cols]), 1, tot, `/`, 1)
  df
}

preg_cell<-paste0("preg_", tropho)

res_deconformer_preg_np<-res_deconformer_preg
res_deconformer_preg_np[preg_cell]<-0
res_deconformer_preg_np<-renorm(res_deconformer_preg_np, setdiff(names(res_deconformer_preg_np),  meta_cols)) %>%
  dplyr::select(-preg_cell)

m_np<-inner_join(res_deconformer_ori, res_deconformer_preg_np, by = "sample") %>%
  mutate(trimester=factor(trimester,levels=c("early","middle","late")), 
          disease=factor(disease,levels=c("GDM","control")), 
          batch=as.character(batch)) %>%
  column_to_rownames("sample")

ccc <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  (2 * cov(x, y)) / (var(x) + var(y) + (mx - my)^2)
}

frac_corr_np <- map_dfr(shared, function(ct) {
  xa <- m_np[[paste0("ori_", ct)]]
  xp <- m_np[[paste0("preg_", ct)]]
  sp <- cor.test(xa, xp, method = "spearman")
  tibble(cell_type = ct,
         spearman_r = sp$estimate, spearman_p = sp$p.value,
         pct_ori  = mean(m_np[[paste0("ori_",  ct)]], na.rm = TRUE),
          pct_preg = mean(m_np[[paste0("preg_", ct)]], na.rm = TRUE),
         ccc = ccc(xa, xp))
})




##### remove blood+preg cell #####
blood_cells <- c(
  "innate_lymphoid_cell",
  "neutrophil",
  "basophil",
  "b_cell",
  "monocyte",
  "plasmablast",
  "plasma_cell",
  "erythrocyte__erythroid_progenitor",
  "mature_conventional_dendritic_cell",
  "myeloid_progenitor",
  "mast_cell",
  "platelet",
  "macrophage",
  "hematopoietic_stem_cell",
  "t_cell",
  "thymocyte",
  "nk_cell"
)



blood_ori<-paste0("ori_", blood_cells)
blood_preg<-paste0("preg_", blood_cells)
preg_cell<-paste0("preg_", tropho)

res_deconformer_ori_nb<-res_deconformer_ori
res_deconformer_ori_nb[blood_ori]<-0
res_deconformer_ori_nb<-renorm(res_deconformer_ori_nb, setdiff(names(res_deconformer_ori_nb),  meta_cols)) %>%
  dplyr::select(-blood_ori)

res_deconformer_preg_nbnp<-res_deconformer_preg
res_deconformer_preg_nbnp[c(blood_preg,preg_cell)]<-0
res_deconformer_preg_nbnp<-renorm(res_deconformer_preg_nbnp, setdiff(names(res_deconformer_preg_nbnp),  meta_cols)) %>%
  dplyr::select(-c(blood_preg,preg_cell))

m_nbnp<-inner_join(res_deconformer_ori_nb, res_deconformer_preg_nbnp, by = "sample") %>%
  mutate(trimester=factor(trimester,levels=c("early","middle","late")), 
          disease=factor(disease,levels=c("GDM","control")), 
          batch=as.character(batch)) %>%
  column_to_rownames("sample")

ccc <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  (2 * cov(x, y)) / (var(x) + var(y) + (mx - my)^2)
}

frac_corr_nbnp <- map_dfr(setdiff(shared, blood_cells), function(ct) {
  xa <- m_nbnp[[paste0("ori_", ct)]]
  xp <- m_nbnp[[paste0("preg_", ct)]]
  sp <- cor.test(xa, xp, method = "spearman")
  tibble(cell_type = ct,
         spearman_r = sp$estimate, spearman_p = sp$p.value,
         ccc = ccc(xa, xp))
})



##### GDM effect #####
blood_ori<-paste0("ori_", blood_cells)
blood_preg<-paste0("preg_", blood_cells)

res_deconformer_ori_nb<-res_deconformer_ori
res_deconformer_ori_nb[blood_ori]<-0
res_deconformer_ori_nb<-renorm(res_deconformer_ori_nb, setdiff(names(res_deconformer_ori_nb),  meta_cols)) %>%
  dplyr::select(-blood_ori)

res_deconformer_preg_nb<-res_deconformer_preg
res_deconformer_preg_nb[blood_preg]<-0
res_deconformer_preg_nb<-renorm(res_deconformer_preg_nb, setdiff(names(res_deconformer_preg_nb),  meta_cols)) %>%
  dplyr::select(-blood_preg)

m_nb<-inner_join(res_deconformer_ori_nb, res_deconformer_preg_nb, by = "sample") %>%
  mutate(trimester=factor(trimester,levels=c("early","middle","late")), 
          disease=factor(disease,levels=c("GDM","control")), 
          batch=as.character(batch)) %>%
  column_to_rownames("sample")

## gdm_effect：在指定子集上拟合 disease 效应（不区分孕周的子集用全部数据，
## 区分孕周的子集只取对应 trimester，公式中去掉 trimester 项）
gdm_effect <- function(col, data, trimester_subset = NULL) {
  sub <- if (is.null(trimester_subset)) {
    data
  } else {
    data %>% filter(trimester == trimester_subset)
  }
  if (nrow(sub) < 5 || var(sub[[col]], na.rm = TRUE) %in% c(0, NA)) return(NA_real_)
  form <- if (is.null(trimester_subset)) {
    as.formula(paste0(col, " ~ disease + trimester"))
  } else {
    as.formula(paste0(col, " ~ disease"))
  }
  fit <- lm(form, data = sub)
  emresult <- tryCatch(
    emmeans::emmeans(fit, specs = pairwise ~ disease),
    error = function(e) {
      warning("emmeans failed for '", col, "' (",
              ifelse(is.null(trimester_subset), "all", trimester_subset),
              "): ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(emresult)) return(NA_real_)
  contrast_tbl <- as.data.frame(emresult[["contrasts"]])
  contrast_tbl[["estimate"]]
}

## 4 种建模整合：early / middle / late / 不区分 (all)
trimester_levels <- c("early", "middle", "late", "all")
eff_all <- map_dfr(trimester_levels, function(tr) {
  ts <- if (tr == "all") NULL else tr
  map_dfr(setdiff(shared, blood_cells), function(ct) {
    ba <- gdm_effect(paste0("ori_", ct),  data = m_nb, trimester_subset = ts)
    bp <- gdm_effect(paste0("preg_", ct), data = m_nb, trimester_subset = ts)
    tibble(trimester = tr,
           cell_type = ct, beta_adult = ba, beta_preg = bp,
           pct_ori  = mean(m_nb[[paste0("ori_",  ct)]], na.rm = TRUE),
           pct_preg = mean(m_nb[[paste0("preg_", ct)]], na.rm = TRUE),
           direction_agree = sign(ba) == sign(bp))
  })
})

write.csv("投稿期刊/cell reports/细胞溯源换ref敏感性分析_GDM_effect.csv",row.names=F)


## 汇总指标（按 4 种 model 分别）
eff_summary <- eff_all %>%
  group_by(model) %>%
  summarise(
    eff_mag_rho  = cor.test(beta_adult, beta_preg, method = "spearman")$estimate,
    sign_concord = mean(direction_agree, na.rm = TRUE),
    rank_rho     = cor.test(rank(-abs(beta_adult)), rank(-abs(beta_preg)),
                            method = "spearman")$estimate
  )


##### tropho pct in trimester #####
preg_cell_new<-data.frame(oir_cell=c("early_VCT","late_VCT","early_SCT","late_SCT","early_EVT","late_EVT"),
                      new_cell=c("The proportion of Early Villous cytotrophoblast (VCT)","The proportion of Late Villous cytotrophoblast (VCT)",
                      "The proportion of Early Syncytiotrophoblast (SCT)","The proportion of Late Syncytiotrophoblast (SCT)",
                      "The proportion of Early Extra-villous trophoblast (EVT)","The proportion of Late Extra-villous trophoblast (EVT)"))

preg_cell_rename<-gsub("preg_","",preg_cell)

preg_cell_new <- preg_cell_new[match(preg_cell_rename, preg_cell_new$oir_cell), ]

p_list<-list()
for (i in sub("^preg_", "", preg_cell)) {
  data<-res_deconformer_preg_nb %>%
    dplyr::select(sample, trimester, disease, batch, ID, paste0("preg_", i)) %>%
    rename(cell_type= paste0("preg_", i)) 
  ymean <- mean(data$cell_type)
  step <- 0.005
  p<-ggplot(data, aes(x = trimester, y = cell_type, fill = trimester)) +
    geom_boxplot(outliers = F,outlier.shape = NA) +
    #geom_jitter(width = 0.2, alpha = 0.5) +
    #coord_cartesian(ylim = c(0, 10)) +
    #ylim(0, 10) +
    #coord_cartesian(ylim = c(0, 10)) +
    #scale_y_break(breaks = c(10,60)) +
    labs(x = "Trimester", y = "Cell Type Proportion") +
    ggtitle(paste("The proportion of ",preg_cell_new$new_cell[match(i, preg_cell_rename)],sep = ""))+
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    #stat_compare_means(method = "kruskal.test") +
    #stat_compare_means(comparisons = list(c("early", "middle"), c("middle", "late"), c("early", "late")),
    #                  method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
    geom_signif(comparisons = list(c("early", "middle"), c("middle", "late"), c("early", "late")),
                  test = "wilcox.test",
                  map_signif_level = F,   # 自动 */**/***
                  y_position = c(ymean + step, ymean + 2*step, ymean + 3*step),
                  tip_length = 0.00005,
                  textsize = 3,
                  color="grey50") +
    theme(legend.position = "none")
  p_list[[i]]<-p
  #ggsave(p, filename = paste0("细胞溯源/add_preg/tropho_pct_", i, ".pdf"), width = 6, height = 4)
}

p_list[["early_VCT"]]<-

( p_list[[1]] + p_list[[2]] + p_list[[3]] ) / ( p_list[[4]] + p_list[[5]] + p_list[[6]]) +
  plot_annotation(title = "Trophoblast cell proportions by trimester") &
  theme(plot.title = element_text(hjust = 0.5,size=9))

ggsave("细胞溯源/add_preg/tropho_all_3x2.pdf", width = 12, height = 8)
