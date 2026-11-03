library(openxlsx)
library(readxl)
library(tidyverse)
library(sjmisc)
library(ggord)
library(FactoMineR)

rm(list=ls())
load("__PATH_001__")
load("select_cfRNA_251205.rda")

interact_gene<-Reduce(intersect, list(rownames(limma_voom_multiRNA), rownames(limma_voom_multiRNA_lyw),rownames(limma_voom_NCBI)))

test_tianjin<-limma_voom_multiRNA %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% interact_gene)
test_lyw<-limma_voom_multiRNA_lyw %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% interact_gene)
test_NCBI<-limma_voom_NCBI %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% interact_gene)

test_all<-cbind(test_tianjin,test_lyw,test_NCBI)

group_list_group<-as.character(c(rep("Tianjin",ncol(test_tianjin)),
                                 rep("LYW",ncol(test_lyw)),
                                 rep("NCBI",ncol(test_NCBI))))

dat <- as.data.frame(t(test_all))
dat.pca <- PCA(dat, graph = FALSE)
col_cohort<-c("Tianjin"="#AED185","LYW"="#F89FA8","NCBI"="#AECDD7")
test<-data.frame(sample=rownames(dat.pca$ind$coord),pc1=dat.pca$ind$coord[,1],pc2=dat.pca$ind$coord[,2],
                 cohort=group_list_group)

p_cohort<-ggplot(test,aes(x=pc1,y=pc2,color=cohort,shape=cohort))+
  geom_point(size=2,alpha=0.8)+
  stat_ellipse(aes(color=cohort),type="norm",geom="polygon",fill=NA)+
  scale_color_manual(values = col_cohort)+
  scale_shape_manual(values = c("Tianjin"=16,"LYW"=18,"NCBI"=17))+
  #scale_color_manual(values = c(NA,NA))+
  xlab(paste0("PC1 ",round(dat.pca$ind$coord[,1],2)," (%)"))+
  ylab(paste0("PC2 ",round(dat.pca$ind$coord[,2],2)," (%)"))+
  ggtitle("mRNA")+
  labs(fill = "cohort")+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = 'gray'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_text(color="black",face="bold"),
        axis.title = element_text(color = "black",face="bold"),
        plot.title = element_text(color="black",face="bold",size=15,hjust = 0.5),
        legend.text = element_text(face="bold")
  )


ggsave(p_cohort,filename = "投稿期刊/cell reports/cfRNA_PCA_motif_QC_3cohort.pdf",width = 10,height = 10)


## ------------------------------------------------------------------
## 新增: 三个数据集 removeBatch (以 cohort 为批次) + 去批次后重新 PCA
## ------------------------------------------------------------------
library(limma)
library(vegan)

## 合并为 基因 x 样本 (行=基因, 列=样本), 与 removeBatchEffect 输入要求一致
## (文档: rows=probes/features, columns=samples)
mat_all      <- cbind(test_tianjin, test_lyw, test_NCBI)
batch_cohort <- factor(group_list_group)

## removeBatchEffect: 去掉 cohort 批次效应, 保留生物学变异
## design 设为 ~1 (仅去除批次); 如需同时保留 GDM/Control 等分组可改为模型矩阵
design_batch <- model.matrix(~ 1, data = data.frame(cohort = batch_cohort))
mat_corrected <- removeBatchEffect(
  as.matrix(mat_all),          # 基因 x 样本 (行=基因, 列=样本)
  batch = batch_cohort,
  design = design_batch
)

## 拆回三个队列 (基因 x 样本, 按列切分), 供后续使用
n_tj <- ncol(test_tianjin); n_lyw <- ncol(test_lyw); n_ncbi <- ncol(test_NCBI)
test_tianjin_rb <- as.data.frame(mat_corrected[, 1:n_tj])
test_lyw_rb     <- as.data.frame(mat_corrected[, (n_tj+1):(n_tj+n_lyw)])
test_NCBI_rb    <- as.data.frame(mat_corrected[, (n_tj+n_lyw+1):(n_tj+n_lyw+n_ncbi)])

## 去批次后重新 PCA: 转置为 样本 x 基因 (行=样本) 再 PCA
dat.rb        <- as.data.frame(t(mat_corrected))
dat.rb.pca    <- PCA(dat.rb, graph = FALSE)

test.rb <- data.frame(
  sample = rownames(dat.rb.pca$ind$coord),
  pc1    = dat.rb.pca$ind$coord[, 1],
  pc2    = dat.rb.pca$ind$coord[, 2],
  cohort = batch_cohort
)

p_cohort_rb <- ggplot(test.rb, aes(x = pc1, y = pc2, color = cohort, shape = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(color = cohort), type = "norm", geom = "polygon", fill = NA) +
  scale_color_manual(values = col_cohort) +
  scale_shape_manual(values = c("Tianjin" = 16, "LYW" = 18, "NCBI" = 17)) +
  xlab(paste0("PC1 ", round(dat.rb.pca$ind$coord[, 1], 2), " (%)")) +
  ylab(paste0("PC2 ", round(dat.rb.pca$ind$coord[, 2], 2), " (%)")) +
  ggtitle("mRNA (batch-corrected)") +
  labs(fill = "cohort") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = 'gray'),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(color = "black", face = "bold", size = 15, hjust = 0.5),
    legend.text = element_text(face = "bold")
  )

ggsave(p_cohort_rb,
       filename = "投稿期刊/cell reports/cfRNA_PCA_motif_QC_3cohort_batchcorrected.pdf",
       width = 10, height = 10)

## 如需分别保存三个去批次后数据集, 取消下方注释
# save(test_tianjin_rb, test_lyw_rb, test_NCBI_rb, file = "interact_gene_3cohort_batchcorrected.rda")


## ------------------------------------------------------------------
## 新增: PERMANOVA 评估 cohort(batch) 解释度, 去批次前后对比
## ------------------------------------------------------------------
library(vegan)

## 距离矩阵: 表达量为 limma-voom 归一化 (近似正态) -> euclidean 合适
## 若为原始 count/相对丰度, 改 method = "bray"
## mat_all / mat_corrected 是 基因 x 样本, vegdist 按行算距离, 需转置为 样本 x 基因
dist_raw <- vegdist(as.matrix(t(mat_all)),   method = "euclidean")  # 去批次前
dist_rb  <- vegdist(as.matrix(t(mat_corrected)), method = "euclidean")  # 去批次后

## 去批次前: cohort 解释度 (应较高, 说明存在 batch 效应)
perma_raw <- adonis2(dist_raw ~ batch_cohort, permutations = 999)
cat("\n===== PERMANOVA (before batch correction) =====\n")
print(perma_raw)

## 去批次后: cohort 解释度 (应大幅下降, 说明 batch 被有效移除)
perma_rb <- adonis2(dist_rb ~ batch_cohort, permutations = 999)
cat("\n===== PERMANOVA (after batch correction) =====\n")
print(perma_rb)

## 汇总对比: 提取 cohort 行的 R2 与 p 值
perma_compare <- data.frame(
  stage    = c("before", "after"),
  R2_cohort = c(perma_raw$R2["batch_cohort"],  perma_rb$R2["batch_cohort"]),
  p_cohort  = c(perma_raw$`Pr(>F)`["batch_cohort"], perma_rb$`Pr(>F)`["batch_cohort"])
)
print(perma_compare)

## 如需保存对比结果
# write.csv(perma_compare, file = "PERMANOVA_batch_explained_3cohort.csv", row.names = FALSE)
