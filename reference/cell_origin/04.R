library(tidyverse)

preg_cells <- c(
  "early SCT",
  "late SCT",
  "early EVT",
  "late EVT",
  "early VCT",
  "late VCT"
)

res_deconformer_ori<-read.table("__PATH_001__") %>%
  mutate(sample=rownames(.)) %>%
  pivot_longer(cols = !sample,names_to = "cell_type",values_to = "score_ori") %>%
  mutate(cell_type=gsub("_"," ",cell_type),cell_type=gsub(" "," ",cell_type),cell_type=gsub("  "," ",cell_type),
         sign=paste0(sample,"_",cell_type))

res_deconformer_preg<-read.table("细胞溯源/res_deconv_349.txt") %>%
  mutate(sample=rownames(.)) %>%
  pivot_longer(cols = !sample,names_to = "cell_type",values_to = "score") %>%
  mutate(is_preg= cell_type %in% preg_cells) %>%
  group_by(sample) %>%
  mutate(sum_nopreg=sum(score[!is_preg])) %>%
  ungroup() %>%
  filter(!is_preg) %>%                    # 如果你不想保留血细胞，直接去掉
  mutate(score_nopreg = score / sum_nopreg,
         cell_type=gsub("_"," ",cell_type),cell_type=gsub(" "," ",cell_type),cell_type=gsub("  "," ",cell_type),
        sign=paste0(sample,"_",cell_type)) %>%
  left_join(res_deconformer_ori[,c("sign","score_ori")],by="sign")

## test
res_cor<-data.frame()
for (i in unique(res_deconformer_preg$cell_type)) {
  if (i %in% preg_cells) next
  tmp<-res_deconformer_preg %>% filter(cell_type==i,!cell_type %in% preg_cells)
  tmp<-cor.test(tmp$score_nopreg,tmp$score_ori,method = "spearman")
  res_cor<-rbind(res_cor,data.frame(cell_type=i,r=tmp$estimate,p=tmp$p.value))
}
