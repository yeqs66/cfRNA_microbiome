library(tidyverse)
library(readxl)

files<-dir("BMI_sensitive/")

res_mRNA<-data.frame()
for (i in files[grep("mRNA.*csv",files)]){
  group<-gsub(".*_","",i)
  group<-gsub(".csv","",group)
  tmp<-read.csv(paste0("BMI_sensitive/", i)) %>%
    mutate(group=group,RNA_type="mRNA")
  res_mRNA<-rbind(res_mRNA, tmp)
}
res_mRNA<-res_mRNA %>%
  mutate(period=case_when(period=="1"~"early",
                  period=="2"~"middle",
                  period=="3"~"late",
                  period=="all"~"all_period"))

res_mRNA_ori<-read.csv("limma_349_batch_only.csv") %>%
    mutate(RNA_type="mRNA") %>%
  mutate(period=case_when(period=="1"~"early",
                  period=="2"~"middle",
                  period=="3"~"late",
                  period=="all"~"all_period"))

res_miRNA<-data.frame()
for (i in files[grep("miRNA.*csv",files)]){
  group<-gsub(".*_","",i)
  group<-gsub(".csv","",group)
  tmp<-read.csv(paste0("BMI_sensitive/", i)) %>%
    mutate(group=group,RNA_type="miRNA")
  res_miRNA<-rbind(res_miRNA, tmp)
}

res_lncRNA<-data.frame()
for (i in files[grep("lncRNA.*csv",files)]){
  group<-gsub(".*_","",i)
  group<-gsub(".csv","",group)
  tmp<-read.csv(paste0("BMI_sensitive/", i)) %>%
    mutate(group=group,RNA_type="lncRNA")
  res_lncRNA<-rbind(res_lncRNA, tmp)
}

res_allRNA<-rbind(res_mRNA,res_miRNA,res_lncRNA)



test_sensitive<-function(df_ori,df_sensitive){
  ## jaccard with ori
  for (tri in unique(df_sensitive$period)){
    df_ori_tri<-df_ori %>% filter(period==tri)
    df_sensitive_tri<-df_sensitive %>% filter(period==tri)
    sig_ori_tri<-df_ori_tri$SYMBOL[df_ori_tri$"adj.P.Val"<0.1]
    sig_sensitive_tri<-df_sensitive_tri$SYMBOL[df_sensitive_tri$"adj.P.Val"<0.1]
    sig_common<-intersect(sig_ori_tri,sig_sensitive_tri)
    sig_share<-union(sig_ori_tri,sig_sensitive_tri)
    jaccard<-length(sig_common)/length(sig_share)
    print(paste0("Jaccard with ori ",tri,":",jaccard))
  }
  
  ## jaccard within sensitive
  sig_sensitive_early<-df_sensitive$SYMBOL[which(df_sensitive$"adj.P.Val"<0.1, df_sensitive$period=="early")]
  sig_sensitive_middle<-df_sensitive$SYMBOL[which(df_sensitive$"adj.P.Val"<0.1, df_sensitive$period=="middle")]
  sig_sensitive_late<-df_sensitive$SYMBOL[which(df_sensitive$"adj.P.Val"<0.1, df_sensitive$period=="late")]

  sig_common_early_middle<-intersect(sig_sensitive_early,sig_sensitive_middle)
  sig_common_early_late<-intersect(sig_sensitive_early,sig_sensitive_late)
  sig_common_middle_late<-intersect(sig_sensitive_middle,sig_sensitive_late)
  sig_share_early_middle<-union(sig_sensitive_early,sig_sensitive_middle)
  sig_share_early_late<-union(sig_sensitive_early,sig_sensitive_late)
  sig_share_middle_late<-union(sig_sensitive_middle,sig_sensitive_late)
  
  jaccard_early_middle<-length(sig_common_early_middle)/length(sig_share_early_middle)
  jaccard_early_late<-length(sig_common_early_late)/length(sig_share_early_late)
  jaccard_middle_late<-length(sig_common_middle_late)/length(sig_share_middle_late)

  print(paste0("Jaccard within sensitive ","early_middle:",jaccard_early_middle))
  print(paste0("Jaccard within sensitive ","early_late:",jaccard_early_late))
  print(paste0("Jaccard within sensitive ","middle_late:",jaccard_middle_late))

  ## log2fc spearman with ori
  order_gene<-df_ori$SYMBOL[order(df_ori$"adj.P.Val")]

  for (tri in unique(df_sensitive$period)){
    for (grp in unique(df_sensitive$group)){
      df_ori_tri<-df_ori %>% dplyr::filter(period==tri)
      df_sensitive_tri<-df_sensitive %>% dplyr::filter(period==tri, group==grp)
      df_ori_tri<-df_ori_tri[match(order_gene,df_ori_tri$SYMBOL),]
      df_sensitive_tri<-df_sensitive_tri[match(order_gene,df_sensitive_tri$SYMBOL),]
      log2fc_spearman<-cor.test(df_ori_tri$logFC,df_sensitive_tri$logFC,method = "spearman")
      print(paste0("log2fc spearman r with ori ",tri,"_",grp,":",log2fc_spearman$estimate))
      print(paste0("log2fc spearman p with ori ",tri,"_",grp,":",log2fc_spearman$p.value))
    }
  }

  ## lof2fc direction with ori
  for (tri in unique(df_sensitive$period)){
    for (grp in unique(df_sensitive$group)){
      df_ori_tri<-df_ori %>% dplyr::filter(period==tri)
      df_sensitive_tri<-df_sensitive %>% dplyr::filter(period==tri, group==grp)
      df_ori_tri<-df_ori_tri[match(order_gene,df_ori_tri$SYMBOL),]
      df_sensitive_tri<-df_sensitive_tri[match(order_gene,df_sensitive_tri$SYMBOL),]
      #log2fc_direction_concord<-sum(sign(df_ori_tri$logFC)==sign(df_sensitive_tri$logFC))/nrow(df_ori_tri)
      #log2fc_direction_flip<-sum(sign(df_ori_tri$logFC)!=sign(df_sensitive_tri$logFC))/nrow(df_ori_tri)

      m<-merge(df_ori_tri[,c("SYMBOL","logFC")],df_sensitive_tri[,c("SYMBOL","logFC")],by="SYMBOL")
      signA <- sign(m$logFC.x)
      signB <- sign(m$logFC.y)
      tab <- table(factor(signA, levels = c(-1, 1), labels = c("-1", "1")),
               factor(signB, levels = c(-1, 1), labels = c("-1", "1")))
      
      N       <- sum(tab)
      concord <- (tab[1,1] + tab[2,2]) / N  
      flip    <- (tab[1,2] + tab[2,1]) / N  
      
      a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
      Po <- (a + d) / N                                          # 观测一致
      pA_up  <- (c + d) / N; pA_down <- (a + b) / N
      pB_up  <- (b + d) / N; pB_down <- (a + c) / N
      Pe <- pA_up*pB_up + pA_down*pB_down                         # 随机一致
      kappa <- (Po - Pe) / (1 - Pe)

      print(paste0("log2fc direction concord with ori ",tri,"_",grp,":",concord))
      print(paste0("log2fc direction flip with ori ",tri,"_",grp,":",flip))
      print(paste0("log2fc direction kappa with ori ",tri,"_",grp,":",kappa))
    }
  }
   

}
