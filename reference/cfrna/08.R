library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)

files<-dir("BMI_sensitive/")
coldata_rna<-read.csv("__PATH_001__", row.names = 1) %>%
  mutate(BMI_group=case_when(
    BMI_mo_Anthropometric < 18.5 ~ "underweight",
    BMI_mo_Anthropometric >= 18.5 & BMI_mo_Anthropometric < 24 ~ "normal",
    BMI_mo_Anthropometric >= 24 & BMI_mo_Anthropometric < 28 ~ "overweight",
    BMI_mo_Anthropometric >= 28 ~ "obese",
    TRUE ~ NA_character_),
    trimester=case_when(trimester=="1"~"1st",
                        trimester=="2"~"2nd",
                        trimester=="3"~"3rd",
                        trimester=="all"~"overall"))

##### BMI sensitive ana  #####

## mRNA
res_mRNA<-data.frame()
for (i in files[grep("mRNA.*csv",files)]){
  group<-gsub(".*_","",i)
  group<-gsub(".csv","",group)
  tmp<-read.csv(paste0("BMI_sensitive/", i)) %>%
    mutate(group=group,RNA_type="mRNA")
  res_mRNA<-rbind(res_mRNA, tmp)
}
res_mRNA<-res_mRNA %>%
  mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall")) %>%
  mutate(group=gsub("wasting","underweight",group))

res_mRNA_ori<-read.csv("limma_349_batch_only.csv") %>%
    mutate(RNA_type="mRNA") %>%
  mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))

## miRNA
res_miRNA<-data.frame()
for (i in files[grep("miRNA.*csv",files)]){
  group<-gsub(".*_","",i)
  group<-gsub(".csv","",group)
  tmp<-read.csv(paste0("BMI_sensitive/", i)) %>%
    mutate(group=group,RNA_type="miRNA")
  res_miRNA<-rbind(res_miRNA, tmp)
}
res_miRNA<-res_miRNA %>%
  mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall")) %>%
  mutate(group=gsub("wasting","underweight",group))

res_miRNA_ori<-read_excel("manu_v2.0/manu_v2.0/Supplementary Table 2_Differential_analysis_results_for_cf-miRNA.xlsx") %>%
  mutate(RNA_type="miRNA") %>%
  dplyr::rename(period=trimester)

## lncRNA
res_lncRNA<-data.frame()
for (i in files[grep("lncRNA.*csv",files)]){
  group<-gsub(".*_","",i)
  group<-gsub(".csv","",group)
  tmp<-read.csv(paste0("BMI_sensitive/", i)) %>%
    mutate(group=group,RNA_type="lncRNA")
  res_lncRNA<-rbind(res_lncRNA, tmp)
}
res_lncRNA<-res_lncRNA %>%
  mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall")) %>%
  mutate(group=gsub("wasting","underweight",group))

res_lncRNA_ori<-data.frame()
for (i in 1:4){
tmp<-read_excel("limma_349_lncRNA_batch_only.xlsx",sheet=i)
  res_lncRNA_ori<-rbind(res_lncRNA_ori, tmp)
}
res_lncRNA_ori<-res_lncRNA_ori %>%
  mutate(RNA_type="lncRNA") %>%
  mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))


#res_allRNA<-rbind(res_mRNA,res_miRNA,res_lncRNA)



test_sensitive_BMI<-function(df_ori,df_sensitive,model_name="design 3 model",RNA_type="mRNA"){
  ## collect per-contrast statistics into a dataframe
  stat_df<-data.frame()

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

  ## Spearman and Kappa between sensitive results and ori for each period
  p<-list()

  for (tri in unique(df_sensitive$period)){
    stat_df_tri<-data.frame()
    m_merge<-data.frame()
    for (grp in unique(df_sensitive$group)){
      df_ori_tri<-df_ori %>%
        dplyr::filter(period==tri, P.Value<0.1) %>%
        dplyr::select(SYMBOL, logFC)
      df_sensitive_tri<-df_sensitive %>%
        dplyr::filter(period==tri, group==grp, P.Value<0.1) %>%
        dplyr::select(SYMBOL, logFC)

      m<-merge(df_ori_tri,df_sensitive_tri,by="SYMBOL")
      m<-m[complete.cases(m),]
      m_tmp<-m %>%
        dplyr::mutate(trimester=tri, group=grp)
      m_merge<-rbind(m_merge, m_tmp)

      ## log2fc Spearman with ori; only retain P.Value < 0.1 in both datasets
      spearman_valid <- nrow(m) >= 3 &&
        length(unique(m$logFC.x)) >= 2 &&
        length(unique(m$logFC.y)) >= 2
      if (spearman_valid){
        log2fc_spearman <- cor.test(m$logFC.x, m$logFC.y, method="spearman")
        spearman_r <- unname(log2fc_spearman$estimate)
        spearman_p <- log2fc_spearman$p.value
        print(paste0("log2fc spearman r with ori ",tri,"_",grp,":",spearman_r))
        print(paste0("log2fc spearman p with ori ",tri,"_",grp,":",spearman_p))
      } else {
        spearman_r <- NA_real_
        spearman_p <- NA_real_
        print(paste0("log2fc spearman skipped with ori ",tri,"_",grp,": insufficient valid data"))
      }

      ## log2fc direction with ori and kappa; only retain P.Value < 0.1 in both datasets
      if (nrow(m)==0){
        print(paste0("log2fc direction skipped with ori ",tri,"_",grp,": no overlapping valid data"))
        next
      }
      signA <- sign(m$logFC.x)
      signB <- sign(m$logFC.y)
      tab <- table(factor(signA, levels = c(-1, 1), labels = c("-1", "1")),
                   factor(signB, levels = c(-1, 1), labels = c("-1", "1")))

      N       <- sum(tab)
      concord <- (tab[1,1] + tab[2,2]) / N
      flip    <- (tab[1,2] + tab[2,1]) / N

      a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
      Po <- (a + d) / N
      pA_up  <- (c + d) / N; pA_down <- (a + b) / N
      pB_up  <- (b + d) / N; pB_down <- (a + c) / N
      Pe <- pA_up*pB_up + pA_down*pB_down
      kappa <- if (Pe==1) NA_real_ else (Po - Pe) / (1 - Pe)

      print(paste0("log2fc direction concord with ori ",tri,"_",grp,":",concord))
      print(paste0("log2fc direction flip with ori ",tri,"_",grp,":",flip))
      print(paste0("log2fc direction kappa with ori ",tri,"_",grp,":",kappa))

      stat_df_tmp<-data.frame(
        RNA_type=RNA_type,
        comparison_type="sensitive_vs_ori",
        trimester=tri,
        group=grp,
        contrast=paste0("ori_vs_",grp),
        spearman_r=spearman_r,
        spearman_p=spearman_p,
        kappa=kappa,
        concordance=concord,
        flip=flip,
        n_overlap=nrow(m),
        stringsAsFactors=FALSE)
      
      stat_df_tri<-rbind(stat_df_tri,stat_df_tmp)

      stat_label <- paste0(
        RNA_type,
        "\nSpearman r = ", ifelse(is.na(spearman_r), "NA", sprintf("%.3f", spearman_r)),
        "\nSpearman p = ", ifelse(is.na(spearman_p), "NA", format.pval(spearman_p, digits=3, eps=0.001)),
        "\nlog2FC direction concord = ", ifelse(is.na(concord), "NA", sprintf("%.3f", concord)),
        "\nKappa = ", ifelse(is.na(kappa), "NA", sprintf("%.3f", kappa))
      )

      p_tmp <- ggplot(m, aes(x=logFC.x, y=logFC.y))+
        geom_point()+
        geom_smooth(method="lm", se=FALSE, color="blue")+
        labs(title=paste0(RNA_type," | ",model_name," vs primary model: ",tri,"_",grp),
             x="log2FC in ori", y="log2FC in sensitive")+
        geom_hline(yintercept=0, color="black", linetype="solid")+
        geom_vline(xintercept=0, color="black", linetype="solid")+
        annotate("text", x=-Inf, y=Inf, label=stat_label,
                 hjust=-0.05, vjust=1.2, size=3.5, color="black")+
        theme(
          panel.background=element_blank(),
          plot.background=element_blank())
      #ggsave(filename = paste0("BMI_sensitive/",RNA_type,"_",tri,"_",grp,".pdf",sep=""), plot = p[[paste(tri, grp, sep="_")]], width = 6, height = 5)
    }
    stat_df<-rbind(stat_df,stat_df_tri)


    lab_list <- list()
      for (grp in unique(stat_df_tri$group)) {
        sub<-stat_df_tri %>% dplyr::filter(trimester==tri,group==grp)
        sr<-sub$spearman_r[1]
        sp<-sub$spearman_p[1]
        sc<-sub$concordance[1]
        sk<-sub$kappa[1]
        lab_list[[grp]] <- data.frame(
          trimester = tri,
          group = grp,
          stat_label = paste0(
            grp, "\n",
            "r = ", sprintf("%.3f", unname(sr)), "\n",
            "p = ", format.pval(sp, digits = 3, eps = 0.001),"\n",
            "concord = ", sprintf("%.3f", unname(sc)), "\n",
            "kappa = ", sprintf("%.3f", unname(sk))
          ),
          stringsAsFactors = FALSE
        )
      }
    lab_df <- bind_rows(lab_list)

    p_all<-ggplot(m_merge, aes(x=logFC.x, y=logFC.y,color=group))+
        geom_point(alpha=0.4)+
        geom_smooth(method="lm", se=FALSE,alpha=0.6)+
        labs(title=paste0(RNA_type," | BMI sensitivity analysis: ",tri),
            x="log2FC in primary model", y=paste0("log2FC in ",model_name))+
        geom_hline(yintercept=0, color="black", linetype="solid")+
        geom_vline(xintercept=0, color="black", linetype="solid")+
        geom_abline(intercept=0, slope=1,color = "grey40", linetype = "dashed")+
        geom_text_repel(
                data=lab_df, 
                aes(x=-Inf,y=Inf,label=stat_label,color=group),
                inherit.aes=FALSE,
                hjust=-0.05, 
                vjust=1.2, 
                size=2.5, 
                segment.color = NA,
                show.legend = FALSE)+
        theme(panel.background=element_blank(), plot.background=element_blank())

    #ggsave(p_all,filename=paste0("BMI_sensitive/",RNA_type,"_",tri,"_allinone.pdf"), width=6, height=5)
    p[[tri]]<-p_all
  }

  return(list(plots=p, stats=stat_df))
}






p_mRNA<-test_sensitive_BMI(res_mRNA_ori,res_mRNA,"design 3 model",RNA_type="mRNA")
p_miRNA<-test_sensitive_BMI(res_miRNA_ori,res_miRNA,"design 3 model",RNA_type="miRNA")
p_lncRNA<-test_sensitive_BMI(res_lncRNA_ori,res_lncRNA,"design 3 model",RNA_type="lncRNA")

### plot merged plots for all contrasts in one figure
#cov<-"design 3 model"
all_plots_mRNA<-list(
  p_mRNA$plots$`1st`, p_mRNA$plots$`2nd`, p_mRNA$plots$`3rd`, p_mRNA$plots$`overall`
)
p_tmp_mRNA <- wrap_plots(all_plots_mRNA, nrow = 2, ncol = 2) +
  plot_annotation(title = paste0("Differential mRNA of BMI stratums")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10))
ggsave(p_tmp_mRNA,filename="BMI_sensitive/BMI_mRNA_all.pdf",width=12,height=10)
  
all_plots_miRNA<-list(
  p_miRNA$plots$`1st`, p_miRNA$plots$`2nd`, p_miRNA$plots$`3rd`, p_miRNA$plots$`overall`
)
p_tmp_miRNA <- wrap_plots(all_plots_miRNA, nrow = 2, ncol = 2) +
  plot_annotation(title = paste0("Differential miRNA of BMI stratums")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10))
ggsave(p_tmp_miRNA,filename="BMI_sensitive/BMI_miRNA_all.pdf",width=12,height=10)

all_plots_lncRNA<-list(
  p_lncRNA$plots$`1st`, p_lncRNA$plots$`2nd`, p_lncRNA$plots$`3rd`, p_lncRNA$plots$`overall`
)
p_tmp_lncRNA <- wrap_plots(all_plots_lncRNA, nrow = 2, ncol = 2) +
  plot_annotation(title = paste0("Differential lncRNA of BMI stratums")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10))
ggsave(p_tmp_lncRNA,filename="BMI_sensitive/BMI_lncRNA_all.pdf",width=12,height=10)


## output all contrast statistics from design 3 (BMI) model to a new file
BMI_stat_all<-rbind(p_mRNA$stats, p_miRNA$stats, p_lncRNA$stats)
write.csv(BMI_stat_all,file="BMI_sensitive/BMI_design3_contrast_statistics.csv",row.names=F)
 
all_plots_mRNA<-list(p_mRNA$plots$`1st_underweight`, p_mRNA$plots$`1st_normal`, p_mRNA$plots$`1st_overweight`, p_mRNA$plots$`1st_obese`,
                p_mRNA$plots$`2nd_underweight`, p_mRNA$plots$`2nd_normal`, p_mRNA$plots$`2nd_overweight`, p_mRNA$plots$`2nd_obese`,
                p_mRNA$plots$`3rd_underweight`, p_mRNA$plots$`3rd_normal`, p_mRNA$plots$`3rd_overweight`, p_mRNA$plots$`3rd_obese`,
                p_mRNA$plots$`overall_underweight`, p_mRNA$plots$`overall_normal`, p_mRNA$plots$`overall_overweight`, p_mRNA$plots$`overall_obese`)
all_plot_miRNA<-list(p_miRNA$plots$`1st_underweight`, p_miRNA$plots$`1st_normal`, p_miRNA$plots$`1st_overweight`, p_miRNA$plots$`1st_obese`,
                p_miRNA$plots$`2nd_underweight`, p_miRNA$plots$`2nd_normal`, p_miRNA$plots$`2nd_overweight`, p_miRNA$plots$`2nd_obese`,
                p_miRNA$plots$`3rd_underweight`, p_miRNA$plots$`3rd_normal`, p_miRNA$plots$`3rd_overweight`, p_miRNA$plots$`3rd_obese`,
                p_miRNA$plots$`overall_underweight`, p_miRNA$plots$`overall_normal`, p_miRNA$plots$`overall_overweight`, p_miRNA$plots$`overall_obese`)
all_plot_lncRNA<-list(p_lncRNA$plots$`1st_underweight`, p_lncRNA$plots$`1st_normal`, p_lncRNA$plots$`1st_overweight`, p_lncRNA$plots$`1st_obese`,
                p_lncRNA$plots$`2nd_underweight`, p_lncRNA$plots$`2nd_normal`, p_lncRNA$plots$`2nd_overweight`, p_lncRNA$plots$`2nd_obese`,
                p_lncRNA$plots$`3rd_underweight`, p_lncRNA$plots$`3rd_normal`, p_lncRNA$plots$`3rd_overweight`, p_lncRNA$plots$`3rd_obese`,
                p_lncRNA$plots$`overall_underweight`, p_lncRNA$plots$`overall_normal`, p_lncRNA$plots$`overall_overweight`, p_lncRNA$plots$`overall_obese`)

cov<-"design 3 model"
p_tmp_mRNA <- wrap_plots(all_plots_mRNA, nrow = 4, ncol = 4) +
  plot_annotation(title = paste0("Differential mRNA ",cov, " vs primary model")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
ggsave(p_tmp_mRNA,filename="BMI_sensitive/BMI_mRNA_all.pdf",width=20,height=15)

p_tmp_miRNA <- wrap_plots(all_plot_miRNA, nrow = 4, ncol = 4) +
  plot_annotation(title = paste0("Differential miRNA ",cov, " vs primary model")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
ggsave(p_tmp_miRNA,filename="BMI_sensitive/BMI_miRNA_all.pdf",width=20,height=15)

p_tmp_lncRNA <- wrap_plots(all_plot_lncRNA, nrow = 4, ncol = 4) +
  plot_annotation(title = paste0("Differential lncRNA ",cov, " vs primary model")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
ggsave(p_tmp_lncRNA,filename="BMI_sensitive/BMI_lncRNA_all.pdf",width=20,height=15)



##### full vs ori #####

## BMI full model
res_mRNA_full<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_mRNA_full.xlsx",sheet = i) %>%
    mutate(RNA_type="mRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_mRNA_full<-rbind(res_mRNA_full,tmp)
}

res_miRNA_full<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_miRNA_full.xlsx",sheet = i) %>%
    mutate(RNA_type="miRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_miRNA_full<-rbind(res_miRNA_full,tmp)
}

res_lncRNA_full<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_lncRNA_full.xlsx",sheet = i) %>%
    mutate(RNA_type="lncRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_lncRNA_full<-rbind(res_lncRNA_full,tmp)
}


## metabolite model
res_mRNA_metabolite<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_mRNA_metabolite.xlsx",sheet = i) %>%
    mutate(RNA_type="mRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_mRNA_metabolite<-rbind(res_mRNA_metabolite,tmp)
}

res_miRNA_metabolite<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_miRNA_metabolite.xlsx",sheet = i) %>%
    mutate(RNA_type="miRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_miRNA_metabolite<-rbind(res_miRNA_metabolite,tmp)
}

res_lncRNA_metabolite<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_lncRNA_metabolite.xlsx",sheet = i) %>%
    mutate(RNA_type="lncRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_lncRNA_metabolite<-rbind(res_lncRNA_metabolite,tmp)
}


## obstetric model
res_mRNA_obstetric<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_mRNA_obstetric.xlsx",sheet = i) %>%
    mutate(RNA_type="mRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_mRNA_obstetric<-rbind(res_mRNA_obstetric,tmp)
}

res_miRNA_obstetric<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_miRNA_obstetric.xlsx",sheet = i) %>%
    mutate(RNA_type="miRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_miRNA_obstetric<-rbind(res_miRNA_obstetric,tmp)
}

res_lncRNA_obstetric<-data.frame()
for (i in 1:4) {
  tmp<-read_excel("covriate_sensitive/limma_349_lncRNA_obstetric.xlsx",sheet = i) %>%
    mutate(RNA_type="lncRNA") %>%
    mutate(period=case_when(period=="1"~"1st",
                  period=="2"~"2nd",
                  period=="3"~"3rd",
                  period=="all"~"overall"))
  res_lncRNA_obstetric<-rbind(res_lncRNA_obstetric,tmp)
}


test_sensitive<-function(df_sensitive,df_ori,model_name="full model",RNA_type="mRNA"){
  ## collect per-contrast statistics into a dataframe
  stat_df<-data.frame()

  ## Jaccard with ori for each period
  for (tri in unique(df_sensitive$period)){
    df_ori_tri <- df_ori %>% dplyr::filter(period==tri)
    df_sensitive_tri <- df_sensitive %>% dplyr::filter(period==tri)
    sig_ori_tri <- df_ori_tri$SYMBOL[df_ori_tri$adj.P.Val<0.1]
    sig_sensitive_tri <- df_sensitive_tri$SYMBOL[df_sensitive_tri$adj.P.Val<0.1]
    sig_common <- intersect(sig_ori_tri, sig_sensitive_tri)
    sig_share <- union(sig_ori_tri, sig_sensitive_tri)
    jaccard <- if (length(sig_share)==0) NA_real_ else {
      length(sig_common) / length(sig_share)
    }
    print(paste0("Jaccard with ori ",tri,":",jaccard))
  }

  ## Jaccard between different periods within sensitive results
  sensitive_periods <- unique(df_sensitive$period)
  if (length(sensitive_periods)>=2){
    period_pairs <- combn(sensitive_periods, 2, simplify=FALSE)
    for (pair in period_pairs){
      period1 <- pair[1]
      period2 <- pair[2]
      sig_period1 <- df_sensitive$SYMBOL[
        df_sensitive$adj.P.Val<0.1 & df_sensitive$period==period1]
      sig_period2 <- df_sensitive$SYMBOL[
        df_sensitive$adj.P.Val<0.1 & df_sensitive$period==period2]
      sig_share <- union(sig_period1, sig_period2)
      jaccard <- if (length(sig_share)==0) NA_real_ else {
        length(intersect(sig_period1, sig_period2)) / length(sig_share)
      }
      print(paste0("Jaccard within sensitive ",period1,"_",period2,":",jaccard))
    }
  }

  ## Spearman and Kappa between sensitive results and ori for each period
  p <- list()
  m_merge<-data.frame()
  for (tri in unique(df_sensitive$period)){
    df_ori_tri <- df_ori %>%
      dplyr::filter(period==tri, P.Value<0.1) %>%
      dplyr::select(SYMBOL, logFC)
    df_sensitive_tri <- df_sensitive %>%
      dplyr::filter(period==tri, P.Value<0.1) %>%
      dplyr::select(SYMBOL, logFC)

    m <- merge(df_ori_tri, df_sensitive_tri, by="SYMBOL")
    m <- m[complete.cases(m),]
    m_tmp<-m %>%
      mutate(trimester=tri)
    m_merge<-rbind(m_merge,m_tmp)

    spearman_valid <- nrow(m)>=3 &&
      length(unique(m$logFC.x))>=2 &&
      length(unique(m$logFC.y))>=2
    if (spearman_valid){
      log2fc_spearman <- cor.test(m$logFC.x, m$logFC.y, method="spearman")
      spearman_r <- unname(log2fc_spearman$estimate)
      spearman_p <- log2fc_spearman$p.value
      print(paste0("log2fc spearman r with ori ",tri,":",spearman_r))
      print(paste0("log2fc spearman p with ori ",tri,":",spearman_p))
    } else {
      spearman_r <- NA_real_
      spearman_p <- NA_real_
      print(paste0("log2fc spearman skipped with ori ",tri,
                   ": insufficient valid data"))
    }

    if (nrow(m)==0){
      print(paste0("log2fc direction skipped with ori ",tri,
                   ": no overlapping valid data"))
      next
    }

    m_kappa <- m[sign(m$logFC.x)!=0 & sign(m$logFC.y)!=0,]
    if (nrow(m_kappa)==0){
      kappa <- NA_real_
      concord <- NA_real_
      flip <- NA_real_
    } else {
      signA <- sign(m_kappa$logFC.x)
      signB <- sign(m_kappa$logFC.y)
      tab <- table(
        factor(signA, levels=c(-1,1), labels=c("-1","1")),
        factor(signB, levels=c(-1,1), labels=c("-1","1"))
      )
      N <- sum(tab)
      concord <- (tab[1,1] + tab[2,2]) / N
      flip    <- (tab[1,2] + tab[2,1]) / N
      a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
      Po <- (a+d) / N
      pA_up <- (c+d) / N; pA_down <- (a+b) / N
      pB_up <- (b+d) / N; pB_down <- (a+c) / N
      Pe <- pA_up*pB_up + pA_down*pB_down
      kappa <- if (Pe==1) NA_real_ else (Po-Pe) / (1-Pe)
    }
    print(paste0("log2fc direction kappa with ori ",tri,":",kappa))
    print(paste0("log2fc direction concord with ori ",tri,"_",":",concord))

    stat_df<-rbind(stat_df, data.frame(
      RNA_type=RNA_type,
      model=model_name,
      trimester=tri,
      contrast="sensitive_vs_ori",
      spearman_r=spearman_r,
      spearman_p=spearman_p,
      kappa=kappa,
      concordance=concord,
      flip=flip,
      n_overlap=nrow(m),
      stringsAsFactors=FALSE))

    stat_label <- paste0(
      RNA_type,
      "\nSpearman r = ", ifelse(is.na(spearman_r), "NA", sprintf("%.3f", spearman_r)),
      "\nSpearman p = ", ifelse(is.na(spearman_p), "NA",format.pval(spearman_p, digits=3, eps=0.001)),
      "\nlog2FC direction concord = ", ifelse(is.na(concord), "NA", sprintf("%.3f", concord)),
      "\nKappa = ", ifelse(is.na(kappa), "NA", sprintf("%.3f", kappa))
    )
    p[[tri]] <- ggplot(m, aes(x=logFC.x, y=logFC.y))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, color="blue")+
      labs(title=paste0(RNA_type," | ",model_name," vs primary model: ",tri),
           x="log2FC in ori", y="log2FC in sensitive")+
      geom_hline(yintercept=0, color="black", linetype="solid")+
      geom_vline(xintercept=0, color="black", linetype="solid")+
      annotate("text", x=-Inf, y=Inf, label=stat_label,
               hjust=-0.05, vjust=1.2, size=3.5, color="black")+
      theme(panel.background=element_blank(), plot.background=element_blank())

    #ggsave(
    #  filename=paste0("covriate_sensitive/",RNA_type,"_",tri,".pdf"),
    #  plot=p[[tri]], width=6, height=5
    #)
  }
  
  lab_list <- list()
    for (tri in unique(stat_df$trimester)) {
      sub<-stat_df %>% dplyr::filter(trimester==tri)
      sr<-sub$spearman_r[1]
      sp<-sub$spearman_p[1]
      sc<-sub$concordance[1]
      sk<-sub$kappa[1]
      lab_list[[tri]] <- data.frame(
        trimester = tri,
        stat_label = paste0(
          tri, "\n",
          "r = ", sprintf("%.3f", unname(sr)), "\n",
          "p = ", format.pval(sp, digits = 3, eps = 0.001),"\n",
          "concord = ", sprintf("%.3f", unname(sc)), "\n",
          "kappa = ", sprintf("%.3f", unname(sk))
        ),
        stringsAsFactors = FALSE
      )
    }
  lab_df <- bind_rows(lab_list)

  p_all<-ggplot(m_merge, aes(x=logFC.x, y=logFC.y,color=trimester))+
      geom_point(alpha=0.4)+
      geom_smooth(method="lm", se=FALSE,alpha=0.6)+
      labs(title=paste0(RNA_type," | ",model_name," vs primary model: "),
           x="log2FC in primary model", y=paste0("log2FC in ",model_name))+
      geom_hline(yintercept=0, color="black", linetype="solid")+
      geom_vline(xintercept=0, color="black", linetype="solid")+
      geom_abline(intercept=0, slope=1,color = "grey40", linetype = "dashed")+
      geom_text_repel(
               data=lab_df, 
               aes(x=-Inf,y=Inf,label=stat_label,color=trimester),
               inherit.aes=FALSE,
               hjust=-0.05, 
               vjust=1.2, 
               size=2.5, 
               segment.color = NA,
               show.legend = FALSE)+
      theme(panel.background=element_blank(), plot.background=element_blank())

  #ggsave(p_all,filename=paste0("covriate_sensitive/",RNA_type,"_",model_name,"_allinone.pdf"), width=6, height=5)
  
  return(list(plots=p_all, stats=stat_df,merge=m_merge))
}



p_mRNA<-test_sensitive(res_mRNA_metabolite,res_mRNA_ori,"design 4 model",RNA_type = "mRNA")
p_miRNA<-test_sensitive(res_miRNA_metabolite,res_miRNA_ori,"design 4 model",RNA_type = "miRNA")
p_lncRNA<-test_sensitive(res_lncRNA_metabolite,res_lncRNA_ori,"design 4 model",RNA_type = "lncRNA")

## output all contrast statistics from design 4 (covariate) models to a new file
cov_stat_all<-rbind(p_mRNA$stats, p_miRNA$stats, p_lncRNA$stats)
write.csv(cov_stat_all,file="covriate_sensitive/covariate_design4_contrast_statistics.csv",row.names=F)


## merged all res



## single plot
cov<-"design 4 model"
all_plots <- list(
  p_mRNA$plots, p_miRNA$plots, p_lncRNA$plots
)

p_tmp <- wrap_plots(all_plots, nrow = 1, ncol = 3) +
  plot_annotation(title = paste0("Differential RNA ",cov," vs primary model")) &
  theme(plot.title = element_text(hjust = 0.5, size = 15))

ggsave(p_tmp,filename=paste0("covriate_sensitive/RNA_",cov,"_vs_design1_model.pdf"),width=18,height=5)


##### sensitive sig genes #####
ori_sig_mRNA<-res_mRNA_ori$SYMBOL[which(res_mRNA_ori$adj.P.Val<0.1 & abs(res_mRNA_ori$logFC)>0.15)]
ori_sig_miRNA<-res_miRNA_ori$SYMBOL[which(res_miRNA_ori$adj.P.Val<0.1 & abs(res_miRNA_ori$logFC)>0.15)]
ori_sig_lncRNA<-res_lncRNA_ori$SYMBOL[which(res_lncRNA_ori$adj.P.Val<0.1 & abs(res_lncRNA_ori$logFC)>0.15)]

full_sig_mRNA<-res_mRNA_full$SYMBOL[which(res_mRNA_full$adj.P.Val<0.1 & abs(res_mRNA_full$logFC)>0.15)]
full_sig_miRNA<-res_miRNA_full$SYMBOL[which(res_miRNA_full$adj.P.Val<0.1 & abs(res_miRNA_full$logFC)>0.15)]
full_sig_lncRNA<-res_lncRNA_full$SYMBOL[which(res_lncRNA_full$adj.P.Val<0.1 & abs(res_lncRNA_full$logFC)>0.15)]
print(paste0("Number of significant mRNA in full model: ", length(full_sig_mRNA)))
print(paste0("Number of significant miRNA in full model: ", length(full_sig_miRNA)))
print(paste0("Number of significant lncRNA in full model: ", length(full_sig_lncRNA)))
print(paste0("same significant mRNA in full model and ori: ", length(intersect(full_sig_mRNA, ori_sig_mRNA))))
print(paste0("same significant miRNA in full model and ori: ", length(intersect(full_sig_miRNA, ori_sig_miRNA))))
print(paste0("same significant lncRNA in full model and ori: ", length(intersect(full_sig_lncRNA, ori_sig_lncRNA))))
print(paste0("new significant mRNA in full model: ", length(setdiff(full_sig_mRNA, ori_sig_mRNA))))
print(paste0("new significant miRNA in full model: ", length(setdiff(full_sig_miRNA, ori_sig_miRNA))))
print(paste0("new significant lncRNA in full model: ", length(setdiff(full_sig_lncRNA, ori_sig_lncRNA))))

bmi_sig_mRNA<-res_mRNA$SYMBOL[which(res_mRNA$adj.P.Val<0.1 & abs(res_mRNA$logFC)>0.15)]
bmi_sig_miRNA<-res_miRNA$SYMBOL[which(res_miRNA$adj.P.Val<0.1 & abs(res_miRNA$logFC)>0.15)]
bmi_sig_lncRNA<-res_lncRNA$SYMBOL[which(res_lncRNA$adj.P.Val<0.1 & abs(res_lncRNA$logFC)>0.15)]
print(paste0("Number of significant mRNA in BMI model: ", length(bmi_sig_mRNA)))
print(paste0("Number of significant miRNA in BMI model: ", length(bmi_sig_miRNA)))
print(paste0("Number of significant lncRNA in BMI model: ", length(bmi_sig_lncRNA)))
print(paste0("same significant mRNA in BMI model and ori: ", length(intersect(bmi_sig_mRNA, ori_sig_mRNA))))
print(paste0("same significant miRNA in BMI model and ori: ", length(intersect(bmi_sig_miRNA, ori_sig_miRNA))) )
print(paste0("same significant lncRNA in BMI model and ori: ", length(intersect(bmi_sig_lncRNA, ori_sig_lncRNA))) )
print(paste0("new significant mRNA in BMI model: ", length(setdiff(bmi_sig_mRNA, ori_sig_mRNA))))
print(paste0("new significant miRNA in BMI model: ", length(setdiff(bmi_sig_miRNA, ori_sig_miRNA))))
print(paste0("new significant lncRNA in BMI model: ", length(setdiff(bmi_sig_lncRNA, ori_sig_lncRNA))))


metabolite_sig_mRNA<-res_mRNA_metabolite$SYMBOL[which(res_mRNA_metabolite$adj.P.Val<0.1 & abs(res_mRNA_metabolite$logFC)>0.15)]
metabolite_sig_miRNA<-res_miRNA_metabolite$SYMBOL[which(res_miRNA_metabolite$adj.P.Val<0.1 & abs(res_miRNA_metabolite$logFC)>0.15)]
metabolite_sig_lncRNA<-res_lncRNA_metabolite$SYMBOL[which(res_lncRNA_metabolite$adj.P.Val<0.1 & abs(res_lncRNA_metabolite$logFC)>0.15)]
print(paste0("Number of significant mRNA in metabolite model: ", length(metabolite_sig_mRNA)))
print(paste0("Number of significant miRNA in metabolite model: ", length(metabolite_sig_miRNA)))
print(paste0("Number of significant lncRNA in metabolite model: ", length(metabolite_sig_lncRNA)))
print(paste0("same significant mRNA in metabolite model and ori: ", length(intersect(metabolite_sig_mRNA, ori_sig_mRNA))))
print(paste0("same significant miRNA in metabolite model and ori: ", length(intersect(metabolite_sig_miRNA, ori_sig_miRNA))) )
print(paste0("same significant lncRNA in metabolite model and ori: ", length(intersect(metabolite_sig_lncRNA, ori_sig_lncRNA))) )
print(paste0("new significant mRNA in metabolite model: ", length(setdiff(metabolite_sig_mRNA, ori_sig_mRNA))))
print(paste0("new significant miRNA in metabolite model: ", length(setdiff(metabolite_sig_miRNA, ori_sig_miRNA))))
print(paste0("new significant lncRNA in metabolite model: ", length(setdiff(metabolite_sig_lncRNA, ori_sig_lncRNA))))

obstetric_sig_mRNA<-res_mRNA_obstetric$SYMBOL[which(res_mRNA_obstetric$adj.P.Val<0.1 & abs(res_mRNA_obstetric$logFC)>0.15)]
obstetric_sig_miRNA<-res_miRNA_obstetric$SYMBOL[which(res_miRNA_obstetric$adj.P.Val<0.1 & abs(res_miRNA_obstetric$logFC)>0.15)]
obstetric_sig_lncRNA<-res_lncRNA_obstetric$SYMBOL[which(res_lncRNA_obstetric$adj.P.Val<0.1 & abs(res_lncRNA_obstetric$logFC)>0.15)]
print(paste0("Number of significant mRNA in obstetric model: ", length(obstetric_sig_mRNA)))
print(paste0("Number of significant miRNA in obstetric model: ", length(obstetric_sig_miRNA)))
print(paste0("Number of significant lncRNA in obstetric model: ", length(obstetric_sig_lncRNA)))
print(paste0("same significant mRNA in obstetric model and ori: ", length(intersect(obstetric_sig_mRNA, ori_sig_mRNA))))
print(paste0("same significant miRNA in obstetric model and ori: ", length(intersect(obstetric_sig_miRNA, ori_sig_miRNA))) )
print(paste0("same significant lncRNA in obstetric model and ori: ", length(intersect(obstetric_sig_lncRNA, ori_sig_lncRNA))) )
print(paste0("new significant mRNA in obstetric model: ", length(setdiff(obstetric_sig_mRNA, ori_sig_mRNA))))
print(paste0("new significant miRNA in obstetric model: ", length(setdiff(obstetric_sig_miRNA, ori_sig_miRNA))))
print(paste0("new significant lncRNA in obstetric model: ", length(setdiff(obstetric_sig_lncRNA, ori_sig_lncRNA))))


res_mRNA_obstetric[which(res_mRNA_obstetric$adj.P.Val<0.1 & abs(res_mRNA_obstetric$logFC)>0.15),]
res_miRNA_obstetric[which(res_miRNA_obstetric$adj.P.Val<0.1 & abs(res_miRNA_obstetric$logFC)>0.15),]
res_lncRNA_obstetric[which(res_lncRNA_obstetric$adj.P.Val<0.1 & abs(res_lncRNA_obstetric$logFC)>0.15),]


###### output #####
tmp<-rbind(res_mRNA,res_miRNA,res_lncRNA)
write.csv(tmp,file="BMI_sensitive/Supplementary Table R2 Differential analysis of mRNA, miRNA and lncRNA of design 3 model.csv",row.names=F)

write.csv(res_sensitive,file="microbiome ko gsea复现/Supplementary Table R3 Differential analysis of microbiome KO of design 3 model.csv",row.names=F)

write.csv(rbind(res_mRNA_metabolite,res_miRNA_metabolite,res_lncRNA_metabolite),file="covriate_sensitive/Supplementary Table R4 Differential analysis of mRNA, miRNA and lncRNA of design 4 model with metabolite covariates.csv",row.names=F)

write.csv(rbind(res_mRNA_obstetric,res_miRNA_obstetric,res_lncRNA_obstetric),file="covriate_sensitive/Supplementary Table R5 Differential analysis of mRNA, miRNA and lncRNA of design 5 model with obstetric covariates.csv",row.names=F)

write.csv(res_microbiome_metabolite,file="microbiome ko gsea复现/Supplementary Table R6 Differential analysis of microbiome KO of design 4 model with metabolite covariates.csv",row.names=F)

write.csv(res_microbiome_obstetric,file="microbiome ko gsea复现/Supplementary Table R7 Differential analysis of microbiome KO of design 5 model with obstetric covariates.csv",row.names=F)


write.csv(rbind(res_mRNA_full,res_miRNA_full,res_lncRNA_full),file="covriate_sensitive/Supplementary Table R3 Differential analysis of mRNA, miRNA and lncRNA of design 3 model.csv",row.names=F)
