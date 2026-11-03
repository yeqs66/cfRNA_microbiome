library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)

load("__PATH_001__")
load("microbiome ko gsea复现/microbiome_sampleinfo_GDM1999s_20240318.rda")
matrix_meta<-read.csv("microbiome ko gsea复现/TJGDM1962s_ko_tri_95_rpk.csv")

df <- t(matrix_meta) %>% as.data.frame()
microbiome_sampleinfo_zinb <- clinical_347_mo_kid_pairs %>%
  dplyr::select(ID,BMI_mo,ethnicity_mo,smoking_mo,drinking_mo_18topreg) %>%
  merge(.,microbiome_sampleinfo, by="ID") %>%
  dplyr::arrange(sample) %>%
  dplyr::select(sample,ID,group,trimester,BMI_mo,ethnicity_mo,smoking_mo,drinking_mo_18topreg) %>%
  dplyr::filter(sample %in% row.names(df)) %>%
  mutate(trimester=factor(trimester,levels = c("early","middle","late")),
         ID=as.character(ID),group=factor(group,levels=c("GDM","Control")),
         BMI_mo=as.numeric(BMI_mo),ethnicity_mo=factor(ethnicity_mo),
         smoking_mo=factor(smoking_mo),drinking_mo_18topreg=factor(drinking_mo_18topreg),
        BMI_group=case_when(
          BMI_mo < 18.5 ~ "wasting",
          BMI_mo >= 18.5 & BMI_mo < 24 ~ "normal",
          BMI_mo >= 24 & BMI_mo < 28 ~ "overweight",
          BMI_mo >= 28 ~ "obese",
          TRUE ~ NA_character_))


res_ori<-data.frame()
for (tri in c("early","middle","late","overall")){
  df<-read.csv(paste("microbiome ko gsea复现/new_deg/TJGDM1962s_ko_tri_95.zinb",tri,"adj_rpk.csv",sep="_")) %>%
    mutate(trimester=tri)
  res_ori<-rbind(res_ori,df)
}

res_sensitive<-data.frame()
for (group in unique(microbiome_sampleinfo_zinb$BMI_group)){
  for (tri in c("early","middle","late","overall")){
    tmp<-read.csv(paste("microbiome ko gsea复现/sensitive/TJGDM1962s_ko_tri_95.zinb",group,tri,"adj_rpk.csv",sep="_")) %>%
      mutate(trimester=tri,BMI_group=group)
    res_sensitive<-rbind(res_sensitive,tmp)
  }
}



##### sensitive within #####
p_group <- list()
  sensitive_trimester <- unique(res_sensitive$trimester)
  sensitive_groups <- unique(res_sensitive$BMI_group)

  if (length(sensitive_groups) >= 2){
    group_pairs <- combn(sensitive_groups, 2, simplify=FALSE)

    for (tri in sensitive_trimester){
      for (pair in group_pairs){
        grp1 <- pair[1]
        grp2 <- pair[2]

        res_sensitive_grp1 <- res_sensitive %>%
          dplyr::filter(trimester==tri, BMI_group==grp1, p.value<0.9) %>%
          dplyr::select(feature, log2MedianRatio)
        res_sensitive_grp2 <- res_sensitive %>%
          dplyr::filter(trimester==tri, BMI_group==grp2, p.value<0.9) %>%
          dplyr::select(feature, log2MedianRatio)

        m_group <- merge(res_sensitive_grp1, res_sensitive_grp2, by="feature")
        m_group <- m_group[complete.cases(m_group),]

        spearman_valid <- nrow(m_group) >= 3 &&
          length(unique(m_group$log2MedianRatio.x)) >= 2 &&
          length(unique(m_group$log2MedianRatio.y)) >= 2
        if (spearman_valid){
          group_spearman <- cor.test(m_group$log2MedianRatio.x, m_group$log2MedianRatio.y,
                                     method="spearman")
          group_r <- unname(group_spearman$estimate)
          group_p <- group_spearman$p.value
        } else {
          group_r <- NA_real_
          group_p <- NA_real_
          print(paste0("group Spearman skipped ",tri,"_",grp1,"_vs_",grp2,
                       ": insufficient valid data"))
        }

        if (nrow(m_group) == 0){
          print(paste0("group Kappa skipped ",tri,"_",grp1,"_vs_",grp2,
                       ": no overlapping valid data"))
          next
        }

        m_group <- m_group[sign(m_group$log2MedianRatio.x) != 0 &
                             sign(m_group$log2MedianRatio.y) != 0,]
        if (nrow(m_group) == 0){
          group_kappa <- NA_real_
        } else {
          signA <- sign(m_group$log2MedianRatio.x)
          signB <- sign(m_group$log2MedianRatio.y)
          tab <- table(
            factor(signA, levels=c(-1, 1), labels=c("-1", "1")),
            factor(signB, levels=c(-1, 1), labels=c("-1", "1"))
          )
          N <- sum(tab)
          a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
          Po <- (a + d) / N
          pA_up <- (c + d) / N; pA_down <- (a + b) / N
          pB_up <- (b + d) / N; pB_down <- (a + c) / N
          Pe <- pA_up*pB_up + pA_down*pB_down
          group_kappa <- if (Pe==1) NA_real_ else (Po - Pe) / (1 - Pe)
        }

        print(paste0("group Spearman r ",tri,"_",grp1,"_vs_",grp2,":",group_r))
        print(paste0("group Spearman p ",tri,"_",grp1,"_vs_",grp2,":",group_p))
        print(paste0("group Kappa ",tri,"_",grp1,"_vs_",grp2,":",group_kappa))

        group_label <- paste0(
          "\nSpearman r = ", ifelse(is.na(group_r), "NA", sprintf("%.3f", group_r)),
          "\nSpearman p = ", ifelse(is.na(group_p), "NA",
                                      format.pval(group_p, digits=3, eps=0.001)),
          "\nKappa = ", ifelse(is.na(group_kappa), "NA",
                               sprintf("%.3f", group_kappa))
        )
        group_key <- paste(tri, grp1, "vs", grp2, sep="_")
        p_group[[group_key]] <- ggplot(m_group, aes(x=log2MedianRatio.x, y=log2MedianRatio.y))+
          geom_point()+
          geom_smooth(method="lm", se=FALSE, color="blue")+
          labs(title=paste0("log2FC between sensitive groups: ",tri,"_",grp1," vs ",grp2),
               x=paste0("log2FC in ",grp1), y=paste0("log2FC in ",grp2))+
          geom_hline(yintercept=0, color="black", linetype="solid")+
          geom_vline(xintercept=0, color="black", linetype="solid")+
          annotate("text", x=-Inf, y=Inf, label=group_label,
                   hjust=-0.05, vjust=1.2, size=3.5, color="black")+
          theme(
            panel.background=element_blank(),
            plot.background=element_blank())

        ggsave(
          filename=paste0("BMI_sensitive_within_sensitive_",tri,"_",grp1,
                          "_vs_",grp2,"_spearman_kappa.pdf"),
          plot=p_group[[group_key]], width=6, height=5
        )
      }
    }
  } else {
    print("group Spearman/Kappa skipped: fewer than two sensitive groups")
  }


##### ori vs sensitive #####
test_sensitive_BMI<-function(res_sensitive,res_ori,model_name="design 3 model"){
  ## Spearman and Kappa between sensitive results and ori for each period 
  p<-list()
  stat_df<-data.frame()
  for (tri in unique(res_sensitive$trimester)){
    stat_df_tri<-data.frame()
    m_merge<-data.frame()
    for (grp in unique(res_sensitive$BMI_group)){
      res_ori_tri<-res_ori %>%
        dplyr::filter(trimester==tri,p.value<0.9) %>%
        dplyr::select(feature, p.value,log2MedianRatio) 
      res_sensitive_tri<-res_sensitive %>%
        dplyr::filter(trimester==tri, BMI_group==grp,p.value<0.9) %>%
        dplyr::select(feature, p.value,log2MedianRatio)

## Jaccard with ori for each period
    sig_ori_tri <- res_ori_tri$feature[res_ori_tri$p.value<0.05 & abs(res_ori_tri$log2MedianRatio)>0.5]
    sig_sensitive_tri <- res_sensitive_tri$feature[res_sensitive_tri$p.value<0.05 & abs(res_sensitive_tri$log2MedianRatio)>0.5]
    sig_common <- intersect(sig_ori_tri, sig_sensitive_tri)
    sig_share <- union(sig_ori_tri, sig_sensitive_tri)
    jaccard <- if (length(sig_share)==0) NA_real_ else {
      length(sig_common) / length(sig_share)
    }
    print(paste0("Jaccard with ori ",tri,"_",grp,":",jaccard))


      
      m<-merge(res_ori_tri,res_sensitive_tri,by="feature")
      m<-m[complete.cases(m),]
      m_tmp<-m %>% mutate(group=grp)
      m_merge<-rbind(m_merge,m_tmp)

      ## log2fc Spearman with ori; only retain p.value < 0.1 in both datasets
      spearman_valid <- nrow(m) >= 3 &&
        length(unique(m$log2MedianRatio.x)) >= 2 &&
        length(unique(m$log2MedianRatio.y)) >= 2
      if (spearman_valid){
        log2MedianRatio_spearman <- cor.test(m$log2MedianRatio.x, m$log2MedianRatio.y, method="spearman")
        spearman_r <- unname(log2MedianRatio_spearman$estimate)
        spearman_p <- log2MedianRatio_spearman$p.value
        print(paste0("log2MedianRatio spearman r with ori ",tri,"_",grp,":",spearman_r))
        print(paste0("log2MedianRatio spearman p with ori ",tri,"_",grp,":",spearman_p))
      } else {
        spearman_r <- NA_real_
        spearman_p <- NA_real_
        print(paste0("log2MedianRatio spearman skipped with ori ",tri,"_",grp,": insufficient valid data"))
      }

      ## log2MedianRatio direction with ori and kappa; only retain p.value < 0.1 in both datasets
      if (nrow(m)==0){
        print(paste0("log2MedianRatio direction skipped with ori ",tri,"_",grp,": no overlapping valid data"))
        next
      }
      signA <- sign(m$log2MedianRatio.x)
      signB <- sign(m$log2MedianRatio.y)
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

      print(paste0("log2MedianRatio direction concord with ori ",tri,"_",grp,":",concord))
      print(paste0("log2MedianRatio direction flip with ori ",tri,"_",grp,":",flip))
      print(paste0("log2MedianRatio direction kappa with ori ",tri,"_",grp,":",kappa))

      stat_df_tmp<-data.frame(
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
        "\nJaccard = ", ifelse(is.na(jaccard), "NA", sprintf("%.3f", jaccard)),
        "\nSpearman r = ", ifelse(is.na(spearman_r), "NA", sprintf("%.3f", spearman_r)),
        "\nSpearman p = ", ifelse(is.na(spearman_p), "NA", format.pval(spearman_p, digits=3, eps=0.001)),
        "\nlog2MedianRatio direction concord = ", ifelse(is.na(concord), "NA", sprintf("%.3f", concord)),
        "\nKappa = ", ifelse(is.na(kappa), "NA", sprintf("%.3f", kappa))
      )

      p[[paste(tri, grp, sep="_")]] <- ggplot(m, aes(x=log2MedianRatio.x, y=log2MedianRatio.y))+
        geom_point()+
        geom_smooth(method="lm", se=FALSE, color="blue")+
        labs(title=paste0("Microbiome KO |", model_name," vs primary model ",tri,"_",grp),
             x="p.value in ori", y="p.value in sensitive")+
        geom_hline(yintercept=0, color="black", linetype="solid")+
        geom_vline(xintercept=0, color="black", linetype="solid")+
        geom_abline(intercept=0, slope=1)+
        annotate("text", x=-Inf, y=Inf, label=stat_label,
                 hjust=-0.05, vjust=1.2, size=3.5, color="black")+
        theme(panel.background=element_blank(),plot.background=element_blank())
      #ggsave(filename = paste0("BMI_sensitive_",tri,"_",grp,".pdf",sep=""), plot = p[[paste(tri, grp, sep="_")]], width = 6, height = 5)
        
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
    
    print(dim(m_merge))

    p_all<-ggplot(m_merge, aes(x=log2MedianRatio.x, y=log2MedianRatio.y,color=group))+
        geom_point(alpha=0.4)+
        geom_smooth(method="lm", se=FALSE,alpha=0.6)+
        labs(title=paste0("Microbiome KO | ",model_name," vs primary model: ",tri),
            x="log2MedianRatio in primary model", y=paste0("log2MedianRatio in ",model_name))+
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
    return(p)
}

p_list<-test_sensitive_BMI(res_sensitive,res_ori,model_name="design 3 model")

## plot merged plots for all contrasts in one figure
cov<-"design 3 model"
all_plots<-list(
  p_list$early, p_list$middle, p_list$late, p_list$overall
)
p_tmp <- wrap_plots(all_plots, nrow = 2, ncol = 2) +
  plot_annotation(title = paste0("Differential microbiome KO ",cov, " vs primary model")) &
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
ggsave(p_tmp,filename="microbiome ko gsea复现/sensitive/BMI_microbiome_sensitive_allinone.pdf",width=15,height=10)
  

## BMI
all_plots<-list(p_list$early_wasting,p_list$middle_wasting,p_list$late_wasting,p_list$overall_wasting,
                p_list$early_normal,p_list$middle_normal,p_list$late_normal,p_list$overall_normal,
                p_list$early_overweight,p_list$middle_overweight,p_list$late_overweight,p_list$overall_overweight,
                p_list$early_obese,p_list$middle_obese,p_list$late_obese,p_list$overall_obese)

p_tmp <- wrap_plots(all_plots, nrow = 4, ncol = 4) +
  plot_annotation(title = "Differential microbiome KO design 3 vs primary model") &
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))

ggsave(p_tmp,filename=paste0("microbiome ko gsea复现/sensitive/BMI_microbiome_sensitive_all.pdf",sep=""),width=20,height=15)



## metabolite & obstetric
res_microbiome_metabolite<-data.frame()
for (tri in c("early","middle","late","overall")){
  df<-read.csv(paste("microbiome ko gsea复现/sensitive/TJGDM1962s_ko_tri_95.zinb_metabolite",tri,"adj_rpk.csv",sep="_")) %>%
    mutate(trimester=tri)
  res_microbiome_metabolite<-rbind(res_microbiome_metabolite,df)
}

res_microbiome_obstetric<-data.frame()
for (tri in c("early","middle","late","overall")){
  df<-read.csv(paste("microbiome ko gsea复现/sensitive/TJGDM1962s_ko_tri_95.zinb_obstetric",tri,"adj_rpk.csv",sep="_")) %>%
    mutate(trimester=tri)
  res_microbiome_obstetric<-rbind(res_microbiome_obstetric,df)
}


test_sensitive<-function(res_sensitive,res_ori,model_name="design 4 model"){
  stat_df<-data.frame()
   ## Jaccard with ori for each period
  for (tri in unique(res_sensitive$trimester)){
    res_ori_tri <- res_ori %>% dplyr::filter(trimester==tri)
    res_sensitive_tri <- res_sensitive %>% dplyr::filter(trimester==tri)
    sig_ori_tri <- res_ori_tri$feature[res_ori_tri$p.value<0.05 & abs(res_ori_tri$log2MedianRatio)>0.5]
    sig_sensitive_tri <- res_sensitive_tri$feature[res_sensitive_tri$p.value<0.05 & abs(res_sensitive_tri$log2MedianRatio)>0.5]
    sig_common <- intersect(sig_ori_tri, sig_sensitive_tri)
    sig_share <- union(sig_ori_tri, sig_sensitive_tri)
    jaccard <- if (length(sig_share)==0) NA_real_ else {
      length(sig_common) / length(sig_share)
    }
    print(paste0("Jaccard with ori ",tri,":",jaccard))
  }



  ## Spearman and Kappa between sensitive results and ori for each period 
  p<-list()
  m_merge<-data.frame()
  for (tri in unique(res_sensitive$trimester)){
      #p_tri<-list()
      res_ori_tri<-res_ori %>%
        dplyr::filter(trimester==tri,p.value<0.9) %>%
        dplyr::select(feature, p.value)
      res_sensitive_tri<-res_sensitive %>%
        dplyr::filter(trimester==tri,p.value<0.9) %>%
        dplyr::select(feature, p.value)

      m<-merge(res_ori_tri,res_sensitive_tri,by="feature")
      m<-m[complete.cases(m),]
      m_tmp<-m %>% mutate(trimester=tri)
      m_merge<-rbind(m_merge,m_tmp)

      ## log2fc Spearman with ori; only retain p.value < 0.1 in both datasets
      spearman_valid <- nrow(m) >= 3 &&
        length(unique(m$p.value.x)) >= 2 &&
        length(unique(m$p.value.y)) >= 2
      if (spearman_valid){
        pvalue_spearman <- cor.test(m$p.value.x, m$p.value.y, method="spearman")
        spearman_r <- unname(pvalue_spearman$estimate)
        spearman_p <- pvalue_spearman$p.value
        print(paste0("pvalue spearman r with ori ",tri,":",spearman_r))
        print(paste0("pvalue spearman p with ori ",tri,":",spearman_p))
      } else {
        spearman_r <- NA_real_
        spearman_p <- NA_real_
        print(paste0("pvalue spearman skipped with ori ",tri,": insufficient valid data"))
      }

      ## pvalue direction with ori and kappa; only retain p.value < 0.1 in both datasets
      if (nrow(m)==0){
        print(paste0("pvalue log2fc direction skipped with ori ",tri,": no overlapping valid data"))
        next
      }
      signA <- sign(m$p.value.x)
      signB <- sign(m$p.value.y)
      tab <- table(factor(signA, levels = c(-1, 1), labels = c("-1", "1")),
                   factor(signB, levels = c(-1, 1), labels = c("-1", "1")))

      #N       <- sum(tab)
      #concord <- (tab[1,1] + tab[2,2]) / N
      #flip    <- (tab[1,2] + tab[2,1]) / N

      #a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
      #Po <- (a + d) / N
      #pA_up  <- (c + d) / N; pA_down <- (a + b) / N
      #pB_up  <- (b + d) / N; pB_down <- (a + c) / N
      #Pe <- pA_up*pB_up + pA_down*pB_down
      #kappa <- if (Pe==1) NA_real_ else (Po - Pe) / (1 - Pe)

      #print(paste0("pvalue direction concord with ori ",tri,":",concord))
      #print(paste0("pvalue direction flip with ori ",tri,":",flip))
      #print(paste0("pvalue direction kappa with ori ",tri,":",kappa))

      stat_df_tmp<-data.frame(
        comparison_type="sensitive_vs_ori",
        trimester=tri,
        contrast=paste0("ori_vs_sensitive"),
        spearman_r=spearman_r,
        spearman_p=spearman_p,
        n_overlap=nrow(m),
        stringsAsFactors=FALSE)
      stat_df<-rbind(stat_df,stat_df_tmp)

      stat_label <- paste0(
        #"\nJaccard = ", ifelse(is.na(jaccard), "NA", sprintf("%.3f", jaccard)),
        "\nSpearman r = ", ifelse(is.na(spearman_r), "NA", sprintf("%.3f", spearman_r)),
        "\nSpearman p = ", ifelse(is.na(spearman_p), "NA", format.pval(spearman_p, digits=3, eps=0.001))
        #"\nlog2MedianRatio direction concord = ", ifelse(is.na(concord), "NA", sprintf("%.3f", concord)),
        #"\nKappa = ", ifelse(is.na(kappa), "NA", sprintf("%.3f", kappa))
      )

      p[[paste(tri, sep="_")]] <- ggplot(m, aes(x=p.value.x, y=p.value.y))+
        geom_point()+
        geom_smooth(method="lm", se=FALSE, color="blue")+
        labs(title=paste0("Microbiome KO |", model_name," vs primary model ",tri),
             x="pvalue in ori", y="pvalue in sensitive")+
        geom_hline(yintercept=0, color="black", linetype="solid")+
        geom_vline(xintercept=0, color="black", linetype="solid")+
        annotate("text", x=-Inf, y=Inf, label=stat_label,
                 hjust=-0.05, vjust=1.2, size=3.5, color="black")+
        theme(panel.background=element_blank(),plot.background=element_blank())
      #ggsave(filename = paste0("BMI_sensitive_",tri,"_",grp,".pdf",sep=""), plot = p[[paste(tri, grp, sep="_")]], width = 6, height = 5)
        
  }
  
  lab_list <- list()
    for (tri in unique(stat_df$trimester)) {
      sub<-stat_df %>% dplyr::filter(trimester==tri)
      sr<-sub$spearman_r[1]
      sp<-sub$spearman_p[1]
      lab_list[[tri]] <- data.frame(
        trimester = tri,
        stat_label = paste0(
          tri, "\n",
          "r = ", sprintf("%.3f", unname(sr)), "\n",
          "p = ", format.pval(sp, digits = 3, eps = 0.001)
        ),
        stringsAsFactors = FALSE
      )
    }
  lab_df <- bind_rows(lab_list)

  p_all<-ggplot(m_merge, aes(x=p.value.x, y=p.value.y,color=trimester))+
      geom_point(alpha=0.4)+
      geom_smooth(method="lm", se=FALSE,alpha=0.6)+
      labs(title=paste0("Microbiome KO | ",model_name," vs primary model: "),
           x="pvalue in primary model", y=paste0("pvalue in ",model_name))+
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

  ggsave(p_all,filename=paste0("microbiome ko gsea复现/sensitive/sensitive_",model_name,"_allinone.pdf"), width=7, height=6)
  
      return(p)
}



p_list<-test_sensitive(res_microbiome_metabolite,res_ori,model_name="design 4 model")
cov<-"desgin 4 model"
all_plots<-list(p_list$early,p_list$middle,p_list$late,p_list$overall)

p_tmp <- wrap_plots(all_plots, nrow = 2, ncol = 2) +
  plot_annotation(title = paste0("Differential microbiome KO ",cov," vs primary model")) &
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(p_tmp,filename=paste0("microbiome ko gsea复现/sensitive/sensitive_microbiome_",cov,"_all.pdf",sep=""),width=20,height=15)



## sig KO retain
ori_sig<-res_ori$feature[which(res_ori$p.value<0.05 & abs(res_ori$log2MedianRatio)>0.5)]
res_sensitive_sig<-res_sensitive %>%
  dplyr::filter(p.value<0.05 & abs(log2MedianRatio)>0.5)
for (grp in unique(res_sensitive_sig$BMI_group)){
  sig<-res_sensitive_sig$feature[res_sensitive_sig$BMI_group==grp]
  print(paste0("sensitive ",grp," sig KO number:",length(unique(sig))))
  share_sig<-intersect(ori_sig,sig)
  print(paste0("sensitive ",grp," shared sig KO number:",length(share_sig)))
  new_sig<-setdiff(sig,ori_sig)
  print(paste0("sensitive ",grp," new sig KO number:",length(new_sig)))
}

sig_metabolite<-res_microbiome_metabolite$feature[which(res_microbiome_metabolite$p.value<0.05 & abs(res_microbiome_metabolite$log2MedianRatio)>0.5)]
share_sig<-intersect(ori_sig,sig_metabolite)
new_sig<-setdiff(sig_metabolite,ori_sig)
print(paste0("sensitive metabolite sig KO number:",length(unique(sig_metabolite))))
print(paste0("sensitive metabolite shared sig KO number:",length(share_sig)))
print(paste0("sensitive metabolite new sig KO number:",length(new_sig)))

sig_obstetric<-res_microbiome_obstetric$feature[which(res_microbiome_obstetric$p.value<0.05 & abs(res_microbiome_obstetric$log2MedianRatio)>0.5)]
share_sig<-intersect(ori_sig,sig_obstetric)
new_sig<-setdiff(sig_obstetric,ori_sig)
print(paste0("sensitive obstetric sig KO number:",length(unique(sig_obstetric))))
print(paste0("sensitive obstetric shared sig KO number:",length(share_sig)))
print(paste0("sensitive obstetric new sig KO number:",length(new_sig)))
