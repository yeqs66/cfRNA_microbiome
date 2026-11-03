# ============================================================
# 代谢组差异分析：线性混合模型（LMM）双模型
#   primary model : feature ~ disease * trimester + (1|ID)
#   adjusted model: primary + batch + BMI_mo_Anthropometric +
#                   ethnicity_mo + smoking_mo + drinking_mo_18topreg
# 数据流/靶代谢物/样本与《代谢组数据wilcoxon.r》完全一致，结果可直接并排比较
# 输出：CSV 结果表 + primary/adjusted 各一张 all 期火山图 PDF
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("__PATH_001__")

# ---------- 1) 数据加载与 349 样本配对（复用 wilcoxon 数据流，扩展带入协变量） ----------
load("代谢组/metabolics_matrix.rda_终极版.rda")
load("代谢组/mergeinfo_终极版.rda")
load("代谢组/metabolites_ID_addclass_addp_addKEGG_终极版.rda")

# adjusted model 额外协变量
adjusted_cov <- c("BMI_mo_Anthropometric","ethnicity_mo","smoking_mo","drinking_mo_18topreg")

coldata_349 <- read.csv("__PATH_002__",row.names = 1) %>%
  dplyr::mutate(sample=gsub("S","",sample),
                ID=as.character(ID),
                trimester=case_when(trimester==1 ~ "early",
                                    trimester==2 ~ "middle",
                                    trimester==3 ~ "late"),
                sign=paste(ID,trimester,sep = "_"),
                batch=as.factor(batch),
                BMI_mo_Anthropometric=as.numeric(BMI_mo_Anthropometric),
                ethnicity_mo=as.factor(ethnicity_mo),
                smoking_mo=as.factor(smoking_mo),
                drinking_mo_18topreg=as.factor(drinking_mo_18topreg)) %>%
  dplyr::rename(sample_cfRNA=sample) %>%
  dplyr::select(-batch) ## 不需要cfRNA的batch变量

# 匹配 349 例 cfRNA 样本，并一次性带入 5 个协变量
metabolics_matrix <- metabolics_matrix %>%
  dplyr::mutate(sign=paste(ID,trimester,sep = "_"),
                ID=as.character(ID)) %>%
  left_join(coldata_349[,c("sign","sample_cfRNA",adjusted_cov)],by="sign") %>%
  dplyr::filter(!is.na(sample_cfRNA))

adjusted_cov <- c("batch","BMI_mo_Anthropometric","ethnicity_mo","smoking_mo","drinking_mo_18topreg")


# ---------- 2) 靶代谢物筛选（BCAA + SCFA，与 wilcoxon 完全一致） ----------
select_meta <- read.csv("代谢组/target_metabolites_三分类.csv")
select_meta_SFCA <- read.csv("代谢组/select_meta_SCFA.csv")
bcaa_id <- mergeinfo$Vname[which(mergeinfo$Compound %in% select_meta$BCAA)]
scfa_id <- mergeinfo$Vname[which(mergeinfo$Vname %in% select_meta_SFCA$Vname)]

trans_meta_name <- mergeinfo %>%
  dplyr::filter(Vname %in% c(bcaa_id,scfa_id)) %>%
  dplyr::select(Vname,Compound,类别) %>%
  dplyr::rename(category=类别) %>%
  mutate(category=case_when(category=="Fatty acids" ~ "SCFA",
                            category=="Amino acids, peptides, and analogues" ~ "BCAA"))

# 建模数据：代谢物值 + 分组 + 受试者 ID（随机效应）+ 协变量
# disease 因子水平 GDM 在前，emmeans 对比方向为 "GDM - Control"
metabolics_matrix_val <- metabolics_matrix %>%
  dplyr::select(sample_cfRNA,ID,trimester,group,all_of(adjusted_cov),c(bcaa_id,scfa_id)) %>%
  dplyr::rename(!!!setNames(trans_meta_name$Vname,trans_meta_name$Compound)) %>%
  dplyr::mutate(disease=factor(group,levels=c("GDM","Control")),
                trimester=factor(trimester,levels=c("early","middle","late")))

# 描述性 log2FC（沿用 wilcoxon 公式：按孕期子集计算，与模型估计值分开报告）
log2fc_fun <- function(feature,tri){
  d <- if (tri=="all") metabolics_matrix_val else metabolics_matrix_val[metabolics_matrix_val$trimester==tri,]
  log2((mean(d[[feature]][d$disease=="GDM"],na.rm=TRUE)+1)/
       (mean(d[[feature]][d$disease=="Control"],na.rm=TRUE)+1))
}

# ---------- 3) 双模型循环拟合 ----------
# Compound 名含连字符/空格（如 "L-Leucine"），公式中必须加反引号，避免被解析为减法
form_primary  <- function(f) as.formula(paste0("`",f,"` ~ disease * trimester + (1|ID)"))
form_adjusted <- function(f) as.formula(paste0("`",f,"` ~ disease * trimester + ",
                                               paste(adjusted_cov,collapse=" + ")," + (1|ID)"))

# 单代谢物 × 单模型：emmeans 一次提取 early/middle/late（分孕期）+ all（三期等权平均）四层对比
# 拟合失败（error，如协变量水平退化）→ 四行 NA 兜底；奇异拟合仅为 warning，不阻断
fit_one <- function(feature,model){
  na_row <- function(tri) data.frame(model=model,Compound=feature,trimester=tri,
                                     contrast="GDM - Control",estimate=NA_real_,pvalue=NA_real_,
                                     log2FC=log2fc_fun(feature,tri),n=NA_integer_)
  tryCatch({
    if (model=="primary") {
      d_fit <- metabolics_matrix_val[!is.na(metabolics_matrix_val[[feature]]),]
      fit <- suppressWarnings(lmer(form_primary(feature),data=d_fit))
    } else {
      d_fit <- metabolics_matrix_val[stats::complete.cases(metabolics_matrix_val[,c(feature,adjusted_cov)]),]
      fit <- suppressWarnings(lmer(form_adjusted(feature),data=d_fit))
    }
    n_fit <- nrow(d_fit)
    # 分孕期对比
    emm_t <- as.data.frame(emmeans(fit,pairwise ~ disease | trimester)$contrasts) %>%
      transmute(model=model,Compound=feature,trimester=as.character(trimester),
                contrast="GDM - Control",estimate,pvalue=p.value,n=n_fit)
    # 总体对比（三个孕期等权平均）
    emm_a <- as.data.frame(emmeans(fit,pairwise ~ disease)$contrasts) %>%
      transmute(model=model,Compound=feature,trimester="all",
                contrast="GDM - Control",estimate,pvalue=p.value,n=n_fit)
    res <- bind_rows(emm_t,emm_a)
    res$log2FC <- sapply(res$trimester,function(tri) log2fc_fun(feature,tri))
    res
  },error=function(e){
    message(sprintf("拟合失败: %s [%s] - %s",feature,model,conditionMessage(e)))
    bind_rows(na_row("early"),na_row("middle"),na_row("late"),na_row("all"))
  })
}

res_lmm_val <- bind_rows(lapply(c("primary","adjusted"),function(m)
  bind_rows(lapply(trans_meta_name$Compound,fit_one,model=m))))

# ---------- 4) 结果整合：按 model × trimester 分组 BH 校正，输出 CSV ----------
res_lmm_val <- res_lmm_val %>%
  group_by(model,trimester) %>%
  mutate(padj=p.adjust(pvalue,method="fdr"),
         logp=-log10(padj)) %>%
  ungroup() %>%
  left_join(trans_meta_name[,c("Compound","category")],by="Compound") %>%
  dplyr::select(model,Compound,category,trimester,contrast,estimate,pvalue,padj,log2FC,logp,n)

write.csv(res_lmm_val,
          file="__PATH_003__",
          row.names=F)

cat("完成：共",nrow(res_lmm_val),"行结果（",n_distinct(res_lmm_val$Compound),"个代谢物 × 2 模型 × 4 时期 ）\n")

# ---------- 5) 火山图（复用 wilcoxon 样式，primary/adjusted 各一张 all 期） ----------
plot_volcano <- function(data_all,model_name,out_pdf){
  data <- data_all %>%
    dplyr::filter(model==model_name,trimester=="all") %>%
    mutate(sign="up",label="")
  data$sign[which(data$log2FC<0)]<-"down"
  data$group<-"Non-significant"
  data$group[which(data$log2FC>0)]<-"Up-regulated"
  data$group[which(data$log2FC<0)]<-"Down-regulated"
  upCompound<-head(data$Compound[which(data$log2FC>0 & data$padj<0.05)],10)
  downCompound<-head(data$Compound[which(data$log2FC<0 & data$padj<0.05)],10)
  dec_sign<-c(as.character(upCompound),as.character(downCompound))
  data$label[match(dec_sign,data$Compound)]<-dec_sign
  up_data<-subset(data,data$Compound %in% upCompound)
  down_data<-subset(data,data$Compound %in% downCompound)

  p<-ggplot(data=data,aes(x=log2FC,y=logp,shape=category,color=group))+
    geom_point(alpha=0.6,size=4.5)+
    geom_text_repel(data=up_data,
                    aes(x=log2FC,y=logp,label=Compound),
                    seed=233,size=3.5,color='black',min.segment.length=0,
                    force=3,force_pull=2,box.padding=0.4,max.overlaps=Inf,
                    segment.linetype=3,segment.color='black',segment.alpha=0.5,
                    fontface="italic",nudge_x=0.5-up_data$log2FC,
                    direction="y",hjust=1)+
    geom_text_repel(data=down_data,
                    aes(x=log2FC,y=logp,label=Compound),
                    seed=233,size=3.5,color='black',min.segment.length=0,
                    force=4,force_pull=2,box.padding=0.4,max.overlaps=Inf,
                    segment.linetype=3,segment.color='black',segment.alpha=0.5,
                    fontface="italic",nudge_x=-0.5-down_data$log2FC,
                    direction="y",hjust=1)+
    scale_shape_manual(values=c("Bile acids"=15,"SCFA"=16,"BCAA"=17))+
    scale_color_manual(values=c("Down-regulated"='#62B197',"Up-regulated"='#E18E6D',"Non-significant"="grey80"))+
    xlab("Log2 fold change")+ylab("-Log10 FDR")+
    theme_classic()+
    ggtitle(paste0("Differential metabolites in overall trimesters (GDM vs Control) - ",model_name," model"))+
    geom_hline(yintercept=1,linetype="dashed")+
    geom_vline(xintercept=c(-0.15,0.15),linetype="dashed")+
    theme(axis.text=element_text(size=14,face="bold",color="black"),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=14,face="bold"),
          legend.title=element_blank(),
          plot.title=element_text(size=14,face="bold",hjust=0.5))
  ggsave(p,filename=out_pdf,width=8,height=6)
}

out_dir <- "__PATH_004__"
plot_volcano(res_lmm_val,"primary",paste0(out_dir,"代谢物_lmm火山图_primary_1019.pdf"))
plot_volcano(res_lmm_val,"adjusted",paste0(out_dir,"代谢物_lmm火山图_adjusted_1019.pdf"))
cat("火山图已输出至 cell reports 目录\n")
