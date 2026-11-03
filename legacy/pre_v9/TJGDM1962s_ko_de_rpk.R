library(tidyverse)
library(dplyr)
library(glmmTMB)
library(emmeans)

# load clinical data
#setwd("/jdfssz1/ST_HEALTH/P17Z10200N0306/P17Z10200N0306/wangkeqing/GDM/AB/AB_Differtial/jdfssz1/ST_HEALTH/P17Z10200N0306/P17Z10200N0306/wangkeqing/GDM/AB/AB_Differtial")
load("/jdfssz1/ST_HEALTH/P17Z10200N0306/P17Z10200N0306/wangkeqing/GDM/GDM_KO_DE/clinical_347_mo_kid_pairs_fromLX_20240424.rda")
# sample info
load("/jdfssz1/ST_HEALTH/P17Z10200N0306/P17Z10200N0306/wangkeqing/GDM/GDM_KO_DE/microbiome_sampleinfo_GDM1999s_20240318.rda")

TJGDM_microbiome_ab_tri_95 <- read.csv("/jdfssz1/ST_HEALTH/P17Z10200N0306/P17Z10200N0306/wangkeqing/GDM/GDM_KO_DE/TJGDM1962s_ko_tri_95_rpk.csv", header = T, row.names = 1)

#if(FALSE){


df <- t(TJGDM_microbiome_ab_tri_95) %>% as.data.frame()
microbiome_sampleinfo_zinb <- clinical_347_mo_kid_pairs %>%
  dplyr::select(ID,BMI_mo,ethnicity_mo,smoking_mo,drinking_mo_18topreg) %>%
  merge(.,microbiome_sampleinfo, by="ID") %>%
  dplyr::arrange(sample) %>%
  dplyr::select(sample,ID,group,trimester,BMI_mo,ethnicity_mo,smoking_mo,drinking_mo_18topreg) %>%
  dplyr::filter(sample %in% row.names(df)) %>%
  mutate(trimester=factor(trimester,levels = c("early","middle","late")),
         ID=as.character(ID),group=factor(group,levels=c("GDM","Control")),
         BMI_mo=as.numeric(BMI_mo),ethnicity_mo=factor(ethnicity_mo),
         smoking_mo=factor(smoking_mo),drinking_mo_18topreg=factor(drinking_mo_18topreg))
indata <- merge(df,microbiome_sampleinfo_zinb,by.x=0,by.y=1)

overall <- data.frame()
early <- data.frame()
middle <- data.frame()
late <- data.frame()

overall.adj <- data.frame()
early.adj <- data.frame()
middle.adj <- data.frame()
late.adj <- data.frame()


for (i in colnames(df)){
#for (i in c("K02424")){

  print(i)  

  overall_fail_out <- data.frame("contrast"="GDM - Control",
                                 "feature"=i,difference=NA,p.value=NA,log2fc=NA,log2MedianRatio=NA)
  early_fail_out <- data.frame("contrast"="GDM early - Control early",
                               "feature"=i,difference=NA,p.value=NA,log2fc=NA,log2MedianRatio=NA)
  middle_fail_out <- data.frame("contrast"="GDM middle - Control middle",
                                "feature"=i,difference=NA,p.value=NA,log2fc=NA,log2MedianRatio=NA)
  late_fail_out <- data.frame("contrast"="GDM late - Control late",
                              "feature"=i,difference=NA,p.value=NA,log2fc=NA,log2MedianRatio=NA)

  early_log2fc <- log2(exp(mean(log(indata[which(indata$trimester=="early" & indata$group=="GDM"),i])))/exp(mean(log(indata[which(indata$trimester=="early" & indata$group=="Control"),i]))))
  middle_log2fc <- log2(exp(mean(log(indata[which(indata$trimester=="middle" & indata$group=="GDM"),i])))/exp(mean(log(indata[which(indata$trimester=="middle" & indata$group=="Control"),i]))))
  late_log2fc <- log2(exp(mean(log(indata[which(indata$trimester=="late" & indata$group=="GDM"),i])))/exp(mean(log(indata[which(indata$trimester=="late" & indata$group=="Control"),i]))))
  all_log2fc <- log2(exp(mean(log(indata[which(indata$group=="GDM"),i])))/exp(mean(log(indata[which(indata$group=="Control"),i]))))


  early_log2mr <- log2(median(indata[which(indata$trimester=="early" & indata$group=="GDM"),i])/median(indata[which(indata$trimester=="early" & indata$group=="Control"),i]))
  middle_log2mr <- log2(median(indata[which(indata$trimester=="middle" & indata$group=="GDM"),i])/median(indata[which(indata$trimester=="middle" & indata$group=="Control"),i]))
  late_log2mr <- log2(median(indata[which(indata$trimester=="late" & indata$group=="GDM"),i])/median(indata[which(indata$trimester=="late" & indata$group=="Control"),i]))
  all_log2mr <- log2(median(indata[which(indata$group=="GDM"),i])/median(indata[which(indata$group=="Control"),i]))



  tryout<-try(zinb.mod.adj <- glmmTMB(reformulate(c("group","trimester","group*trimester",
                                                    "BMI_mo","ethnicity_mo","smoking_mo",
                                                    "drinking_mo_18topreg","(1|ID)"),i), 
                                      ziformula=~., data = indata),silent=TRUE)
  if('try-error' %in% class(tryout)){

    overall.adj <- rbind(overall.adj,overall_fail_out)
    early.adj <- rbind(early.adj, early_fail_out)
    middle.adj <- rbind(middle.adj, middle_fail_out)
    late.adj <- rbind(late.adj, late_fail_out)
    
  }else{
    
    zinb.mod.adj <- glmmTMB(reformulate(c("group","trimester","group*trimester",
                                                    "BMI_mo","ethnicity_mo","smoking_mo",
                                                    "drinking_mo_18topreg","(1|ID)"),i), 
                                      ziformula=~., data = indata)

    tryout2 <- try(zinb.mod.adj <- update(zinb.mod.adj,family=nbinom2),silent=TRUE)
    if('try-error' %in% class(tryout2)){    

      overall.adj <- rbind(overall.adj,overall_fail_out)
      early.adj <- rbind(early.adj, early_fail_out)
      middle.adj <- rbind(middle.adj, middle_fail_out)
      late.adj <- rbind(late.adj, late_fail_out)
    
    }else{

      zinb.mod.adj <- update(zinb.mod.adj,family=nbinom2)

      terms <- c('group', 'trimester')
      lsmeansOut <- emmeans(zinb.mod.adj, as.formula(sprintf('~%s', paste(terms, collapse = '*'))), type="response", component="response")
      lsmeans_output<-summary(lsmeansOut) #%>%
        #dplyr::mutate(emmean=emmean)
      lsmeans_contrast<-pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast,conf.int = T) %>%
        dplyr::mutate(contrast1=contrast) %>%
        tidyr::separate(contrast, c("a","b","c","d")) %>%
        dplyr::filter(b==d) %>%
        dplyr::mutate(difference=estimate, feature=i) %>%
        dplyr::select(contrast1,feature,difference,contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"),contrast=contrast1)

      contrast.group[grep("early", contrast.group$contrast),"log2fc"] <- early_log2fc
      contrast.group[grep("middle", contrast.group$contrast),"log2fc"] <- middle_log2fc
      contrast.group[grep("late", contrast.group$contrast),"log2fc"] <- late_log2fc

      contrast.group[grep("early", contrast.group$contrast),"log2MedianRatio"] <- early_log2mr
      contrast.group[grep("middle", contrast.group$contrast),"log2MedianRatio"] <- middle_log2mr
      contrast.group[grep("late", contrast.group$contrast),"log2MedianRatio"] <- late_log2mr

      early.adj <- rbind(early.adj, contrast.group[grep("early", contrast.group$contrast),])
      middle.adj <- rbind(middle.adj, contrast.group[grep("middle",contrast.group$contrast),])
      late.adj <- rbind(late.adj, contrast.group[grep("late", contrast.group$contrast),])
    
      lsmeansOut <- emmeans(zinb.mod.adj, ~group, type="response", component="response")
      lsmeans_output<-summary(lsmeansOut) #%>%
        #dplyr::mutate(emmean=emmean)
      lsmeans_contrast<-pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast,conf.int = T) %>%
        dplyr::mutate(difference=estimate, feature=i) %>%
        dplyr::select(contrast, feature,difference,contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"))

      contrast.group[1,"log2fc"] <- all_log2fc
      contrast.group[1,"log2MedianRatio"] <- all_log2mr
      overall.adj <- rbind(overall.adj,contrast.group[1,])
    }    
  }
  
  
  
  tryout<-try(zinb.mod <- glmmTMB(reformulate(c("group + trimester + group*trimester"),i), 
                                  ziformula=~., data = indata),silent=TRUE)
  
  if('try-error' %in% class(tryout)){
    
    overall <- rbind(overall,overall_fail_out)
    early <- rbind(early, early_fail_out)
    middle <- rbind(middle, middle_fail_out)
    late <- rbind(late, late_fail_out)
    
  }else{
    
    zinb.mod <- glmmTMB(reformulate(c("group + trimester + group*trimester"),i), 
                        ziformula=~., data = indata)

    tryout2<- try(zinb.mod <- update(zinb.mod,family=nbinom2),silent = TRUE)
    if('try-error' %in% class(tryout2)){

      overall <- rbind(overall,overall_fail_out)
      early <- rbind(early, early_fail_out)
      middle <- rbind(middle, middle_fail_out)
      late <- rbind(late, late_fail_out)

    }else{

      zinb.mod <- update(zinb.mod,family=nbinom2)

      terms <- c('group', 'trimester')
      lsmeansOut <- emmeans(zinb.mod, as.formula(sprintf('~%s', paste(terms, collapse = '*'))), type="response", component="response")
      lsmeans_output<-summary(lsmeansOut) #%>%
        #dplyr::mutate(emmean=emmean)
      lsmeans_contrast<-pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast,conf.int = T) %>%
        dplyr::mutate(contrast1=contrast) %>%
        tidyr::separate(contrast, c("a","b","c","d")) %>%
        dplyr::filter(b==d) %>%
        dplyr::mutate(difference=estimate, feature=i) %>%
        dplyr::select(contrast1,feature,difference,contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"),contrast=contrast1)

      contrast.group[grep("early", contrast.group$contrast),"log2fc"] <- early_log2fc
      contrast.group[grep("middle", contrast.group$contrast),"log2fc"] <- middle_log2fc
      contrast.group[grep("late", contrast.group$contrast),"log2fc"] <- late_log2fc

      contrast.group[grep("early", contrast.group$contrast),"log2MedianRatio"] <- early_log2mr
      contrast.group[grep("middle", contrast.group$contrast),"log2MedianRatio"] <- middle_log2mr
      contrast.group[grep("late", contrast.group$contrast),"log2MedianRatio"] <- late_log2mr
    
      early <- rbind(early, contrast.group[grep("early", contrast.group$contrast),])
      middle <- rbind(middle, contrast.group[grep("middle", contrast.group$contrast),])
      late <- rbind(late, contrast.group[grep("late", contrast.group$contrast),])
    
      lsmeansOut <- emmeans(zinb.mod, ~group, type="response", component="response")
      lsmeans_output<-summary(lsmeansOut) #%>%
        #dplyr::mutate(emmean=emmean)
      lsmeans_contrast<-pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast,conf.int = T) %>%
        dplyr::mutate(difference=estimate, feature=i) %>%
        dplyr::select(contrast, feature,difference,contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"))

      contrast.group[1,"log2fc"] <- all_log2fc  
      contrast.group[1,"log2MedianRatio"] <- all_log2mr
      overall <- rbind(overall,contrast.group[1,])

    }
  }
}

overall$p.adj <- p.adjust(p = overall$p.value, method = "BH")
overall.adj$p.adj <- p.adjust(p = overall.adj$p.value, method = "BH")

early$p.adj <- p.adjust(p = early$p.value, method = "BH")
middle$p.adj <- p.adjust(p = middle$p.value, method = "BH")
late$p.adj <- p.adjust(p = late$p.value, method = "BH")

early.adj$p.adj <- p.adjust(p = early.adj$p.value, method = "BH")
middle.adj$p.adj <- p.adjust(p = middle.adj$p.value, method = "BH")
late.adj$p.adj <- p.adjust(p = late.adj$p.value, method = "BH")


write.csv(overall, file="TJGDM1962s_ko_tri_95.zinb_overall_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(overall.adj, file="TJGDM1962s_ko_tri_95.zinb_overall_adj_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(early, file="TJGDM1962s_ko_tri_95.zinb_early_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(early.adj, file="TJGDM1962s_ko_tri_95.zinb_early_adj_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(middle, file="TJGDM1962s_ko_tri_95.zinb_middle_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(middle.adj, file="TJGDM1962s_ko_tri_95.zinb_middle_adj_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(late, file="TJGDM1962s_ko_tri_95.zinb_late_rpk.csv", quote=FALSE, row.names = FALSE)
write.csv(late.adj, file="TJGDM1962s_ko_tri_95.zinb_late_adj_rpk.csv", quote=FALSE, row.names = FALSE)
  

