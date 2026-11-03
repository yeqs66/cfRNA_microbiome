## microbiome 敏感性分析建模 —— furrr 多进程加速版
## 用法:
##   前台:   Rscript microbiome敏感性分析建模_furrr.r [BMI_group]
##   后台:   nohup Rscript microbiome敏感性分析建模_furrr.r normal > nohup_zinb.out 2>&1 &
##   看进度: tail -f nohup_zinb.out
## 依赖: tidyverse, glmmTMB, emmeans, furrr, future
## 注意: Linux/macOS 下自动使用 multicore (fork, 共享内存零拷贝); Windows 下用 multisession

library(tidyverse)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(furrr)

## ------------------------------------------------------------------
## 1. 解析命令行参数 (group)
## ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
group <- if (length(args) >= 1 && !is.na(args[1]) && nzchar(trimws(args[1]))) {
  trimws(args[1])
} else {
  "normal"
}
message("BMI_group = ", group)

## ------------------------------------------------------------------
## 2. 读入数据 + 子集筛选 (与原脚本一致)
## ------------------------------------------------------------------
load("__PATH_001__")
load("microbiome ko gsea复现/microbiome_sampleinfo_GDM1999s_20240318.rda")
matrix_meta <- read.csv("microbiome ko gsea复现/TJGDM1962s_ko_tri_95_rpk.csv")

df <- t(matrix_meta) %>% as.data.frame()

microbiome_sampleinfo_zinb <- clinical_347_mo_kid_pairs %>%
  dplyr::select(ID, BMI_mo, ethnicity_mo, smoking_mo, drinking_mo_18topreg) %>%
  merge(., microbiome_sampleinfo, by = "ID") %>%
  dplyr::arrange(sample) %>%
  dplyr::select(sample, ID, group, trimester, BMI_mo, ethnicity_mo,
                smoking_mo, drinking_mo_18topreg) %>%
  dplyr::filter(sample %in% row.names(df)) %>%
  mutate(trimester = factor(trimester, levels = c("early", "middle", "late")),
         ID = as.character(ID),
         group = factor(group, levels = c("GDM", "Control")),
         BMI_mo = as.numeric(BMI_mo),
         ethnicity_mo = factor(ethnicity_mo),
         smoking_mo = factor(smoking_mo),
         drinking_mo_18topreg = factor(drinking_mo_18topreg),
         BMI_group = case_when(
           BMI_mo < 18.5 ~ "wasting",
           BMI_mo >= 18.5 & BMI_mo < 24 ~ "normal",
           BMI_mo >= 24 & BMI_mo < 28 ~ "overweight",
           BMI_mo >= 28 ~ "obese",
           TRUE ~ NA_character_)) %>%
  dplyr::filter(BMI_group == group)

indata <- merge(df, microbiome_sampleinfo_zinb, by.x = 0, by.y = 1)

genes <- colnames(df)

## ------------------------------------------------------------------
## 3. 单基因建模函数 (在 worker 内执行, 需自包含)
##    返回 list(overall, early, middle, late, overall.adj, early.adj,
##               middle.adj, late.adj)
## ------------------------------------------------------------------
fit_one_gene <- function(i, indata) {

  library(glmmTMB); library(emmeans); library(broom)
  library(tidyr); library(dplyr)

  overall_fail_out <- data.frame("contrast" = "GDM - Control",
                                 "feature" = i, difference = NA, p.value = NA,
                                 log2fc = NA, log2MedianRatio = NA)
  early_fail_out  <- data.frame("contrast" = "GDM early - Control early",
                                "feature" = i, difference = NA, p.value = NA,
                                log2fc = NA, log2MedianRatio = NA)
  middle_fail_out <- data.frame("contrast" = "GDM middle - Control middle",
                                "feature" = i, difference = NA, p.value = NA,
                                log2fc = NA, log2MedianRatio = NA)
  late_fail_out   <- data.frame("contrast" = "GDM late - Control late",
                                "feature" = i, difference = NA, p.value = NA,
                                log2fc = NA, log2MedianRatio = NA)

  early_log2fc  <- log2(exp(mean(log(indata[which(indata$trimester == "early"  & indata$group == "GDM"),     i]))) /
                        exp(mean(log(indata[which(indata$trimester == "early"  & indata$group == "Control"), i]))))
  middle_log2fc <- log2(exp(mean(log(indata[which(indata$trimester == "middle" & indata$group == "GDM"),     i]))) /
                        exp(mean(log(indata[which(indata$trimester == "middle" & indata$group == "Control"), i]))))
  late_log2fc   <- log2(exp(mean(log(indata[which(indata$trimester == "late"   & indata$group == "GDM"),     i]))) /
                        exp(mean(log(indata[which(indata$trimester == "late"   & indata$group == "Control"), i]))))
  all_log2fc    <- log2(exp(mean(log(indata[which(indata$group == "GDM"),     i]))) /
                        exp(mean(log(indata[which(indata$group == "Control"), i]))))

  early_log2mr  <- log2(median(indata[which(indata$trimester == "early"  & indata$group == "GDM"),     i]) /
                        median(indata[which(indata$trimester == "early"  & indata$group == "Control"), i]))
  middle_log2mr <- log2(median(indata[which(indata$trimester == "middle" & indata$group == "GDM"),     i]) /
                        median(indata[which(indata$trimester == "middle" & indata$group == "Control"), i]))
  late_log2mr   <- log2(median(indata[which(indata$trimester == "late"   & indata$group == "GDM"),     i]) /
                        median(indata[which(indata$trimester == "late"   & indata$group == "Control"), i]))
  all_log2mr    <- log2(median(indata[which(indata$group == "GDM"),     i]) /
                        median(indata[which(indata$group == "Control"), i]))

  ## ---- 调整模型 (adj) ----
  overall.adj  <- early.adj  <- middle.adj  <- late.adj  <- NULL

  tryout <- try(zinb.mod.adj <- glmmTMB(reformulate(c("group", "trimester", "group*trimester",
                                                      "BMI_mo", "ethnicity_mo", "smoking_mo",
                                                      "drinking_mo_18topreg", "(1|ID)"), i),
                                        ziformula = ~., data = indata), silent = TRUE)
  if ('try-error' %in% class(tryout)) {
    overall.adj <- overall_fail_out
    early.adj   <- early_fail_out
    middle.adj  <- middle_fail_out
    late.adj    <- late_fail_out
  } else {
    zinb.mod.adj <- glmmTMB(reformulate(c("group", "trimester", "group*trimester",
                                          "BMI_mo", "ethnicity_mo", "smoking_mo",
                                          "drinking_mo_18topreg", "(1|ID)"), i),
                            ziformula = ~., data = indata)
    tryout2 <- try(zinb.mod.adj <- update(zinb.mod.adj, family = nbinom2), silent = TRUE)
    if ('try-error' %in% class(tryout2)) {
      overall.adj <- overall_fail_out
      early.adj   <- early_fail_out
      middle.adj  <- middle_fail_out
      late.adj    <- late_fail_out
    } else {
      zinb.mod.adj <- update(zinb.mod.adj, family = nbinom2)

      terms <- c('group', 'trimester')
      lsmeansOut <- emmeans(zinb.mod.adj,
                            as.formula(sprintf('~%s', paste(terms, collapse = '*'))),
                            type = "response", component = "response")
      lsmeans_contrast <- pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast, conf.int = T) %>%
        dplyr::mutate(contrast1 = contrast) %>%
        tidyr::separate(contrast, c("a", "b", "c", "d")) %>%
        dplyr::filter(b == d) %>%
        dplyr::mutate(difference = estimate, feature = i) %>%
        dplyr::select(contrast1, feature, difference, contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"), contrast = contrast1)

      contrast.group[grep("early",  contrast.group$contrast), "log2fc"] <- early_log2fc
      contrast.group[grep("middle", contrast.group$contrast), "log2fc"] <- middle_log2fc
      contrast.group[grep("late",   contrast.group$contrast), "log2fc"] <- late_log2fc

      contrast.group[grep("early",  contrast.group$contrast), "log2MedianRatio"] <- early_log2mr
      contrast.group[grep("middle", contrast.group$contrast), "log2MedianRatio"] <- middle_log2mr
      contrast.group[grep("late",   contrast.group$contrast), "log2MedianRatio"] <- late_log2mr

      early.adj  <- contrast.group[grep("early",  contrast.group$contrast), ]
      middle.adj <- contrast.group[grep("middle", contrast.group$contrast), ]
      late.adj   <- contrast.group[grep("late",   contrast.group$contrast), ]

      lsmeansOut <- emmeans(zinb.mod.adj, ~group, type = "response", component = "response")
      lsmeans_contrast <- pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast, conf.int = T) %>%
        dplyr::mutate(difference = estimate, feature = i) %>%
        dplyr::select(contrast, feature, difference, contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"))

      contrast.group[1, "log2fc"] <- all_log2fc
      contrast.group[1, "log2MedianRatio"] <- all_log2mr
      overall.adj <- contrast.group[1, ]
    }
  }

  ## ---- 非调整模型 (unadj) ----
  overall <- early <- middle <- late <- NULL

  tryout <- try(zinb.mod <- glmmTMB(reformulate(c("group + trimester + group*trimester"), i),
                                    ziformula = ~., data = indata), silent = TRUE)
  if ('try-error' %in% class(tryout)) {
    overall <- overall_fail_out
    early   <- early_fail_out
    middle  <- middle_fail_out
    late    <- late_fail_out
  } else {
    zinb.mod <- glmmTMB(reformulate(c("group + trimester + group*trimester"), i),
                        ziformula = ~., data = indata)
    tryout2 <- try(zinb.mod <- update(zinb.mod, family = nbinom2), silent = TRUE)
    if ('try-error' %in% class(tryout2)) {
      overall <- overall_fail_out
      early   <- early_fail_out
      middle  <- middle_fail_out
      late    <- late_fail_out
    } else {
      zinb.mod <- update(zinb.mod, family = nbinom2)

      terms <- c('group', 'trimester')
      lsmeansOut <- emmeans(zinb.mod,
                            as.formula(sprintf('~%s', paste(terms, collapse = '*'))),
                            type = "response", component = "response")
      lsmeans_contrast <- pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast, conf.int = T) %>%
        dplyr::mutate(contrast1 = contrast) %>%
        tidyr::separate(contrast, c("a", "b", "c", "d")) %>%
        dplyr::filter(b == d) %>%
        dplyr::mutate(difference = estimate, feature = i) %>%
        dplyr::select(contrast1, feature, difference, contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"), contrast = contrast1)

      contrast.group[grep("early",  contrast.group$contrast), "log2fc"] <- early_log2fc
      contrast.group[grep("middle", contrast.group$contrast), "log2fc"] <- middle_log2fc
      contrast.group[grep("late",   contrast.group$contrast), "log2fc"] <- late_log2fc

      contrast.group[grep("early",  contrast.group$contrast), "log2MedianRatio"] <- early_log2mr
      contrast.group[grep("middle", contrast.group$contrast), "log2MedianRatio"] <- middle_log2mr
      contrast.group[grep("late",   contrast.group$contrast), "log2MedianRatio"] <- late_log2mr

      early  <- contrast.group[grep("early",  contrast.group$contrast), ]
      middle <- contrast.group[grep("middle", contrast.group$contrast), ]
      late   <- contrast.group[grep("late",   contrast.group$contrast), ]

      lsmeansOut <- emmeans(zinb.mod, ~group, type = "response", component = "response")
      lsmeans_contrast <- pairs(lsmeansOut)
      contrast.group <- broom::tidy(lsmeans_contrast, conf.int = T) %>%
        dplyr::mutate(difference = estimate, feature = i) %>%
        dplyr::select(contrast, feature, difference, contains("p.value")) %>%
        dplyr::rename(p.value = contains("p.value"))

      contrast.group[1, "log2fc"] <- all_log2fc
      contrast.group[1, "log2MedianRatio"] <- all_log2mr
      overall <- contrast.group[1, ]
    }
  }

  list(overall = overall, early = early, middle = middle, late = late,
       overall.adj = overall.adj, early.adj = early.adj,
       middle.adj = middle.adj, late.adj = late.adj)
}

## ------------------------------------------------------------------
## 4. 并行调度 + 实时进度打印 (nohup 友好)
## ------------------------------------------------------------------
n_workers <- min(parallel::detectCores() - 1, 8)
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = n_workers)
} else {
  # Linux/macOS 用 fork, copy-on-write 共享内存, 不会每个进程都吃 13G
  plan(multicore, workers = n_workers)
}
message("Parallel workers = ", n_workers)

# 把基因分成 chunk(每块 = workers*10), 主进程跑完一块就 print 进度
chunk_size <- n_workers * 10
chunks     <- split(genes, ceiling(seq_along(genes) / chunk_size))

t0 <- Sys.time()
results <- vector("list", length(chunks))

for (k in seq_along(chunks)) {
  results[[k]] <- future_map(
    chunks[[k]],
    function(g) fit_one_gene(g, indata),
    .options = furrr_options(seed = 123, globals = c(indata = indata))
  )
  done      <- min(k * chunk_size, length(genes))
  pct       <- 100 * done / length(genes)
  elapsed   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  rate      <- done / max(elapsed, 1e-6)
  eta_mins  <- (length(genes) - done) / max(rate, 1e-6)

  msg <- sprintf("[%.1f%%] %d/%d genes done | %.1f mins elapsed | ETA %.1f mins\n",
                 pct, done, length(genes), elapsed, eta_mins)
  cat(msg)   # 实时输出到 stdout (nohup.out 或终端都能看到)
  message(msg)
}

results <- unlist(results, recursive = FALSE)

## ------------------------------------------------------------------
## 5. 合并结果
## ------------------------------------------------------------------
bind_part <- function(name) bind_rows(map(results, name))

overall     <- bind_part("overall")
early       <- bind_part("early")
middle      <- bind_part("middle")
late        <- bind_part("late")
overall.adj <- bind_part("overall.adj")
early.adj   <- bind_part("early.adj")
middle.adj  <- bind_part("middle.adj")
late.adj    <- bind_part("late.adj")

## ------------------------------------------------------------------
## 6. 多重检验校正 + 写出
## ------------------------------------------------------------------
overall$p.adj     <- p.adjust(p = overall$p.value, method = "BH")
overall.adj$p.adj <- p.adjust(p = overall.adj$p.value, method = "BH")
early$p.adj       <- p.adjust(p = early$p.value, method = "BH")
middle$p.adj      <- p.adjust(p = middle$p.value, method = "BH")
late$p.adj        <- p.adjust(p = late$p.value, method = "BH")
early.adj$p.adj   <- p.adjust(p = early.adj$p.value, method = "BH")
middle.adj$p.adj  <- p.adjust(p = middle.adj$p.value, method = "BH")
late.adj$p.adj    <- p.adjust(p = late.adj$p.value, method = "BH")

write.csv(overall,     file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_overall_rpk.csv"),         quote = FALSE, row.names = FALSE)
write.csv(overall.adj, file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_overall_adj_rpk.csv"),     quote = FALSE, row.names = FALSE)
write.csv(early,       file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_early_rpk.csv"),           quote = FALSE, row.names = FALSE)
write.csv(early.adj,   file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_early_adj_rpk.csv"),       quote = FALSE, row.names = FALSE)
write.csv(middle,      file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_middle_rpk.csv"),          quote = FALSE, row.names = FALSE)
write.csv(middle.adj,  file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_middle_adj_rpk.csv"),      quote = FALSE, row.names = FALSE)
write.csv(late,        file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_late_rpk.csv"),            quote = FALSE, row.names = FALSE)
write.csv(late.adj,    file = paste("TJGDM1962s_ko_tri_95.zinb", group, "_late_adj_rpk.csv"),        quote = FALSE, row.names = FALSE)

message("Done. Output files written for group = ", group)
