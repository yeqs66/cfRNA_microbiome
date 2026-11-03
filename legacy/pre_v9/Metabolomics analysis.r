# Metabolite differential analysis (Wilcoxon) + volcano plot
library(ggplot2)
library(dplyr)
library(ggrepel)

# ------------------------------ 1. Simulate metabolite data ------------------------------
set.seed(456)
n_met <- 50
met_names <- paste0("Met", 1:n_met)
n_samples <- 100

# Simulate abundance for BCAA, SCFA, Bile acids
category <- sample(c("BCAA", "SCFA", "Bile acids"), n_met, replace = TRUE, prob = c(0.2, 0.3, 0.5))

# Create data frame: sample-wise metabolite levels
met_data <- matrix(runif(n_met * n_samples, 0, 10), nrow = n_met, ncol = n_samples,
                   dimnames = list(met_names, paste0("sample", 1:n_samples)))
# Add group labels
coldata <- data.frame(
  sample = colnames(met_data),
  group = sample(c("GDM", "Control"), n_samples, replace = TRUE),
  trimester = sample(c("early", "middle", "late"), n_samples, replace = TRUE),
  stringsAsFactors = FALSE
)

# Introduce some differential effects (e.g., GDM up in 2nd/3rd for some metabolites)
for (i in 1:n_met) {
  if (runif(1) < 0.3) {  # 30% metabolites are differential
    effect <- runif(1, 0.5, 2)
    idx <- which(coldata$group == "GDM" & coldata$trimester %in% c("middle", "late"))
    met_data[i, idx] <- met_data[i, idx] + effect
  }
}

# ------------------------------ 2. Wilcoxon test per trimester + overall ------------------------------
run_wilcox <- function(met_vec, group_vec) {
  test <- wilcox.test(met_vec ~ group_vec)
  return(test$p.value)
}

compute_log2FC <- function(met_vec, group_vec) {
  log2((mean(met_vec[group_vec == "GDM"], na.rm = TRUE) + 1) /
       (mean(met_vec[group_vec == "Control"], na.rm = TRUE) + 1))
}

results <- data.frame()
for (met in met_names) {
  for (tri in c("early", "middle", "late", "all")) {
    if (tri == "all") {
      idx <- 1:n_samples
    } else {
      idx <- which(coldata$trimester == tri)
    }
    if (length(unique(coldata$group[idx])) < 2) next
    pval <- run_wilcox(met_data[met, idx], coldata$group[idx])
    logFC <- compute_log2FC(met_data[met, idx], coldata$group[idx])
    results <- rbind(results, data.frame(
      metabolite = met,
      trimester = tri,
      pvalue = pval,
      log2FC = logFC,
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values per trimester
results <- results %>%
  group_by(trimester) %>%
  mutate(padj = p.adjust(pvalue, "BH"),
         logp = -log10(padj)) %>%
  ungroup()

# Add category info
results <- results %>%
  left_join(data.frame(metabolite = met_names, category = category), by = "metabolite")

write.csv(results, "metabolite_results.csv", row.names = FALSE)

# ------------------------------ 3. Volcano plot (overall trimester) ------------------------------
volc_data <- results %>% filter(trimester == "all")
volc_data <- volc_data %>%
  mutate(sign = ifelse(padj < 0.05 & log2FC > 0.15, "Up-regulated",
                       ifelse(padj < 0.05 & log2FC < -0.15, "Down-regulated", "Non-significant")),
         label = ifelse(padj < 0.05 & abs(log2FC) > 0.3, metabolite, ""))

p_volcano <- ggplot(volc_data, aes(x = log2FC, y = logp, color = sign, shape = category)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_text_repel(aes(label = label), size = 3, color = "black", max.overlaps = 20) +
  scale_color_manual(values = c("Up-regulated" = "#E18E6D",
                                "Down-regulated" = "#62B197",
                                "Non-significant" = "grey80")) +
  scale_shape_manual(values = c("BCAA" = 15, "SCFA" = 16, "Bile acids" = 17)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed") +
  labs(x = "Log2 fold change", y = "-Log10 FDR", title = "Metabolites (Overall)") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_volcano)

# Optionally save
# ggsave("metabolite_volcano.pdf", p_volcano, width = 8, height = 6)

# ------------------------------ 4. Optional: boxplots for top hits ------------------------------
top_hits <- volc_data %>% filter(padj < 0.05) %>% arrange(padj) %>% head(6)
if (nrow(top_hits) > 0) {
  # Prepare boxplot data for top metabolites
  box_data <- coldata
  for (met in top_hits$metabolite) {
    box_data[[met]] <- met_data[met, box_data$sample]
  }
  box_long <- box_data %>%
    pivot_longer(cols = all_of(top_hits$metabolite), names_to = "metabolite", values_to = "abundance")
  
  p_box <- ggplot(box_long, aes(x = trimester, y = abundance, fill = group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~ metabolite, scales = "free") +
    scale_fill_manual(values = c("GDM" = "#E18E6D", "Control" = "#62B197")) +
    labs(y = "Abundance") +
    theme_minimal()
  print(p_box)
}