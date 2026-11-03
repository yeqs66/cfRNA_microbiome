# Sensitivity analysis: compare alternative models (mock)
library(ggplot2)
library(dplyr)
library(patchwork)

# ------------------------------ 1. Simulate differential results ------------------------------
set.seed(123)
n_features <- 500
features <- paste0("feature", 1:n_features)

# Primary model (design 1: disease + trimester + batch)
res_primary <- data.frame(
  feature = features,
  logFC = rnorm(n_features, 0, 0.8),
  P.Value = runif(n_features, 0, 0.2),
  stringsAsFactors = FALSE
)
res_primary$P.Value[sample(1:n_features, 30)] <- runif(30, 1e-6, 0.05)
res_primary$adj.P.Val <- p.adjust(res_primary$P.Value, "BH")

# Alternative model 1 (design 3: add BMI group)
res_alt1 <- res_primary
res_alt1$logFC <- res_alt1$logFC + rnorm(n_features, 0, 0.2)  # small shift
res_alt1$P.Value <- runif(n_features, 0, 0.2)
res_alt1$P.Value[sample(1:n_features, 25)] <- runif(25, 1e-6, 0.05)
res_alt1$adj.P.Val <- p.adjust(res_alt1$P.Value, "BH")

# Alternative model 2 (design 5: add obstetric covariates)
res_alt2 <- res_primary
res_alt2$logFC <- res_alt2$logFC + rnorm(n_features, 0, 0.3)
res_alt2$P.Value <- runif(n_features, 0, 0.2)
res_alt2$P.Value[sample(1:n_features, 20)] <- runif(20, 1e-6, 0.05)
res_alt2$adj.P.Val <- p.adjust(res_alt2$P.Value, "BH")

# ------------------------------ 2. Sensitivity metrics ------------------------------
compare_models <- function(primary, alt, name) {
  # Jaccard of significant genes (adj.P.Val < 0.1)
  sig_pri <- primary$feature[primary$adj.P.Val < 0.1]
  sig_alt <- alt$feature[alt$adj.P.Val < 0.1]
  jaccard <- length(intersect(sig_pri, sig_alt)) / length(union(sig_pri, sig_alt))
  
  # Spearman correlation of logFC (for features with P<0.1 in both)
  merged <- merge(primary, alt, by = "feature", suffixes = c(".pri", ".alt"))
  merged <- merged[merged$adj.P.Val.pri < 0.1 & merged$adj.P.Val.alt < 0.1, ]
  if (nrow(merged) >= 3) {
    sp <- cor.test(merged$logFC.pri, merged$logFC.alt, method = "spearman")
    spearman_r <- sp$estimate
    spearman_p <- sp$p.value
  } else {
    spearman_r <- NA
    spearman_p <- NA
  }
  
  # Kappa for direction concordance
  if (nrow(merged) > 0) {
    sign_pri <- sign(merged$logFC.pri)
    sign_alt <- sign(merged$logFC.alt)
    tab <- table(sign_pri, sign_alt)
    if (sum(dim(tab)) == 4) {
      N <- sum(tab)
      Po <- (tab[1,1] + tab[2,2]) / N
      p_up <- (tab[2,1] + tab[2,2]) / N
      p_down <- (tab[1,1] + tab[1,2]) / N
      Pe <- p_up^2 + p_down^2
      kappa <- (Po - Pe) / (1 - Pe)
    } else {
      kappa <- NA
    }
  } else {
    kappa <- NA
  }
  
  data.frame(
    model = name,
    jaccard = jaccard,
    spearman_r = spearman_r,
    spearman_p = spearman_p,
    kappa = kappa,
    n_overlap = nrow(merged)
  )
}

stats1 <- compare_models(res_primary, res_alt1, "BMI_model")
stats2 <- compare_models(res_primary, res_alt2, "Obstetric_model")
stats <- rbind(stats1, stats2)
print(stats)

# ------------------------------ 3. Visualization: scatter plots ------------------------------
# Prepare merged data for plotting (only significant in both)
merged1 <- merge(res_primary, res_alt1, by = "feature", suffixes = c(".pri", ".alt"))
merged1 <- merged1[merged1$adj.P.Val.pri < 0.1 & merged1$adj.P.Val.alt < 0.1, ]
merged1$comparison <- "BMI_model"

merged2 <- merge(res_primary, res_alt2, by = "feature", suffixes = c(".pri", ".alt"))
merged2 <- merged2[merged2$adj.P.Val.pri < 0.1 & merged2$adj.P.Val.alt < 0.1, ]
merged2$comparison <- "Obstetric_model"

plot_data <- rbind(merged1, merged2)

p <- ggplot(plot_data, aes(x = logFC.pri, y = logFC.alt, color = comparison)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~ comparison) +
  labs(x = "log2FC (primary model)", y = "log2FC (alternative model)") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)

# ------------------------------ 4. Summary barplot of Jaccard/Kappa ------------------------------
stats_long <- stats %>%
  pivot_longer(cols = c(jaccard, kappa), names_to = "metric", values_to = "value")

p_bar <- ggplot(stats_long, aes(x = model, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Metric value", title = "Sensitivity metrics") +
  theme_minimal()
print(p_bar)

# save combined if needed
# ggsave("sensitivity_plots.pdf", p / p_bar, width = 10, height = 8)