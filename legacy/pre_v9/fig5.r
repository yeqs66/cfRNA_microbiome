# Figure 5: MOFA multi-omics integration
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(grid)

# -------------------- 5A: Variance explained --------------------
set.seed(20)
views <- c("KO", "mRNA", "miRNA", "lncRNA")
factors <- paste0("Factor", 1:8)
# Simulate variance explained matrix (views x factors)
var_mat <- matrix(runif(length(views) * length(factors), 0, 0.3),
                  nrow = length(views), ncol = length(factors),
                  dimnames = list(views, factors))
# Scale to sum to 1 per view (like real MOFA)
var_mat <- var_mat / rowSums(var_mat) * 0.8  # total variance per view ~0.8

# Convert to long format for ggplot
var_df <- as.data.frame(var_mat) %>%
  rownames_to_column("View") %>%
  pivot_longer(-View, names_to = "Factor", values_to = "Variance")

pA <- ggplot(var_df, aes(x = Factor, y = Variance, fill = View)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Variance explained", title = "Variance explained per factor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")

# -------------------- 5B: Factor-trait correlation heatmap --------------------
set.seed(21)
traits <- paste0("Trait", 1:12)
factors <- paste0("Factor", 1:5)  # only top 5 factors
cor_mat <- matrix(runif(length(traits) * length(factors), -0.4, 0.4),
                  nrow = length(traits), ncol = length(factors),
                  dimnames = list(traits, factors))
# Make some stronger correlations
cor_mat[1:3, 1:2] <- runif(6, 0.3, 0.6)
cor_mat[4:6, 3:4] <- runif(6, -0.5, -0.3)
pval_mat <- matrix(runif(length(traits) * length(factors), 0.001, 0.2),
                   nrow = length(traits), ncol = length(factors))

# Add annotation: trait class (simulated)
trait_class <- rep(c("Anthropometric", "Glycemic", "Lipid", "Liver"), each = 3)[1:length(traits)]
ha_row <- rowAnnotation(Class = trait_class,
                        col = list(Class = c("Anthropometric" = "#E41A1C",
                                             "Glycemic" = "#377EB8",
                                             "Lipid" = "#4DAF4A",
                                             "Liver" = "#984EA3")),
                        show_annotation_name = TRUE)

htB <- Heatmap(cor_mat, name = "Correlation",
               col = colorRamp2(c(-0.5, 0, 0.5), c("#6D3580", "white", "#FFE26F")),
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10),
               right_annotation = ha_row,
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if (pval_mat[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 12))
               },
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10)))

# -------------------- 5C: miRNA-target KEGG enrichment (back-to-back lollipop) --------------------
set.seed(22)
pathways <- paste0("Path", 1:20)
factor1_ratio <- runif(20, 0.05, 0.25)
factor2_ratio <- runif(20, 0.05, 0.25)
log10padj1 <- -log10(runif(20, 0.001, 0.1))
log10padj2 <- -log10(runif(20, 0.001, 0.1))
# Order pathways by combined score
order_idx <- order(factor1_ratio + factor2_ratio, decreasing = TRUE)
pathways <- pathways[order_idx]
factor1_ratio <- factor1_ratio[order_idx]
factor2_ratio <- factor2_ratio[order_idx]
log10padj1 <- log10padj1[order_idx]
log10padj2 <- log10padj2[order_idx]

# Factor 1 (left): negative bars, Factor 2 (right): positive bars
df_left <- data.frame(Pathway = pathways, log10padj = log10padj1, Ratio = factor1_ratio,
                      Factor = "Factor1")
df_right <- data.frame(Pathway = pathways, log10padj = log10padj2, Ratio = factor2_ratio,
                       Factor = "Factor2")

# Use ggplot with geom_segment to create lollipop
# Factor1: plot negative x ( -log10padj ), Factor2: positive x
df_left$x <- -log10padj1
df_right$x <- log10padj2
df_comb <- rbind(df_left, df_right)
df_comb$Factor <- factor(df_comb$Factor, levels = c("Factor2", "Factor1"))  # left first

pC <- ggplot(df_comb, aes(x = x, y = Pathway, color = Factor)) +
  geom_segment(aes(x = 0, xend = x, y = Pathway, yend = Pathway),
               linetype = "dashed", color = "grey50", size = 0.3) +
  geom_point(aes(size = Ratio), shape = 19) +
  scale_x_continuous(breaks = seq(-4, 4, 2), labels = abs(seq(-4, 4, 2))) +
  scale_color_manual(values = c("Factor1" = "#a52a2a", "Factor2" = "#dc968d")) +
  scale_size_continuous(range = c(2, 6), name = "GeneRatio") +
  labs(x = "-log10(p.adj)", y = "", title = "miRNA-target KEGG enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 8),
        legend.position = "bottom",
        panel.grid.major.y = element_blank())

# Combine factor1 and factor2 side by side with negative/positive axes
# Actually the above already does it with negative/positive x.
# To make back-to-back with separate panels we can use facet_grid, but above is simpler.

# If want separate panels, use:
pC_left <- ggplot(df_left, aes(x = -log10padj, y = Pathway)) +
  geom_segment(aes(x = 0, xend = -log10padj, y = Pathway, yend = Pathway),
               color = "#a52a2a", size = 0.5) +
  geom_point(aes(size = Ratio), color = "#a52a2a") +
  scale_x_continuous(limits = c(-4, 0), breaks = c(-4, -2, 0), labels = c(4, 2, 0)) +
  labs(x = "-log10(p.adj)", y = "", title = "Factor 1") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 8))

pC_right <- ggplot(df_right, aes(x = log10padj, y = Pathway)) +
  geom_segment(aes(x = 0, xend = log10padj, y = Pathway, yend = Pathway),
               color = "#dc968d", size = 0.5) +
  geom_point(aes(size = Ratio), color = "#dc968d") +
  scale_x_continuous(limits = c(0, 4), breaks = c(0, 2, 4), labels = c(0, 2, 4)) +
  labs(x = "-log10(p.adj)", y = "", title = "Factor 2") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # hide y labels on right side
        axis.text.x = element_text(size = 8))

# Combine left and right with patchwork
pC_combined <- pC_left + pC_right + plot_layout(widths = c(1, 1))

# -------------------- Output --------------------
print("Figure 5A: Variance explained")
print(pA)

print("Figure 5B: Factor-trait correlation heatmap")
draw(htB)

print("Figure 5C: miRNA-target KEGG enrichment (back-to-back)")
print(pC_combined)