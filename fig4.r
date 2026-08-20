# Figure 4: Microbiome species + KO GSEA + pathway boxplots
library(ggplot2)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(grid)

# -------------------- 4A: Species dotplot --------------------
set.seed(10)
species <- paste0("s__", c("Bacteroides_vulgatus", "Prevotella_copri", "Faecalibacterium_prausnitzii",
                           "Ruminococcus_bromii", "Akkermansia_muciniphila", "Bifidobacterium_longum",
                           "Eubacterium_rectale", "Roseburia_intestinalis", "Alistipes_putredinis",
                           "Parabacteroides_distasonis"))
trimesters <- c("early", "middle", "late")
dot_data <- expand.grid(Species = species, Trimester = trimesters, stringsAsFactors = FALSE)
dot_data$log10padj <- runif(nrow(dot_data), 0.5, 3)
dot_data$direction <- sample(c("enriched in GDM", "enriched in control"), nrow(dot_data), replace = TRUE)
dot_data$size <- dot_data$log10padj * 1.5

pA <- ggplot(dot_data, aes(x = Species, y = Trimester, size = size, fill = direction)) +
  geom_point(shape = 21, stroke = 0.3) +
  scale_size_continuous(range = c(2, 8), name = "-log10(FDR)") +
  scale_fill_manual(values = c("enriched in GDM" = "#E18E6D", "enriched in control" = "#62B197"),
                    name = "Direction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_blank(),
        legend.position = "bottom") +
  labs(title = "Species-level differences")

# -------------------- 4B: GSEA circular barplot (fan) --------------------
# Simulate GSEA results for KO pathways across trimesters
set.seed(11)
ko_pathways <- paste0("KO", 1:20, ":", c("Valine degradation", "Fatty acid degradation", "Butanoate metabolism",
                                         "Glycolysis", "TCA cycle", "Oxidative phosphorylation",
                                         "Pyruvate metabolism", "Propanoate metabolism", "Cysteine metabolism",
                                         "Arginine biosynthesis", "Lysine degradation", "Phenylalanine metabolism",
                                         "Tryptophan metabolism", "Histidine metabolism", "Tyrosine metabolism",
                                         "Beta-alanine metabolism", "Glycerolipid metabolism", "Glycerophospholipid metabolism",
                                         "Arachidonic acid metabolism", "Linoleic acid metabolism"))
periods <- c("Early", "Middle", "Late", "All")
# Simulate NES and pvalue
fan_data <- expand.grid(Pathway = ko_pathways, Period = periods, stringsAsFactors = FALSE)
fan_data <- fan_data %>%
  group_by(Period) %>%
  mutate(NES = rnorm(n(), 0, 1.5),
         pvalue = runif(n(), 0.001, 0.2),
         pvalue = ifelse(runif(n()) < 0.4, runif(n(), 0.001, 0.05), pvalue),
         log10padj = -log10(pvalue),
         # leading edge ratio (simulated)
         ratio = ifelse(pvalue < 0.1, runif(n(), 0.1, 0.8), 0),
         ratio = ifelse(ratio < 0, 0, ratio)) %>%
  ungroup() %>%
  filter(pvalue < 0.1)  # keep only significant for display

# Prepare for circlize: we need a matrix of ratios with pathways as rows and periods as columns
fan_mat <- fan_data %>%
  select(Pathway, Period, ratio, NES) %>%
  pivot_wider(names_from = Period, values_from = c(ratio, NES), values_fill = list(ratio = 0, NES = 0))
fan_mat <- as.data.frame(fan_mat)
rownames(fan_mat) <- fan_mat$Pathway
fan_mat <- fan_mat[, -1]

# Order pathways by overall ratio
order_path <- order(rowMeans(fan_mat[, grep("ratio", colnames(fan_mat))], na.rm = TRUE), decreasing = TRUE)
fan_mat <- fan_mat[order_path, , drop = FALSE]
path_names <- rownames(fan_mat)

# Prepare for circlize: each period is a track; we need a matrix of ratios and NES colors
ratio_cols <- grep("ratio", colnames(fan_mat))
nes_cols <- grep("NES", colnames(fan_mat))
ratio_mat <- as.matrix(fan_mat[, ratio_cols])
nes_mat <- as.matrix(fan_mat[, nes_cols])
colnames(ratio_mat) <- gsub("ratio_", "", colnames(ratio_mat))
colnames(nes_mat) <- gsub("NES_", "", colnames(nes_mat))
# Ensure periods order
period_order <- c("Early", "Middle", "Late", "All")
ratio_mat <- ratio_mat[, period_order, drop = FALSE]
nes_mat <- nes_mat[, period_order, drop = FALSE]

# Define colors for NES (red positive, green negative)
nes_color <- function(x) {
  if (x > 0) "#E18E6D" else "#62B197"
}
color_mat <- matrix(apply(nes_mat, 2, function(col) sapply(col, nes_color)),
                    nrow = nrow(nes_mat), ncol = ncol(nes_mat),
                    dimnames = dimnames(nes_mat))

# Draw circular plot using circlize
circos_plot <- function() {
  circos.par(start.degree = 270, gap.degree = 180,
             track.height = 0.18, circle.margin = c(0.1, 0.1))
  circos.initialize(letters[1], xlim = c(0, nrow(ratio_mat)))
  # Tracks from outside to inside: All, Late, Middle, Early
  track_order <- c("All", "Late", "Middle", "Early")
  for (tr in track_order) {
    idx <- which(colnames(ratio_mat) == tr)
    values <- ratio_mat[, idx]
    cols <- color_mat[, idx]
    circos.track(ylim = c(0, max(values, na.rm = TRUE) * 1.2),
                 panel.fun = function(x, y) {
                   circos.barplot(values, 1:nrow(ratio_mat) - 0.5,
                                  col = cols, border = NA)
                   # Add axis and label
                   circos.yaxis(side = "left", at = c(0, max(values)*0.6, max(values)*1.2),
                                labels.cex = 0.4)
                   circos.text(nrow(ratio_mat) + 2, 0, tr,
                               facing = "reverse.clockwise", adj = c(0.5, -0.5), cex = 0.6)
                 }, bg.border = NA)
  }
  circos.clear()
}

# Save fan plot as PDF (since circlize doesn't integrate with ggplot easily)
pdf("Figure4B_fan.pdf", width = 6, height = 6)
circos_plot()
dev.off()
# To view, we'll open the saved file (or show inline in RStudio)
# Since we are in a script, we'll just print a message.

# -------------------- 4C-D: Pathway difference boxplots --------------------
# Simulate abundance data for two pathways (e.g., Valine degradation and Fatty acid degradation)
set.seed(12)
pathways_of_interest <- c("Valine degradation (KO00280)", "Fatty acid degradation (KO00071)")
n_samples <- 100
box_data <- expand.grid(Sample = 1:n_samples, Trimester = c("1st", "2nd", "3rd"),
                        Pathway = pathways_of_interest, stringsAsFactors = FALSE)
box_data$Abundance <- runif(nrow(box_data), 0, 10)
# Add GDM effect: increase in GDM for both pathways in 2nd and 3rd trimester
box_data$Group <- sample(c("GDM", "Control"), nrow(box_data), replace = TRUE)
box_data$Abundance <- box_data$Abundance +
  ifelse(box_data$Group == "GDM" & box_data$Trimester %in% c("2nd", "3rd"), runif(nrow(box_data), 1, 3), 0)

pC <- ggplot(box_data %>% filter(Pathway == pathways_of_interest[1]),
             aes(x = Trimester, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = c("GDM" = "#E18E6D", "Control" = "#62B197")) +
  labs(y = "Relative abundance", title = pathways_of_interest[1]) +
  theme_minimal() +
  theme(legend.position = "bottom")

pD <- ggplot(box_data %>% filter(Pathway == pathways_of_interest[2]),
             aes(x = Trimester, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = c("GDM" = "#E18E6D", "Control" = "#62B197")) +
  labs(y = "Relative abundance", title = pathways_of_interest[2]) +
  theme_minimal() +
  theme(legend.position = "bottom")

# -------------------- Output --------------------
print("Figure 4A (Species dotplot):")
print(pA)

print("Figure 4B (Fan plot) saved as Figure4B_fan.pdf in working directory. Open to view.")

print("Figure 4C & 4D (Pathway boxplots):")
print(pC)
print(pD)

# If you want to combine 4C and 4D side by side:
library(patchwork)
pCD <- pC + pD + plot_layout(guides = "collect")
print("Combined C-D:")
print(pCD)