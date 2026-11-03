# Figure 3: Cell origin + WGCNA (mock data)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -------------------- A-B: Donut plots --------------------
set.seed(1)
cell_types <- c("Hepatocyte", "Pancreas", "Muscle", "Immune", "Other")
prop_gdm <- c(0.40, 0.25, 0.15, 0.12, 0.08)
prop_ctrl <- c(0.35, 0.20, 0.18, 0.15, 0.12)
donut_data <- rbind(
  data.frame(Group = "GDM", Cell = cell_types, Prop = prop_gdm),
  data.frame(Group = "Control", Cell = cell_types, Prop = prop_ctrl)
)
make_donut <- function(dat, grp) {
  d <- dat[dat$Group == grp, ]
  d$fraction <- d$Prop / sum(d$Prop)
  d$ymax <- cumsum(d$fraction)
  d$ymin <- c(0, head(d$ymax, -1))
  d$label_pos <- (d$ymax + d$ymin) / 2
  ggplot(d, aes(fill = Cell, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect() + coord_polar(theta = "y") + xlim(c(2, 4)) +
    scale_fill_brewer(palette = "Set3") +
    theme_void() + ggtitle(grp) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
}
p_donut <- make_donut(donut_data, "GDM") + make_donut(donut_data, "Control") +
  plot_layout(guides = "collect")

# -------------------- C-H: Trophoblast boxplots --------------------
set.seed(2)
tropho_cells <- c("early_SCT", "late_SCT", "early_EVT", "late_EVT", "early_VCT", "late_VCT")
samples <- 1:100
tropho_data <- expand.grid(Sample = samples, Trimester = c("1st", "2nd", "3rd"),
                           Cell = tropho_cells, stringsAsFactors = FALSE)
tropho_data$Prop <- runif(nrow(tropho_data), 0, 0.05)
# add group effect for demo
tropho_data$Group <- ifelse(tropho_data$Trimester == "2nd" & tropho_data$Cell == "early_EVT",
                            "GDM", "Control")  # dummy for coloring
tropho_data$Group <- sample(c("GDM", "Control"), nrow(tropho_data), replace = TRUE)

p_box <- ggplot(tropho_data, aes(x = Trimester, y = Prop, fill = Trimester)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ Cell, nrow = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Pastel2") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 8)) +
  labs(x = "Trimester", y = "Cell Proportion")

# -------------------- I: Module pathway enrichment heatmap --------------------
set.seed(3)
modules <- paste0("M", 1:10)
pathways <- paste0("Path", 1:20)
mat_I <- matrix(runif(200, 0, 12), nrow = 20, ncol = 10,
                dimnames = list(pathways, modules))
ht_I <- Heatmap(mat_I, name = "-log10(FDR)", col = colorRamp2(c(0, 6, 12), c("white", "#6BAED6", "#08519C")),
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10)))

# -------------------- J-K: Module gene count + Module-trait correlation --------------------
set.seed(4)
n_mod <- length(modules)
gene_counts <- sample(50:400, n_mod)
deg_p <- runif(n_mod, 0.001, 0.2)
deg_sym <- ifelse(deg_p < 0.001, "***", ifelse(deg_p < 0.01, "**", ifelse(deg_p < 0.05, "*", "")))

traits <- paste0("Trait", 1:6)
mat_K <- matrix(runif(n_mod * length(traits), -0.5, 0.5), nrow = length(traits), ncol = n_mod,
                dimnames = list(traits, modules))
pval_K <- matrix(runif(n_mod * length(traits), 0, 0.1), nrow = length(traits), ncol = n_mod)

# bottom annotation: gene counts + DEG significance
ha_bottom <- HeatmapAnnotation(
  GeneNum = anno_barplot(gene_counts, gp = gpar(fill = "#74C476", col = NA),
                         axis_param = list(gp = gpar(fontsize = 8))),
  DEG = anno_simple(-log10(deg_p), col = colorRamp2(c(0, 1, 2), c("#EFEDF5", "#DADAEB", "#9E9AC8")),
                    pch = deg_sym, pt_size = unit(0.4, "npc"),
                    border = TRUE),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 8)
)

ht_JK <- Heatmap(mat_K, name = "Correlation",
                 col = colorRamp2(c(-0.5, 0, 0.5), c("#6D3580", "white", "#FFE26F")),
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                 bottom_annotation = ha_bottom,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if (pval_K[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 12))
                 },
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 10)))

# -------------------- Output --------------------
print("Figure 3A-3H (Donuts + Boxplots):")
print(p_donut / p_box + plot_layout(heights = c(1, 2)))

print("Figure 3I (Pathway enrichment):")
draw(ht_I)

print("Figure 3J-3K (Module gene count + trait correlation):")
draw(ht_JK)