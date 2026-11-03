# Standalone legends for the Cytoscape networks created on 2026-10-26.
# No source data are read. Required packages: ggplot2 and cowplot.
# Usage in R: source("network_legend.R")
# Usage in a terminal: Rscript network_legend.R optional_output_directory

needed <- c("ggplot2", "cowplot")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}
library(ggplot2)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else file.path(getwd(), "network_legend_R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# These are mapping anchors, not categorical cutoffs or observed maxima.
omics <- c("Microbiome KO", "cfRNA", "Metabolite")
node_colors <- setNames(c("#009E73", "#0072B2", "#E69F00"), omics)
node_shapes <- setNames(c(22, 21, 23), omics)  # square, circle, diamond
edge_colors <- c(Positive = "#D77B58", Negative = "#8064A2")
edge_types <- c(Positive = "solid", Negative = "longdash")
edge_alpha <- 155 / 255

# Preserve the Cytoscape diameter ratio 28:44:60 and width ratio 1.3:2.8:4.3.
# R symbol sizes and Cytoscape coordinates use different units.
# These scale factors control the printed size of the standalone legend.
node_size_factor <- 0.125
edge_width_factor <- 0.30
degree_sizes <- setNames(c(28, 44, 60) * node_size_factor, c("1", "7", "13"))
rho_widths <- setNames(c(1.3, 2.8, 4.3) * edge_width_factor,
                      c("0.10", "0.18", "0.26"))

legend_theme <- theme_void(base_family = "sans", base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(10, "mm"),
    legend.key.width = grid::unit(17, "mm"),
    legend.margin = margin(1, 2, 1, 2, unit = "mm"),
    legend.background = element_blank(),
    plot.background = element_rect(fill = NA, colour = NA)
  )

# 1. Merge node color and shape into one guide.
node_data <- data.frame(x = 1, y = 1:3, omics = factor(omics, levels = omics))
p_omics <- ggplot(node_data, aes(x, y, shape = omics, fill = omics)) +
  geom_point(size = 4.6, colour = "#444444", stroke = 0.35) +
  scale_shape_manual(name = "Omics", values = node_shapes) +
  scale_fill_manual(name = "Omics", values = node_colors) +
  legend_theme

# 2. Neutral circles explain node degree without assigning it to one omics.
degree_data <- data.frame(x = 1, y = 1:3,
                          degree = factor(names(degree_sizes), levels = names(degree_sizes)))
p_degree <- ggplot(degree_data, aes(x, y, size = degree)) +
  geom_point(shape = 16, colour = "#777777") +
  scale_size_manual(name = "Node degree", values = degree_sizes) +
  legend_theme

# 3. Merge correlation color and line type into one guide.
direction_data <- data.frame(x = 0, xend = 1, y = 1:2,
                              direction = factor(names(edge_colors), levels = names(edge_colors)))
p_direction <- ggplot(direction_data,
                      aes(x = x, xend = xend, y = y, yend = y,
                          colour = direction, linetype = direction)) +
  geom_segment(linewidth = 0.9, alpha = edge_alpha) +
  scale_colour_manual(name = "Correlation", values = edge_colors) +
  scale_linetype_manual(name = "Correlation", values = edge_types) +
  legend_theme

# 4. Neutral lines explain absolute Spearman correlation, not its sign.
rho_data <- data.frame(x = 0, xend = 1, y = 1:3,
                       rho = factor(names(rho_widths), levels = names(rho_widths)))
p_rho <- ggplot(rho_data,
                aes(x = x, xend = xend, y = y, yend = y, linewidth = rho)) +
  geom_segment(colour = "#666666", alpha = edge_alpha) +
  scale_linewidth_manual(name = expression(abs(rho)), values = rho_widths) +
  legend_theme

legends <- lapply(list(p_omics, p_degree, p_direction, p_rho), cowplot::get_legend)
stopifnot(all(!vapply(legends, is.null, logical(1))))

# Anchor each legend at the top, including the two-row correlation guide.
top_align <- function(g) {
  box <- gtable::gtable(widths = grid::unit(1, "null"),
                       heights = grid::unit.c(sum(g$heights), grid::unit(1, "null")))
  gtable::gtable_add_grob(box, g, t = 1, l = 1, clip = "off")
}
legend_horizontal <- cowplot::plot_grid(
  plotlist = lapply(legends, top_align), nrow = 1,
  rel_widths = c(1.45, 1.05, 1.25, 0.85)
)
legend_compact <- cowplot::plot_grid(
  plotlist = lapply(legends, top_align), ncol = 2,
  rel_widths = c(1.25, 1)
)

save_legend <- function(plot, stem, width, height) {
  # White preview; transparent PNG and vector files for figure assembly.
  ggsave(file.path(out_dir, paste0(stem, "_preview.png")), plot,
         width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 600, bg = "transparent")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
         device = grDevices::pdf, width = width, height = height,
         units = "in", bg = "transparent", useDingbats = FALSE)
  if (isTRUE(capabilities("cairo"))) {
    ggsave(file.path(out_dir, paste0(stem, ".svg")), plot,
           device = grDevices::svg, width = width, height = height,
           units = "in", bg = "transparent")
  }
}

save_legend(legend_horizontal, "legend_horizontal", width = 8.4, height = 1.8)
save_legend(legend_compact, "legend_compact", width = 4.6, height = 3.6)

# Verify the numerical encodings used to construct the guides.
stopifnot(
  isTRUE(all.equal(unname(degree_sizes / node_size_factor), c(28, 44, 60))),
  isTRUE(all.equal(unname(rho_widths / edge_width_factor), c(1.3, 2.8, 4.3)))
)
expected <- file.path(out_dir, c("legend_horizontal.png", "legend_horizontal.pdf",
                                "legend_compact.png", "legend_compact.pdf"))
stopifnot(all(file.exists(expected)), all(file.info(expected)$size > 0))
if (interactive()) print(legend_horizontal)
message("Legends exported to: ", normalizePath(out_dir, winslash = "/"))
