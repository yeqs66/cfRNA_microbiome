# Figure 2: Volcano + GSEA bubble (mock data)
library(ggplot2)
library(ggrepel)
library(dplyr)

# -------------------- Volcano (2A) --------------------
set.seed(123)
n <- 2000
logFC <- rnorm(n, 0, 0.8)
logFC[sample(1:n, 50)] <- runif(50, -3, -1.5)
logFC[sample(1:n, 50)] <- runif(50, 1.5, 3)
pval <- runif(n, 0, 0.2)
pval[sample(1:n, 30)] <- runif(30, 1e-6, 0.05)
logp <- -log10(pval)
SYMBOL <- paste0("Gene_", 1:n)
deg <- data.frame(SYMBOL, logFC, logp, adj.P.Val = pval)
deg <- deg %>%
  mutate(group = case_when(
    logFC > 0.15 & logp > 1 ~ "Up-regulated",
    logFC < -0.15 & logp > 1 ~ "Down-regulated",
    TRUE ~ "Non-significant"
  ))
up <- head(deg$SYMBOL[deg$group == "Up-regulated"], 10)
down <- head(deg$SYMBOL[deg$group == "Down-regulated"], 10)
deg$label <- ifelse(deg$SYMBOL %in% c(up, down), deg$SYMBOL, "")

p1 <- ggplot(deg, aes(x = logFC, y = logp, color = group)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_text_repel(data = subset(deg, group == "Up-regulated" & SYMBOL %in% up),
                  aes(label = SYMBOL), size = 3.5, color = "black",
                  segment.linetype = 3, segment.color = "grey50",
                  nudge_x = 0.5, direction = "y", hjust = 0) +
  geom_text_repel(data = subset(deg, group == "Down-regulated" & SYMBOL %in% down),
                  aes(label = SYMBOL), size = 3.5, color = "black",
                  segment.linetype = 3, segment.color = "grey50",
                  nudge_x = -0.5, direction = "y", hjust = 1) +
  scale_color_manual(values = c("Down-regulated" = "#62B197",
                                "Up-regulated"   = "#E18E6D",
                                "Non-significant"= "grey80")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed", color = "grey30") +
  labs(x = "Log2 fold change", y = "-Log10 FDR", title = "GDM VS Control") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# -------------------- GSEA bubble (2B) --------------------
set.seed(456)
periods <- c("1st", "2nd", "3rd", "All Period")
pathways <- c("Ribosome", "Oxidative phosphorylation", "Parkinson disease",
              "Huntington disease", "Alzheimer disease", "Non-alcoholic fatty liver disease",
              "Thermogenesis", "Chemical carcinogenesis", "Drug metabolism",
              "Metabolism of xenobiotics", "Retinol metabolism", "Fatty acid degradation",
              "Valine, leucine and isoleucine degradation", "Butanoate metabolism",
              "Propanoate metabolism")
gsea <- expand.grid(period = periods, pathway = pathways, stringsAsFactors = FALSE) %>%
  group_by(period) %>%
  mutate(NES = rnorm(n(), 0, 1.2),
         p.adjust = runif(n(), 0.001, 0.2),
         p.adjust = ifelse(runif(n()) < 0.3, runif(n(), 0.001, 0.05), p.adjust),
         log10padj = -log10(p.adjust)) %>%
  ungroup() %>%
  filter(p.adjust < 0.05) %>%
  group_by(period) %>%
  slice_max(order_by = log10padj, n = 5) %>%
  ungroup()
pathway_order <- gsea %>%
  group_by(pathway) %>%
  summarise(mean_log10 = mean(log10padj)) %>%
  arrange(desc(mean_log10)) %>%
  pull(pathway)
gsea$pathway <- factor(gsea$pathway, levels = rev(pathway_order))

p2 <- ggplot(gsea, aes(x = period, y = pathway, size = log10padj, fill = NES)) +
  geom_point(shape = 21, stroke = 0.2) +
  scale_size_continuous(range = c(3, 8), name = "-log10(p.adj)") +
  scale_fill_gradient(low = "#62B197", high = "#E18E6D", name = "NES") +
  ylab("Pathways") + xlab("Trimester") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_line(color = "grey92", linetype = 1, linewidth = 0.2),
        panel.grid.major.x = element_line(color = "grey92", linetype = 1, linewidth = 0.2),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.key = element_rect(fill = "white", colour = NA)) +
  guides(size = guide_legend(title = "-log10(p.adj)"),
         fill = guide_colorbar(title = "NES"))

# -------------------- Output --------------------
print(p1)
print(p2)
# Save if needed: ggsave("Figure2A.pdf", p1); ggsave("Figure2B.pdf", p2)