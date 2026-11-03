library(patchwork)

load("res_pred_miRNA.rda")

pathway_sig001<-unique(res_target_miRNA_kegg$Description[res_target_miRNA_kegg$p.adjust<0.01 & res_target_miRNA_kegg$factor %in% c(1,2)])

data_plot<-res_target_miRNA_kegg %>%
  dplyr::filter(factor %in% c(1,2),Description %in% pathway_sig001) %>%
  arrange(p.adjust) %>%
  separate(GeneRatio, into = c("a", "b"), sep = "/", convert = TRUE) %>% # 分割为 "a" 和 "b"，并转为数值型
  mutate(GeneRatio = a / b,
         log10padj= -log10(p.adjust)) %>%
  dplyr::select(-a,-b) %>%
  pivot_wider(names_from = factor, values_from = c(GeneRatio, log10padj), names_prefix = "factor",values_fill = 0)

p1 <- ggplot(data_plot, aes(x = GeneRatio_factor1, y = Description)) +
  geom_col(fill = "#a52a2a", width = 0.1) +
  geom_point(aes(size = log10padj_factor1), color = "#a52a2a") +
  scale_size_continuous(range = c(0, 8)) +
  scale_x_continuous(limits = c(0, 0.2),  breaks = seq(0, 5, 0.2)) + #x轴刻度自定义
  labs(title = 'factor1') + #添加标题
  theme_bw() + 
  labs(size="-log10(p.adjust)") +
  theme(
    plot.title = element_text(face = 'bold',  color = '#a42929', size = 8), # #修改标题字体、字号、颜色
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.2,size = 0.3),
    axis.text.y = element_text(size = 8,colour = "black",hjust = 0.5),
    axis.text.x = element_text(size = 8,face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p2 <- ggplot(data_plot, aes(x = GeneRatio_factor2, y = Description)) +
  geom_col(fill = "#dc968d", width = 0.1) +
  geom_point(aes(size = log10padj_factor2), color = "#dc968d") +
  scale_size_continuous(range = c(0, 8)) +
  scale_y_discrete(position ="right") +
  scale_x_reverse(limits = c(0.2, 0), breaks = seq(0, 5, 0.2)) +
  labs(title = 'factor2') + #添加标题
  theme_bw() + 
  theme(
    plot.title = element_text(face ='bold', color ='#dc968d', size = 8, hjust = 1),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.2,size = 0.3),
    axis.text.x = element_text(size = 8,face ="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = "None"
  )


p <- (p2 | p1) + plot_layout(guides ='collect') &
  theme( legend.justification = c("right","bottom"),
         legend.text = element_text(face ="bold",size=15)# 移除图例标题
  )

cairo_pdf("enrichment_kegg_miRNA.pdf", width = 9, height = 5, family = "Arial")
p
dev.off()




