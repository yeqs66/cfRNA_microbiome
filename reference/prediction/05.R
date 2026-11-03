library(ggplot2)
library(tidyverse)

rm(list = ls())
setwd("__PATH_001__")

load("BMI_kid/predictingchild2yearBMI_usingEARLY.rda")
load("BMI_kid/predictingchild2yearBMI_usingMIDDLE.rda")
load("BMI_kid/predictingchild2yearBMI_usingLATE.rda")


rmse_vec <- sapply(results, function(x) x$rmse_mean)
best_idx <- which.min(rmse_vec)

best_model <- results[[best_idx]]

cat("\nBest model index:", best_idx, "\n")
print(best_model$params)



df_plot <- data.frame(
  model_id = 1:length(results),
  MSE = sapply(results, function(x) x$rmse_mean^2),
  R2  = sapply(results, function(x) x$r2_mean),
  units1 = sapply(results, function(x) x$params$units1),
  units2 = sapply(results, function(x) x$params$units2),
  dropout1 = sapply(results, function(x) x$params$dropout1)
)

# 构建标签（不要太长！）
df_plot$label <- sprintf(
  "MSE=%.3f\nR2=%.2f\nu1=%d u2=%d d=%.2f",
  df_plot$MSE,
  df_plot$R2,
  df_plot$units1,
  df_plot$units2,
  df_plot$dropout1
)

best_idx <- which.min(df_plot$MSE)


cairo_pdf(filename = "BMI_kid/figure_feature_select_middle.pdf",width = 6,height=6,family="Arial")
ggplot(df_plot, aes(x = model_id, y = MSE)) +
  geom_point(size = 3) +
  geom_line(group = 1) +
  geom_text(
    data = df_plot[best_idx, ],
    aes(label = label),
    nudge_y = 0.06,  
    nudge_x = 0.5,
    vjust = -1,
    size = 4
  ) +
  labs(
    title = "2nd Trimester",
    subtitle = "Model Selection (MSE with R² & Parameters)",
    x = "Model Index",
    y = "MSE"
  ) +
  theme_minimal()
  #theme(legend.position = c(0.8, 0.5))
dev.off()





df_pred <- data.frame(
  truth = best_model$trues,
  pred  = best_model$preds
)

cairo_pdf(filename = "BMI_kid/figure_effect_early.pdf",width = 6,height=6,family="Arial")
ggplot(df_pred, aes(x = truth, y = pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "1st trimester: Predicted offspring BMI zscrore vs True",
       subtitle = paste0("MSE=",round(df_plot[best_idx, "MSE"],3),"   R2=",round(df_plot[best_idx, "R2"],3)),
       x = "True",
       y = "Predicted(using 1st trimester cfRNA)") +
  theme_minimal()
dev.off()
