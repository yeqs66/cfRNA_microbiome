# Figure 6: Predictive models (ROC curves + BMI prediction)
library(ggplot2)
library(dplyr)
library(pROC)
library(patchwork)

# -------------------- Helper: generate ROC data --------------------
gen_roc_data <- function(n = 200, auc_true = 0.75, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Simulate scores and labels with given AUC (via probit)
  z <- rnorm(n)
  score <- z + rnorm(n, 0, 0.5)  # noise
  # threshold to achieve desired AUC (approx)
  # use logistic to convert to probability, then sample labels
  p <- plogis(score)
  label <- rbinom(n, 1, p)
  # Adjust to desired AUC by reweighting (simple method: use roc function to compute, we'll just return)
  # Return data frame
  data.frame(label = label, score = score)
}

# Or simpler: directly generate ROC curves with known AUC via binormal model
sim_roc <- function(n = 200, auc = 0.75, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # binormal: control mean=0 sd=1, case mean=d sd=1, d = sqrt(2)*qnorm(auc)
  d <- sqrt(2) * qnorm(auc)
  labels <- c(rep(0, n/2), rep(1, n/2))
  scores <- c(rnorm(n/2, 0, 1), rnorm(n/2, d, 1))
  list(labels = labels, scores = scores)
}

# -------------------- 6A: Internal cfRNA ROC --------------------
rocA <- sim_roc(n = 300, auc = 0.82, seed = 123)
roc_objA <- roc(rocA$labels, rocA$scores)
ciA <- ci.auc(roc_objA)

pA <- ggroc(roc_objA, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "A: cfRNA (Internal)",
       subtitle = paste0("AUC = ", round(auc(roc_objA), 3),
                         " (95% CI: ", round(ciA[1], 3), "-", round(ciA[3], 3), ")")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# -------------------- 6B: External validation (combined) --------------------
rocB <- sim_roc(n = 200, auc = 0.78, seed = 456)
roc_objB <- roc(rocB$labels, rocB$scores)
ciB <- ci.auc(roc_objB)
pB <- ggroc(roc_objB, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "B: cfRNA (External combined)",
       subtitle = paste0("AUC = ", round(auc(roc_objB), 3),
                         " (95% CI: ", round(ciB[1], 3), "-", round(ciB[3], 3), ")")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# -------------------- 6C: cfRNA + BMI --------------------
rocC <- sim_roc(n = 300, auc = 0.88, seed = 789)
roc_objC <- roc(rocC$labels, rocC$scores)
ciC <- ci.auc(roc_objC)
pC <- ggroc(roc_objC, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "C: cfRNA + BMI",
       subtitle = paste0("AUC = ", round(auc(roc_objC), 3),
                         " (95% CI: ", round(ciC[1], 3), "-", round(ciC[3], 3), ")")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# -------------------- 6D: cfRNA + Microbiome KO --------------------
rocD <- sim_roc(n = 300, auc = 0.91, seed = 101)
roc_objD <- roc(rocD$labels, rocD$scores)
ciD <- ci.auc(roc_objD)
pD <- ggroc(roc_objD, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "D: cfRNA + Microbiome KO",
       subtitle = paste0("AUC = ", round(auc(roc_objD), 3),
                         " (95% CI: ", round(ciD[1], 3), "-", round(ciD[3], 3), ")")) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# Combine 6A-6D into a 2x2 grid
p_roc_grid <- (pA + pB) / (pC + pD) +
  plot_annotation(title = "Figure 6A-D: ROC curves for GDM prediction") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -------------------- 6E-F: BMI z-score prediction (scatter + model selection) --------------------
# 6E: Model selection (MSE vs model index)
set.seed(202)
n_models <- 20
model_df <- data.frame(
  model_id = 1:n_models,
  MSE = runif(n_models, 0.2, 0.8) + (1:n_models)/50,  # increasing trend
  R2 = runif(n_models, 0.2, 0.6)
)
# Add some structure so best is around model 8
model_df$MSE[8] <- 0.25
model_df$R2[8] <- 0.65
best_idx <- which.min(model_df$MSE)

pE <- ggplot(model_df, aes(x = model_id, y = MSE)) +
  geom_line(color = "grey70") +
  geom_point(size = 3, color = "steelblue") +
  geom_point(data = model_df[best_idx, ], aes(x = model_id, y = MSE),
             size = 4, color = "red") +
  geom_text(data = model_df[best_idx, ],
            aes(label = paste0("MSE=", round(MSE,3), "\nR2=", round(R2,3), "\nu1=32,u2=8,d=0.2")),
            nudge_x = 2, nudge_y = 0.05, size = 3) +
  labs(title = "E: Model selection (2nd trimester)",
       x = "Model Index", y = "MSE") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# 6F: Predicted vs true BMI z-score
set.seed(303)
n_test <- 100
true_bmi <- rnorm(n_test, 0, 1)
pred_bmi <- true_bmi + rnorm(n_test, 0, 0.4)
# Ensure some correlation
cor_est <- cor(true_bmi, pred_bmi)
r2 <- cor_est^2

pF <- ggplot(data.frame(truth = true_bmi, pred = pred_bmi),
             aes(x = truth, y = pred)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "F: Predicted vs True BMI z-score",
       subtitle = paste0("MSE = ", round(mean((true_bmi - pred_bmi)^2), 3),
                         ", R2 = ", round(r2, 3)),
       x = "True BMI z-score", y = "Predicted BMI z-score") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# Combine E and F side by side
p_ef <- pE + pF + plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Figure 6E-F: BMI prediction") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -------------------- Output --------------------
print("Figure 6A-D (ROC curves):")
print(p_roc_grid)

print("Figure 6E-F (BMI prediction):")
print(p_ef)