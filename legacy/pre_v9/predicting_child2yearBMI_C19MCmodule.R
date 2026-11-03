rm(list = ls())

library(tidyverse)
library(sjmisc)
library(xgboost)
library(keras)
library(pROC)
library(tensorflow)
library(conflicted)

conflicted::conflicts_prefer(pROC::roc)
conflicts_prefer(dplyr::filter)

load("~/349/349/test_bmi_z_kid_prediction_260320.rda")


forprediction <- limma_voom_module13 %>%
  rotate_df() %>%
  mutate(id =row.names(.) ) %>%
  mutate(sjid_mo = sub("[^0-9].*", "", id) ) %>%
  filter(str_detect(id,"late")) %>%
  dplyr::select(-id)


coldata_1476_bmizkid <- coldata_1476_bmizkid %>%
  dplyr::select(sjid_mo ,bmi_z_kid) %>%
  inner_join(forprediction,by ="sjid_mo") %>%
  dplyr::select(-sjid_mo)


df_model <- coldata_1476_bmizkid


df_model$label <-coldata_1476_bmizkid$bmi_z_kid


X_all <- df_model[, setdiff(colnames(df_model), c("bmi_z_kid", "label"))]
Y_all <- df_model$label



X_all <- as.matrix(X_all)

############################################################
# CV setup
############################################################
K <- 5
folds <- sample(rep(1:K, length.out = nrow(X_all)))

############################################################
# Hyperparameter grid
############################################################
grid <- expand.grid(
  units1   = c(32, 64),
  units2   = c(8, 16),
  dropout1 = c(0.1, 0.2),
  dropout2 = c(0.1, 0.2),
  l2reg    = c(1e-4, 1e-3),
  lr       = c(1e-3)
)

############################################################
# Model builder (REGRESSION)
############################################################
build_model <- function(input_dim, units1, units2, dropout1, dropout2, l2reg, lr) {
  
  keras_model_sequential() %>%
    layer_dense(units = units1,
                activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg),
                input_shape = input_dim) %>%
    
    layer_batch_normalization() %>%
    layer_dropout(dropout1) %>%
    
    layer_dense(units = units2,
                activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg)) %>%
    
    layer_batch_normalization() %>%
    layer_dropout(dropout2) %>%
    
    layer_dense(units = 1, activation = "linear") %>%   # regression output
    
    compile(
      optimizer = optimizer_adam(learning_rate = lr),
      loss = "mse",
      metrics = list("mae")
    )
}

############################################################
# Training loop
############################################################
results <- list()

for (i in 1:nrow(grid)) {
  
  params <- grid[i, ]
  cat("\n========== Grid", i, "/", nrow(grid), "==========\n")
  print(params)
  
  rmse_list <- c()
  mae_list  <- c()
  r2_list   <- c()
  
  all_preds <- c()
  all_trues <- c()
  
  for (k in 1:K) {
    
    test_idx  <- which(folds == k)
    train_idx <- setdiff(1:nrow(X_all), test_idx)
    
    X_train <- X_all[train_idx, ]
    X_test  <- X_all[test_idx, ]
    y_train <- Y_all[train_idx]
    y_test  <- Y_all[test_idx]
    
    ########################################################
    # Scaling (within fold, leakage-safe)
    ########################################################
    m <- colMeans(X_train)
    s <- apply(X_train, 2, sd)
    s[s == 0] <- 1
    
    X_train <- scale(X_train, center = m, scale = s)
    X_test  <- scale(X_test,  center = m, scale = s)
    
    ########################################################
    # Model
    ########################################################
    model <- build_model(
      input_dim = ncol(X_train),
      units1   = params$units1,
      units2   = params$units2,
      dropout1 = params$dropout1,
      dropout2 = params$dropout2,
      l2reg    = params$l2reg,
      lr       = params$lr
    )
    
    ########################################################
    # Train (with early stopping)
    ########################################################
    history <- model %>% fit(
      X_train, y_train,
      epochs = 200,
      batch_size = 16,
      validation_split = 0.2,
      verbose = 0,
      callbacks = list(
        callback_early_stopping(
          patience = 15,
          restore_best_weights = TRUE
        )
      )
    )
    
    ########################################################
    # Predict
    ########################################################
    pred <- as.numeric(model %>% predict(X_test))
    
    ########################################################
    # Metrics
    ########################################################
    rmse <- sqrt(mean((y_test - pred)^2))
    mae  <- mean(abs(y_test - pred))
    r2   <- cor(y_test, pred)^2
    
    rmse_list <- c(rmse_list, rmse)
    mae_list  <- c(mae_list, mae)
    r2_list   <- c(r2_list, r2)
    
    all_preds <- c(all_preds, pred)
    all_trues <- c(all_trues, y_test)
  }
  
  ##########################################################
  # Store results
  ##########################################################
  results[[i]] <- list(
    params = params,
    rmse_mean = mean(rmse_list),
    rmse_sd   = sd(rmse_list),
    mae_mean  = mean(mae_list),
    r2_mean   = mean(r2_list),
    preds     = all_preds,
    trues     = all_trues
  )
  
  cat(sprintf("RMSE = %.4f ± %.4f | MAE = %.4f | R2 = %.4f\n",
              mean(rmse_list), sd(rmse_list),
              mean(mae_list), mean(r2_list)))
}

############################################################
# Select best model (by RMSE)
############################################################
rmse_vec <- sapply(results, function(x) x$rmse_mean)
best_idx <- which.min(rmse_vec)

best_model <- results[[best_idx]]

cat("\nBest model index:", best_idx, "\n")
print(best_model$params)



save(results,file="predictingchild2yearBMI_usingLATE.rda")



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

ggplot(df_plot, aes(x = model_id, y = MSE)) +
  geom_point(size = 3) +
  geom_line(group = 1) +
  
  # 只标最优点（关键！）
  geom_text(
    data = df_plot[best_idx, ],
    aes(label = label),
    vjust = -1,
    size = 4
  ) +
  
  labs(
    title = "Model Selection (MSE with R² & Parameters)",
    x = "Model Index",
    y = "MSE"
  ) +
  theme_minimal()






df_pred <- data.frame(
  truth = best_model$trues,
  pred  = best_model$preds
)

ggplot(df_pred, aes(x = truth, y = pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Predicted offspring BMI zscrore vs True",
       subtitle = paste0("MSE=",round(df_plot[best_idx, "MSE"],3),"   R2=",round(df_plot[best_idx, "R2"],3)),
       x = "True",
       y = "Predicted(using middle pregnancy cfRNA)") +
  theme_minimal()