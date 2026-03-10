rm(list = ls())

# ============================================================
# cfRNA GDM Prediction Model
# ============================================================
# This script builds and evaluates deep learning models to predict
# gestational diabetes mellitus (GDM) using cell-free RNA (cfRNA)
# expression data from first-trimester plasma samples.
#
# Three model variants are evaluated:
#   1. cfRNA features only
#   2. cfRNA features + BMI
#   3. cfRNA features + gut microbiome KO abundances
#
# Workflow:
#   1. Feature selection: differentially expressed cfRNA genes
#   2. 5-fold cross-validated hyperparameter grid search
#   3. Final model training on all data with best hyperparameters
#   4. External validation on two independent cohorts
# ============================================================

library(tidyverse)
library(sjmisc)
library(xgboost)
library(keras)
library(pROC)
library(tensorflow)

# ============================================================
# NOTE: Update the paths below to match your local data location
# before running this script.
# ============================================================
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/eQTL/select_cfRNA_251205.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/eQTL/limma_voom_logTMM_normed_NCBI&lyw_matrix_251217.rda")



# ============================================================
# Section 1: cfRNA-only model (first trimester)
# ============================================================

# Select top differentially expressed genes from first trimester
# with nominal p-value < 0.01
selecttop <-  deg_muiltRNA_withcoef %>%
  filter( period=="1st" & P.Value<0.01 ) 

# Build the gene expression matrix for the training cohort:
# filter to selected genes, transpose, and merge with clinical metadata
filtered_limma_voom_multiRNA <-limma_voom_multiRNA %>%
  rownames_to_column("SYMBOL") %>%
  filter(SYMBOL %in% selecttop$SYMBOL) %>%
  column_to_rownames("SYMBOL") %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  left_join(coldata_1562,by="sample") %>%
  filter(trimester==1)


# Identify common genes between the expression matrix and the selected gene list
selected_genes <- intersect(colnames(filtered_limma_voom_multiRNA), selecttop$SYMBOL)

# Build the model input: selected gene columns + disease outcome
df_model <- filtered_limma_voom_multiRNA %>%
  dplyr:: select(all_of(selected_genes), disease)

# Convert outcome to ordered factor required by the model
df_model$disease <- factor(df_model$disease, levels = c("control", "GDM"))

# Encode outcome as binary integer (0 = control, 1 = GDM)
df_model$label <- ifelse(df_model$disease == "GDM", 1, 0)

# Separate feature matrix and label vector
X_all <- df_model[, setdiff(colnames(df_model), c("disease", "label"))]
Y_all <- df_model$label




################################################################
# Section 2: 5-fold cross-validated hyperparameter grid search
################################################################

# Define the number of cross-validation folds
K <- 5
# Assign each sample to a fold (stratification is random here)
folds <- sample(rep(1:K, length.out = nrow(df_model)))

# Define hyperparameter grid: all combinations will be evaluated
grid <- expand.grid(
  units1 = c(16, 32),      # units in the first hidden layer
  units2 = c(4, 8),        # units in the second hidden layer
  dropout1 = c(0.3, 0.4),  # dropout rate after layer 1
  dropout2 = c(0.2, 0.3),  # dropout rate after layer 2
  l2reg = c(1e-4, 1e-3)    # L2 regularisation strength
)


# Helper: build a two-hidden-layer fully connected network with dropout
# and L2 regularisation for binary classification
build_model <- function(input_dim, units1, units2, dropout1, dropout2, l2reg) {
  keras_model_sequential() %>%
    layer_dense(units = units1, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg),
                input_shape = input_dim) %>%
    layer_dropout(dropout1) %>%
    layer_dense(units = units2, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg)) %>%
    layer_dropout(dropout2) %>%
    layer_dense(units = 1, activation = "sigmoid") %>%
    compile(
      optimizer = optimizer_adam(),
      loss = "binary_crossentropy",
      metrics = c("AUC")
    )
}

################################################################
# Section 3: Train and evaluate each hyperparameter combination
#            (skip if loading pre-computed results below)
################################################################
results <- list()
for (i in 1:nrow(grid)) {
  params <- grid[i, ]
  cat(paste0("\n========== Testing params ", i, "/", nrow(grid), " ==========\n"))
  print(params)
  auc_list <- c()
  # Lists to store predictions and true labels for the current parameter set across all folds
  all_preds <- list()
  all_trues <- list()
  
  for (k in 1:K) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:nrow(df_model), test_idx)
    X_train <- X_all[train_idx, ]
    X_test <- X_all[test_idx, ]
    y_train <- Y_all[train_idx]
    y_test <- Y_all[test_idx]
    
    # Standardise features using training-fold mean and SD only
    # (avoids data leakage from test fold)
    m <- apply(X_train, 2, mean)
    s <- apply(X_train, 2, sd)
    X_train <- scale(X_train, center = m, scale = s)
    X_test <- scale(X_test, center = m, scale = s)
    X_train <- as.matrix(X_train)
    X_test <- as.matrix(X_test)
    
    # Build and train the model for this fold
    model <- build_model(
      input_dim = ncol(X_train),
      units1 = params$units1,
      units2 = params$units2,
      dropout1 = params$dropout1,
      dropout2 = params$dropout2,
      l2reg = params$l2reg
    )
    model %>% fit(
      X_train, y_train,
      epochs = 50,
      batch_size = 16,
      verbose = 0
    )
    pred <- as.numeric(model %>% predict(X_test))
    auc_fold <- roc(y_test, pred)$auc
    auc_list <- c(auc_list, auc_fold)
    
    # Accumulate out-of-fold predictions for final ROC curve
    all_preds[[k]] <- pred
    all_trues[[k]] <- y_test
  }
  results[[i]] <- list(
    params = params,
    auc_mean = mean(auc_list),
    auc_sd = sd(auc_list),
    all_preds = unlist(all_preds),
    all_trues = unlist(all_trues)
  )
  cat(paste0("AUC MEAN = ", round(mean(auc_list), 4), " SD = ", round(sd(auc_list), 4), "\n"))
}

#########################################################
# Section 4: Load pre-computed results and select the best model
# (start here if you skipped the training loop above)
#########################################################

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/veryimportantparametersfortraincfRNA")
# Find the hyperparameter combination with the highest mean cross-validated AUC
auc_values <- sapply(results, function(x) x$auc_mean)
best_idx <- which.max(auc_values)
cat("\n=========== BEST MODEL ===========\n")
print(results[[best_idx]]$params)
best_model_info <- results[[best_idx]]
cat("\n=========== BEST MODEL ===========\n")
print(best_model_info$params)

#save(results,file="C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/veryimportantparametersfortraincfRNA")

# Plot the cross-validated ROC curve for the best hyperparameter set
# using the accumulated out-of-fold predictions
final_trues <- best_model_info$all_trues
final_preds <- best_model_info$all_preds

# Generate the ROC curve object
best_roc <- roc(final_trues, final_preds)
# 95% CI for AUC using the DeLong method
auc_ci <- ci.auc(best_roc)


# Plot the ROC curve
plot.roc(best_roc, main = "ROC Curve for Best Model (Cross-Validated)",
     lwd = 2,
     xlim = c(1, 0), 
     ylim = c(0, 1),
     asp = NA
     )

legend(
  "bottomright",
  legend = paste0(
    "Only cfRNA first trimester \n AUC = ", round(best_model_info[["auc_mean"]], 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "–",
    round(auc_ci[3], 3)
  ),
  bty = "n"
)

best_model_info <- results[[best_idx]]
params <- best_model_info$params


final_model <- build_model(
  input_dim = ncol(X_all),
  units1   = results[[best_idx]]$params$units1,
  units2   = results[[best_idx]]$params$units2,
  dropout1 = results[[best_idx]]$params$dropout1,
  dropout2 = results[[best_idx]]$params$dropout2,
  l2reg    = results[[best_idx]]$params$l2reg
)


save_model_tf(
  final_model,
  filepath = "final_cfRNA_model"
)



# Standardise the full training dataset for final model training
m_all <- apply(X_all, 2, mean)
s_all <- apply(X_all, 2, sd)

X_all_scaled <- scale(X_all, center = m_all, scale = s_all)



#############################################################
# Section 5: External validation (overall)
# Retrain the final model with fixed random seed for reproducibility,
# then evaluate on two independent external cohorts.
#############################################################

  tensorflow::tf$keras$backend$clear_session()  
  
  set.seed(1997)
  Sys.setenv("PYTHONHASHSEED" = "0")
  Sys.setenv("TF_DETERMINISTIC_OPS" = "1")
  Sys.setenv("TF_CUDNN_DETERMINISTIC" = "1")
  Sys.setenv("CUDA_VISIBLE_DEVICES" = "")
  
  tf$random$set_seed(1997)
  tf$keras$utils$set_random_seed(as.integer(1997))
  # Train the final model on all training data using the best hyperparameters
  final_model <- build_model(
  input_dim = ncol(X_all),
  units1   = results[[11]]$params$units1,
  units2   = results[[11]]$params$units2,
  dropout1 = results[[11]]$params$dropout1,
  dropout2 = results[[11]]$params$dropout2,
  l2reg    = results[[11]]$params$l2reg
)

final_model %>% fit(
  as.matrix(X_all_scaled),
  Y_all,
  epochs = 50,
  batch_size = 16,
  verbose = 0
) 

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/eQTL/limma_voom_logTMM_normed_NCBI&lyw_matrix_251217.rda")

coldata_verify <- coldata_lyw %>%
  dplyr::select(sample,condition,trimester) 

lyw <- as.data.frame(limma_voom_multiRNA_lyw) %>%
  rownames_to_column("SYMBOL") %>%
  filter(SYMBOL %in% selecttop$SYMBOL) %>%
  column_to_rownames("SYMBOL") %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  left_join(coldata_verify,by="sample") %>%
  filter(trimester=="1st") %>%
  dplyr::select(-trimester)




train_features <- colnames(X_all)

X_ext_raw <- limma_voom_multiRNA   # matrix, samples x genes
## outcome (convert 1/2 → 0/1)
y_ext_raw <- lyw$condition
Y_ext_lyw <- ifelse(y_ext_raw == "GDM", 1, 0)

## external feature matrix (drop Group)
X_ext_raw <- lyw[, setdiff(colnames(lyw), "condition"), drop = FALSE]

## training feature names
train_features <- colnames(X_all)

## keep only features present in both training and external data
X_ext <- X_ext_raw[, intersect(colnames(X_ext_raw), train_features), drop = FALSE]

## fill missing training features with zero (gene not detected in external cohort)
missing <- setdiff(train_features, colnames(X_ext))
for(m in missing){
  X_ext[[m]] <- 0
}

## reorder columns to match training feature order
X_ext <- X_ext[, train_features]

## scale using training mean and SD (must use training statistics to avoid leakage)
X_ext_scaled_lyw <- scale(X_ext, center = m_all, scale = s_all)



pred_ext_lyw <- as.numeric(final_model %>% predict(as.matrix(X_ext_scaled_lyw)))

roc_ext <- roc(
  response = Y_ext_lyw,
  predictor = pred_ext_lyw,
  quiet = TRUE,
  direction = "<"
)

auc_ext <- auc(roc_ext)
auc_ci <- ci.auc(roc_ext)

plot.roc(
  roc_ext,
  lwd = 2,
  xlim = c(1, 0), 
  ylim = c(0, 1),
  asp = NA,
  main = paste0(
    "Yuwei Liu et al. validation \nAUC = ",
    round(auc_ext, 3)
  )
)


legend(
  "bottomright",
  legend = paste0(
    "Yuwei Liu  et al. validation,\n AUC = ", round(auc_ext, 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "–",
    round(auc_ci[3], 3)
  ),
  bty = "n"
)



##############################################
# External validation: Giorgia Del Vecchio et al. (NCBI) cohort
##############################################

coldata_verify <- coldata_NCBI %>%
  dplyr::select(sample,condition,trimester) 

ncbi <- as.data.frame(limma_voom_NCBI) %>%
  rownames_to_column("SYMBOL") %>%
  filter(SYMBOL %in% selecttop$SYMBOL) %>%
  column_to_rownames("SYMBOL") %>%
  rotate_df() %>%
  rownames_to_column("sample") %>%
  left_join(coldata_verify,by="sample") %>%
  filter(trimester=="1st") %>%
  dplyr::select(-trimester)



train_features <- colnames(X_all)

## outcome
y_ext_raw <- ncbi$condition
Y_ext_ncbi <- ifelse(y_ext_raw == "GDM", 1, 0)

## external feature matrix (drop outcome column)
X_ext_raw <- ncbi[, setdiff(colnames(ncbi), "condition"), drop = FALSE]

## training feature names
train_features <- colnames(X_all)

## keep only features present in both training and external data
X_ext <- X_ext_raw[, intersect(colnames(X_ext_raw), train_features), drop = FALSE]

## fill missing training features with zero
missing <- setdiff(train_features, colnames(X_ext))
for(m in missing){
  X_ext[[m]] <- 0
}

## reorder columns to match training feature order
X_ext <- X_ext[, train_features]

## scale using training mean and SD
X_ext_scaled_ncbi <- scale(X_ext, center = m_all, scale = s_all)


pred_ext_ncbi <- as.numeric(final_model %>% predict(as.matrix(X_ext_scaled_ncbi)))

roc_ext <- roc(
  response = Y_ext_ncbi,
  predictor = pred_ext_ncbi,
  quiet = TRUE,
  direction = "<"
)

auc_ext <- auc(roc_ext)
auc_ci <- ci.auc(roc_ext)

plot.roc(
  roc_ext,
  lwd = 2,
  xlim = c(1, 0), 
  ylim = c(0, 1),
  asp = NA,
  main = paste0(
    "Giorgia Del Vecchio  et al. validation \nAUC = ",
    round(auc_ext, 3)
  )
)


legend(
  "bottomright",
  legend = paste0(
    "Giorgia Del Vecchio  et al.,\n AUC = ", round(auc_ext, 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "–",
    round(auc_ci[3], 3)
  ),
  bty = "n"
)



#####################################################
# Combined external validation (both cohorts pooled)
#####################################################

x_ext_combine <-rbind(as.matrix(X_ext_scaled_lyw),as.matrix(X_ext_scaled_ncbi))
Y_ext_combine <- c(Y_ext_lyw,Y_ext_ncbi)


pred_ext_combine <- as.numeric(final_model %>% predict(as.matrix(x_ext_combine)))

roc_ext <- roc(
  response = Y_ext_combine,
  predictor = pred_ext_combine,
  quiet = TRUE,
  direction = "<"
)

auc_ext <- auc(roc_ext)
auc_ci <- ci.auc(roc_ext)

plot.roc(
  roc_ext,
  lwd = 2,
  xlim = c(1, 0), 
  ylim = c(0, 1),
  asp = NA,
  main = paste0(
    "Combine two external \nAUC = ",
    round(auc_ext, 3)
  )
)


legend(
  "bottomright",
  legend = paste0(
    "Combine external data,\n AUC = ", round(auc_ext, 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "–",
    round(auc_ci[3], 3)
  ),
  bty = "n"
)






####################################################
# Section 6: cfRNA + BMI model (internal cross-validation)
####################################################

# Rebuild the model input including BMI as an additional feature
selected_genes <- intersect(colnames(filtered_limma_voom_multiRNA), selecttop$SYMBOL)

df_model <- filtered_limma_voom_multiRNA %>%
  dplyr::select(all_of(selected_genes),BMI_mo_Anthropometric, disease)

df_model$disease <- factor(df_model$disease, levels = c("control", "GDM"))

df_model$label <- ifelse(df_model$disease == "GDM", 1, 0)

# Feature matrix now includes both gene expression and BMI
X_all_BMI <- df_model[, setdiff(colnames(df_model), c("disease", "label"))]
Y_all_BMI <- df_model$label

K <- 5
folds <- sample(rep(1:K, length.out = nrow(df_model)))

# Re-use the same hyperparameter grid
grid <- expand.grid(
  units1 = c(16, 32),
  units2 = c(4, 8),
  dropout1 = c(0.3, 0.4),
  dropout2 = c(0.2, 0.3),
  l2reg = c(1e-4, 1e-3)
)

results <- list()

# build_model is unchanged — the same architecture applies regardless of input dimension
build_model <- function(input_dim, units1, units2, dropout1, dropout2, l2reg) {
  keras_model_sequential() %>%
    layer_dense(units = units1, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg),
                input_shape = input_dim) %>%
    layer_dropout(dropout1) %>%
    layer_dense(units = units2, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg)) %>%
    layer_dropout(dropout2) %>%
    layer_dense(units = 1, activation = "sigmoid") %>%
    compile(
      optimizer = optimizer_adam(),
      loss = "binary_crossentropy",
      metrics = c("AUC")
    )
}

####################################################
# Skip this loop if loading pre-computed bmi_results below
####################################################
bmi_results <- list()


for (i in 1:nrow(grid)) {
  params <- grid[i, ]
  cat(paste0("\n========== Testing params ", i, "/", nrow(grid), " ==========\n"))
  print(params)
  auc_list <- c()
  all_preds <- list()
  all_trues <- list()
  
  for (k in 1:K) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:nrow(df_model), test_idx)
    X_train <- X_all_BMI[train_idx, ]
    X_test <- X_all_BMI[test_idx, ]
    y_train <- Y_all_BMI[train_idx]
    y_test <- Y_all_BMI[test_idx]
    
    # Standardise within fold to prevent data leakage
    m <- apply(X_train, 2, mean)
    s <- apply(X_train, 2, sd)
    X_train <- scale(X_train, center = m, scale = s)
    X_test <- scale(X_test, center = m, scale = s)
    X_train <- as.matrix(X_train)
    X_test <- as.matrix(X_test)
    
    model <- build_model(
      input_dim = ncol(X_train),
      units1 = params$units1,
      units2 = params$units2,
      dropout1 = params$dropout1,
      dropout2 = params$dropout2,
      l2reg = params$l2reg
    )
    model %>% fit(
      X_train, y_train,
      epochs = 50,
      batch_size = 16,
      verbose = 0
    )
    pred <- as.numeric(model %>% predict(X_test))
    auc_fold <- roc(y_test, pred)$auc
    auc_list <- c(auc_list, auc_fold)
    
    all_preds[[k]] <- pred
    all_trues[[k]] <- y_test
  }
  
  bmi_results[[i]] <- list(
    params = params,
    auc_mean = mean(auc_list),
    auc_sd = sd(auc_list),
    all_preds = unlist(all_preds),
    all_trues = unlist(all_trues)
  )
  cat(paste0("AUC MEAN = ", round(mean(auc_list), 4), " SD = ", round(sd(auc_list), 4), "\n"))
}

save(bmi_results, file="C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/veryimportantparametersfortraincfRNA_BMI.rda")




##############################################################
# Section 7: Load pre-computed BMI model results and plot ROC
##############################################################
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/veryimportantparametersfortraincfRNA_BMI.rda")
# Find the best hyperparameter combination for the cfRNA + BMI model
auc_values <- sapply(bmi_results, function(x) x$auc_mean)
best_idx <- which.max(auc_values)
cat("\n=========== BEST MODEL ===========\n")
print(bmi_results[[best_idx]]$params)
best_model_info <- bmi_results[[best_idx]]
cat("\n=========== BEST MODEL ===========\n")
print(best_model_info$params)


# Retrieve out-of-fold predictions accumulated during cross-validation
final_trues <- best_model_info$all_trues
final_preds <- best_model_info$all_preds

# Generate the ROC curve object
best_roc <- roc(final_trues, final_preds)

roc_ext <- roc(
  response = final_trues,
  predictor = final_preds,
  quiet = TRUE,
  direction = "<"
)

auc_ext <- auc(roc_ext)
# ci.auc() requires a roc object (not a numeric AUC value)
auc_ci <- ci.auc(roc_ext)

plot.roc(
  roc_ext,
  lwd = 2,
  xlim = c(1, 0), 
  ylim = c(0, 1),
  asp = NA,
  main = paste0(
    "ROC Curve for Best Model (Cross-Validated)"
  )
)


legend(
  "bottomright",
  legend = paste0(
    "cfRNA+BMI first trimester","\n AUC = ",round(auc_ext, 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "–",
    round(auc_ci[3], 3)
  ),
  bty = "n"
)


####################################################
# Section 8: cfRNA + microbiome KO abundance model
####################################################
rm(list = ls())

load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/eQTL/select_cfRNA_251205.rda")
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/test_prediction_260202.rda")


selecttop <-  deg_muiltRNA_withcoef %>%
  filter( period=="1st" & P.Value<0.01 ) 

# Select significantly associated microbiome KEGG orthologs (KOs)
selected <- deg_ko_merge %>%
  dplyr::filter( p.adj_overall <0.001)

# Subset the KO abundance matrix to significant features
rpm_KO_1476_sig <- rpm_KO_1476_sig[selected$feature,]

selected_genes <- intersect(rownames(limma_voom_multiRNA_1476), selecttop$SYMBOL)

DISEASEinfo <- coldata_1476 %>%
  dplyr::select(ID,disease)
# Combine cfRNA and microbiome KO features, subset to early-trimester samples
selectlimma_voom_multiRNA_1476 <- limma_voom_multiRNA_1476[selected_genes,]
df_model <- rbind(selectlimma_voom_multiRNA_1476,rpm_KO_1476_sig) %>%
  dplyr::select(contains("early")) %>%
  rotate_df() %>%
  rownames_to_column("ID") %>%
  mutate(ID = str_remove(ID,"early")) %>%
  left_join(DISEASEinfo,by="ID")

df_model$disease <- factor(df_model$disease, levels = c("control", "GDM"))

df_model$label <- ifelse(df_model$disease == "GDM", 1, 0)

# Feature matrix combines cfRNA expression and microbiome KO abundances
X_all_ko_fcRNA <- df_model[, setdiff(colnames(df_model), c("ID","disease", "label"))]
Y_all_ko_fcRNA  <- df_model$label

K <- 5
folds <- sample(rep(1:K, length.out = nrow(df_model)))

# Re-use the same hyperparameter grid
grid <- expand.grid(
  units1 = c(16, 32),
  units2 = c(4, 8),
  dropout1 = c(0.3, 0.4),
  dropout2 = c(0.2, 0.3),
  l2reg = c(1e-4, 1e-3)
)


build_model <- function(input_dim, units1, units2, dropout1, dropout2, l2reg) {
  keras_model_sequential() %>%
    layer_dense(units = units1, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg),
                input_shape = input_dim) %>%
    layer_dropout(dropout1) %>%
    layer_dense(units = units2, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg)) %>%
    layer_dropout(dropout2) %>%
    layer_dense(units = 1, activation = "sigmoid") %>%
    compile(
      optimizer = optimizer_adam(),
      loss = "binary_crossentropy",
      metrics = c("AUC")
    )
}

KOcfRNA_results <- list()
for (i in 1:nrow(grid)) {
  params <- grid[i, ]
  cat(paste0("\n========== Testing params ", i, "/", nrow(grid), " ==========\n"))
  print(params)
  auc_list <- c()
  all_preds <- list()
  all_trues <- list()
  
  for (k in 1:K) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:nrow(df_model), test_idx)
    X_train <- X_all_ko_fcRNA[train_idx, ]
    X_test <- X_all_ko_fcRNA[test_idx, ]
    y_train <- Y_all_ko_fcRNA[train_idx]
    y_test <- Y_all_ko_fcRNA[test_idx]
    
    # Standardise within fold to prevent data leakage
    m <- apply(X_train, 2, mean)
    s <- apply(X_train, 2, sd)
    X_train <- scale(X_train, center = m, scale = s)
    X_test <- scale(X_test, center = m, scale = s)
    X_train <- as.matrix(X_train)
    X_test <- as.matrix(X_test)
    
    model <- build_model(
      input_dim = ncol(X_train),
      units1 = params$units1,
      units2 = params$units2,
      dropout1 = params$dropout1,
      dropout2 = params$dropout2,
      l2reg = params$l2reg
    )
    model %>% fit(
      X_train, y_train,
      epochs = 50,
      batch_size = 16,
      verbose = 0
    )
    pred <- as.numeric(model %>% predict(X_test))
    auc_fold <- roc(y_test, pred)$auc
    auc_list <- c(auc_list, auc_fold)
    
    all_preds[[k]] <- pred
    all_trues[[k]] <- y_test
  }
  
  KOcfRNA_results[[i]] <- list(
    params = params,
    auc_mean = mean(auc_list),
    auc_sd = sd(auc_list),
    all_preds = unlist(all_preds),
    all_trues = unlist(all_trues)
  )
  cat(paste0("AUC MEAN = ", round(mean(auc_list), 4), " SD = ", round(sd(auc_list), 4), "\n"))
}

save(KOcfRNA_results, file="C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/microbiome_cfRNA_BMI.rda")


######################################################################
# Section 9: Load pre-computed KO+cfRNA results and plot ROC
######################################################################
load("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/microbiome_cfRNA_BMI.rda")
# Find the best hyperparameter combination for the cfRNA + microbiome model
auc_values <- sapply(KOcfRNA_results, function(x) x$auc_mean)
best_idx <- which.max(auc_values)
cat("\n=========== BEST MODEL ===========\n")
print(KOcfRNA_results[[best_idx]]$params)
best_model_info <- KOcfRNA_results[[best_idx]]
cat("\n=========== BEST MODEL ===========\n")
print(best_model_info$params)


# Retrieve out-of-fold predictions accumulated during cross-validation
final_trues <- best_model_info$all_trues
final_preds <- best_model_info$all_preds

# Generate the ROC curve object
best_roc <- roc(final_trues, final_preds)

roc_ext <- roc(
  response = final_trues,
  predictor = final_preds,
  quiet = TRUE,
  direction = "<"
)

auc_ext <- auc(roc_ext)
# ci.auc() requires a roc object (not a numeric AUC value)
auc_ci <- ci.auc(roc_ext)

plot.roc(
  roc_ext,
  lwd = 2,
  xlim = c(1, 0), 
  ylim = c(0, 1),
  asp = NA,
  main = paste0(
    "ROC Curve for Best Model (Cross-Validated)"
  )
)


legend(
  "bottomright",
  legend = paste0(
    "cfRNA+microbiome first trimester","\n AUC = ",round(auc_ext, 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "–",
    round(auc_ci[3], 3)
  ),
  bty = "n"
)
