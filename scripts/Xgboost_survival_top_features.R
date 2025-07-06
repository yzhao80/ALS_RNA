# Survival prediction using features selected by XGBoost
# Author: Minghua Li
# Date: 2025-1-20

# Load required libraries
library(caret)
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(timeROC)
library(xgboost)
library(survcomp)
library(cv)
library(biomaRt)


# setwd("/Put/Your/Working/Directory/Here")
load("UM20.rda")

#####################################################################################################
# Preprocessing
# UM data
head(protein_gene_cpm)
dim(protein_gene_cpm) # 14270   422

# Calculate the threshold for the number of samples (90% of total samples) in UM data
threshold <- 0.9 * ncol(protein_gene_cpm)

# Identify genes meeting both conditions:
# 1. Average CPM >= 10
# 2. 90% of samples with CPM >= 1
genes_to_keep <- rowMeans(protein_gene_cpm) >= 10 & rowSums(protein_gene_cpm >= 1) >= threshold

# Filter the dataframe to retain only the desired genes
filtered_genes_cpm <- protein_gene_cpm[genes_to_keep, , drop = FALSE]

# View the filtered UM dataframe
head(filtered_genes_cpm)
dim(filtered_genes_cpm) # 8589  422

# Calculate the threshold for the number of samples (90% of total samples) in external data
threshold2 <- 0.9 * ncol(cpm_matrix)

# Identify genes meeting both conditions:
# 1. Average RPKM >= 10
# 2. 90% of samples with RPKM >= 1
genes_to_keep2 <- rowMeans(cpm_matrix) >= 10 & rowSums(cpm_matrix >= 1) >= threshold2

# Filter the data frame to retain only the desired genes
ext_cpm_matrix <- cpm_matrix[genes_to_keep2, , drop = FALSE]

# View the filtered external dataframe
head(ext_cpm_matrix)
dim(ext_cpm_matrix) # 5721  132

# Match the row names of UM data with external data
matched_gene_list <- intersect(rownames(filtered_genes_cpm), rownames(ext_cpm_matrix))
length(matched_gene_list) # 5450 genes

# Filter for valid gene IDs
matching_ids <- matched_gene_list %in% colnames(exprs)
valid_gene_ids <- matched_gene_list[matching_ids] 
length(valid_gene_ids) # 5449

# exprs is the residual matrix based on UM case only, 422 samples, 18070 genes, 
# already filtered globin RNAs, and RNAs that only in one library type.
kept_gene_data <- exprs[, valid_gene_ids, drop = FALSE]
head(kept_gene_data)
dim(kept_gene_data) # 422 5449

# clinical data
head(case_df_stepwise)
dim(case_df_stepwise) # 420  32

# match gene expression data with clinical data
ordered_gene_matrix <- kept_gene_data[match(rownames(case_df_stepwise), rownames(kept_gene_data)), ]

head(ordered_gene_matrix)
dim(ordered_gene_matrix) # 420 5449

# only use the logrank selected genes
logrank_gene_list <- read_table('single_gene_logrank_with_pval0.01_cClinVar_575genes.txt')
logrank_gene_list <- logrank_gene_list$GeneID
length(logrank_gene_list) # 575

ordered_gene_matrix <- ordered_gene_matrix[, colnames(ordered_gene_matrix) %in% logrank_gene_list]

umdata <- cbind(case_df_stepwise, ordered_gene_matrix)
dim(umdata) # 420 607

#####################################################################################################
# xgboost model
umdata_xgboost <- umdata[,-c(1:8,10:12,14,16:32)]
datax_train <- as.matrix(umdata_xgboost)
datay_train <- setNames(umdata[, c("interval_time", "is_dead")], c("time", "status"))

# case only external data
external_data_xgboost <- external_dataset1[,match(colnames(datax_train), colnames(external_dataset1))]
datax_test <- as.matrix(external_data_xgboost)
# Convert character matrix to numeric
datax_test <- apply(datax_test, 2, as.numeric)
datay_test <- setNames(external_dataset1[, c("interval_time", "is_dead")], c("time", "status"))
datay_test$status <- ifelse(datay_test$status == "TRUE", 1, 0)

# Convert to DMatrix for XGBoost
time <- umdata$interval_time
y <- ifelse(umdata$is_dead == 1, time, -time) 
dtrain <- xgb.DMatrix(data = datax_train, label = y)

set.seed(314)
# survival object
surv_xgboost <- Surv(time = datay_train$time, event = datay_train$status)

# Grid search for hyperparameter tuning
search_grid <- expand.grid(
  max_depth = c(6, 7, 8),
  eta = c(0.01, 0.1),
  alpha = c(0.01, 0.1, 1),
  gamma = c(0, 0.1)
)

# 5 fold cross validation
folds <- createFolds(surv_xgboost[,1], k = 5, list = TRUE)

evaluate_model <- function(params) {
  cindices <- c()  
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(datax_train)), test_idx)
    
    dtrain <- xgb.DMatrix(data = datax_train[train_idx, ], label = surv_xgboost[train_idx, 1])
    dtest <- xgb.DMatrix(data = datax_train[test_idx, ], label = surv_xgboost[test_idx, 1])
    
    model <- xgb.train(
      data = dtrain,
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      max_depth = params$max_depth,
      eta = params$eta,
      alpha = params$alpha,
      gamma = params$gamma,
      nrounds = 100,
      verbose = 1
    )
    
    pred <- predict(model, dtest)
    cindices[i] <- concordance.index(x = pred, surv_xgboost[test_idx, 1], surv_xgboost[test_idx, 2])$c.index
  }
  
  return(mean(cindices))  
}

cv_results <- search_grid
cv_results$C_index <- NA  

for (i in 1:nrow(search_grid)) {
  params <- search_grid[i, ]
  cv_results$C_index[i] <- evaluate_model(params)
}

# find the best parameter
best_params <- cv_results[which.max(cv_results$C_index), ]
print(best_params)

best_params <- best_params %>%
  dplyr::select(-C_index)
  
best_params_list <- as.list(best_params[1,])

extra_param_list <- list(
  objective = "survival:cox",
  eval_metric = "cox-nloglik"
)

params <- c(best_params_list, extra_param_list)
print(params)

# Train a new XGBoost survival model with the best parameters
set.seed(314)

survival_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  verbose = 1
)

# Prepare the external test data as DMatrix
dtest <- xgb.DMatrix(data = datax_test)

# Predict risk scores using the trained model
risk_scores <- predict(survival_model, newdata = dtest)

# Evaluate model performance using concordance index (C-index)
cindex <- concordance.index(
  x = risk_scores,                       # Predicted risk scores
  surv.time = datay_test$time,           # Survival time
  surv.event = datay_test$status         # Survival event (dead/alive)
)

print(paste("Concordance Index:", cindex$c.index))

# Evaluate model performance using AUC 
tROC_xgboost_test <- timeROC(
  T = external_dataset1$interval_time,
  delta = as.logical(external_dataset1$is_dead),
  marker = risk_scores,
  cause = 1,
  weighting = "marginal",
  ROC = TRUE,
  times = c(1, 2, 3, 4, 5, 6, 7, 8),
  iid = TRUE
)

tROC_xgboost_df <-as.data.frame(tROC_xgboost_test[["AUC"]])

print(tROC_xgboost_df)

#####################################################################################################
# Feature selection
# get feature importance
importance_matrix <- xgb.importance(model = survival_model)
median_gain <- median(importance_matrix$Gain, na.rm = TRUE)
threshold_75 <- quantile(importance_matrix$Gain, probs = 0.75, na.rm = TRUE)
top_25_percent_features <- importance_matrix[importance_matrix$Gain > threshold_75, ]
top_significant_features  <- top_25_percent_features
top_significant_features  <- importance_matrix[importance_matrix$Gain > median_gain, ]

# top 10 features
top_10_feature <- importance_matrix[1:10,]
print(top_10_feature$Feature)

# top 20 features
top_20_feature <- importance_matrix[1:20,]
print(top_20_feature$Feature)

# top 50 features
top_50_feature <- importance_matrix[1:50,]
print(top_50_feature$Feature)

#####################################################################################################
# Use selected features to train new xgboost models and evaluate their performance on external data
xgboost_eval <- function(datax_train_in, datay_train_in, datax_test_in,
                        datay_test_in, features, parameter) {
  datax_train_sub <- datax_train_in[, features]
  datax_test_sub <- datax_test_in[, features]
  
  set.seed(314)
  # Convert to DMatrix for XGBoost
  time <- umdata$interval_time
  y <- ifelse(umdata$is_dead == 1, time, -time) 
  dtrain <- xgb.DMatrix(data = datax_train_sub, label = y)
  
  # Train the XGBoost survival model
  survival_model <- xgb.train(
    params = parameter,
    data = dtrain,
    nrounds = 100,
    verbose = 1
  )
  
  # Prepare the test data as DMatrix
  dtest <- xgb.DMatrix(data = datax_test_sub)
  
  # Predict risk scores using the trained model
  risk_scores <- predict(survival_model, newdata = dtest)
  
  # Evaluate model performance using concordance index (C-index)
  cindex <- concordance.index(
    x = risk_scores,                       # Predicted risk scores
    surv.time = datay_test_in$time,           # Survival time
    surv.event = datay_test_in$status         # Survival event (dead/alive)
  )
  
  print(paste("Concordance Index:", cindex$c.index))
  
  # Evaluate model performance using AUC 
  tROC_xgboost_test <- timeROC(
    T = external_dataset1$interval_time,
    delta = as.logical(external_dataset1$is_dead),
    marker = risk_scores,
    cause = 1,
    weighting = "marginal",
    ROC = TRUE,
    times = c(1, 2, 3, 4, 5, 6, 7, 8),
    iid = TRUE
  )
  
  tROC_xgboost_df <-as.data.frame(tROC_xgboost_test[["AUC"]])
  
  print(tROC_xgboost_df)
}

# top 10 features
xgboost_eval(datax_train, datay_train, datax_test, 
            datay_test, top_10_feature$Feature, params) 

# top 20 features
xgboost_eval(datax_train, datay_train, datax_test, 
             datay_test, top_20_feature$Feature, params) 

# top 50 features
xgboost_eval(datax_train, datay_train, datax_test, 
             datay_test, top_50_feature$Feature, params) 

#####################################################################################################
# Use selected features to train new CoxPH models and evaluate their performance on external data
coxph_eval <- function(features) {
  formula_xgboost <- as.formula(
    paste(
      "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment  + as.numeric(symptom_onset_date_age_at) + sex +",
      paste(features, collapse = " + ")
    )
  )
  model_xgboost <- coxph(formula_xgboost, data = umdata, method = "breslow")
  
  risk_scores <- predict(model_xgboost, newdata = external_dataset1, type = "risk")
  
  # Evaluate model performance using concordance index (C-index)
  cindex <- concordance.index(
    x = risk_scores,                       # Predicted risk scores
    surv.time = datay_test$time,           # Survival time
    surv.event = datay_test$status         # Survival event (dead/alive)
  )
  
  print(paste("Concordance Index:", cindex$c.index))
  
  # Evaluate model performance using AUC
  tROC_xgboost_test <- timeROC(
    T = external_dataset1$interval_time,
    delta = as.logical(external_dataset1$is_dead),
    marker = risk_scores,
    cause = 1,
    weighting = "marginal",
    ROC = TRUE,
    times = c(1, 2, 3, 4, 5, 6, 7, 8),
    iid = TRUE
  )
  
  tROC_xgboost_df <-as.data.frame(tROC_xgboost_test[["AUC"]])
  
  print(tROC_xgboost_df)
}

# Clinical features + the genes from the top 10 xgboost features
coxph_eval(top_10_feature$Feature[3:length(top_10_feature$Feature)])

# Clinical features + the genes from the top 20 xgboost features
coxph_eval(top_20_feature$Feature[3:length(top_20_feature$Feature)])

# Clinical features + the genes from the top 50 xgboost features
coxph_eval(top_50_feature$Feature[3:length(top_50_feature$Feature)])

#####################################################################################################
# convert the ensembl ids of the genes from the top10 xgboost features to gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_in_top10 <- top_10_feature$Feature[3:length(top_10_feature$Feature)]

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = genes_in_top10,
  mart = ensembl
)

print(gene_info)













