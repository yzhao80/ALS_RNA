# Title: ALS RNA-seq Survival and WGCNA Analysis
# Author: Bo Li, Yue Zhao
# Date: 2025-06-29
# Description: RNA-seq gene expression filtering, survival modeling, WGCNA, stepwise Cox regression, and external validation.
#              Also contains the codes for Figure 4b,c,d, e,f,g,h,i,j; Supplementary Figure S6, S7.
# ------------------------------
# 1. Setup 
# ------------------------------

library(caret)
library(leaps)
library(MASS)
library(qqman)
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(pROC)
library(timeROC)
library(survivalAnalysis)
library(finalfit)
library(biomaRt)
library(edgeR)
library(impute)
library(WGCNA)
library(gplots)
library(dplyr)
library(forcats)
library(grid)
library(gridExtra)
library(VennDiagram)
library(survcomp)
library(compareC)
library(rstatix)

setwd('/home/yuedz/ALS_RNA/')

# ------------------------------
# 2. Data Loading & Preprocessing 
# ------------------------------
## Load UM dataset (case_only)
load(file='residual_exprs_0527.rda')
exprs=residual_exprs

load(file='deg.case.control.sample694_nobatch.rda')
case_count <- deg.case.control.sample694.nobatch$df[,colnames(deg.case.control.sample694.nobatch$df) %in% rownames(exprs)]
dge_case <- DGEList(counts = case_count)

# Calculate CPM and add gene annotation
cpm_matrix_case <- cpm(dge_case)
rpkm <- (deg.case.control.sample694.nobatch$rpkm_genes)
rpkm_annotation <- rpkm[, (ncol(rpkm) - 3):ncol(rpkm)]
cpm_case_df <- as.data.frame(cpm_matrix_case)
cpm_case_df$geneid <- rownames(cpm_matrix_case)
cpm_merged_df <- cpm_case_df %>%
    dplyr::left_join(rpkm_annotation, by = "geneid")
rownames(cpm_merged_df) <- cpm_merged_df$geneid
cpm_merged_df$geneid <- NULL 

# Filter rows where the biotype column is equal to "protein_coding"
cpm_merged_df_filter <- cpm_merged_df[cpm_merged_df$biotype == "protein_coding", ]
protein_gene_cpm <- cpm_merged_df_filter[,1:(ncol(cpm_merged_df_filter) - 3)]

# Calculate the threshold for the number of samples (90% of total samples)
threshold <- 0.9 * ncol(protein_gene_cpm)

# Identify genes meeting both conditions:
# 1. Average CPM >= 10
# 2. 90% of samples with CPM >= 1
genes_to_keep <- rowMeans(protein_gene_cpm) >= 10 & rowSums(protein_gene_cpm >= 1) >= threshold
filtered_genes_cpm <- protein_gene_cpm[genes_to_keep, , drop = FALSE]

## Load external dataset
load("GSE234297_set.rda")
filtered_annot <- annot[annot$GeneType == "protein-coding" & annot$EnsemblGeneID != "", ]
matched_tbl <- tbl[rownames(tbl) %in% rownames(filtered_annot), ] 
rownames(matched_tbl) <- filtered_annot[rownames(matched_tbl), "EnsemblGeneID"]
dge <- DGEList(counts = matched_tbl)
cpm_matrix <- cpm(dge)

# Calculate the threshold for the number of samples (90% of total samples)
threshold2 <- 0.9 * ncol(cpm_matrix)

# Identify genes meeting both conditions:
# 1. Average RPKM >= 10
# 2. 90% of samples with RPKM >= 1
genes_to_keep2 <- rowMeans(cpm_matrix) >= 10 & rowSums(cpm_matrix >= 1) >= threshold2
ext_cpm_matrix <- cpm_matrix[genes_to_keep2, , drop = FALSE]
matched_gene_list <- intersect(rownames(filtered_genes_cpm), rownames(ext_cpm_matrix))
matching_ids <- matched_gene_list %in% colnames(exprs)
valid_gene_ids <- matched_gene_list[matching_ids] #5449
kept_gene_data <- exprs[, valid_gene_ids, drop = FALSE]

## Load clinical data
case_df_stepwise <- read.table("clinical_data_case.txt", header = TRUE, sep = "\t", row.names = 1)
case_df_stepwise$symptom_onset_date_age_at <- as.numeric(case_df_stepwise$symptom_onset_date_age_at)
case_df_stepwise$death_date_age_at <- as.numeric(case_df_stepwise$death_date_age_at)
case_df_stepwise$last_contact_date_age_at <- as.numeric(case_df_stepwise$last_contact_date_age_at)
case_df_stepwise$followup_time_from_symptom_onset_days <- as.numeric(case_df_stepwise$followup_time_from_symptom_onset_days)
case_df_stepwise$followup_time_from_symptom_onset_years <- case_df_stepwise$followup_time_from_symptom_onset_days/365
case_df_stepwise <- case_df_stepwise %>% mutate(interval_time = if_else(as.logical(is_dead), 
                                             as.numeric(death_date_age_at)-as.numeric(symptom_onset_date_age_at), 
                                             as.numeric(last_contact_date_age_at)-as.numeric(symptom_onset_date_age_at)))

# remove 2 samples with NA in interval time
case_df_stepwise <- case_df_stepwise %>% filter(!is.na(interval_time))
# format clinical data
case_df_stepwise$is_dead <- ifelse(case_df_stepwise$is_dead == "TRUE", 1, 0)
case_df_stepwise$sex <- ifelse(case_df_stepwise$sex == "M", 1, 0)
case_df_stepwise$onset_segment <- ifelse(case_df_stepwise$onset_segment == "Bulbar", 1, 0)
dim(case_df_stepwise)
ordered_gene_matrix <- kept_gene_data[match(rownames(case_df_stepwise), rownames(kept_gene_data)), ]

# ------------------------------
# 3. Train/Test Split 
# ------------------------------
## use 100% UM data as training
n_samples <- nrow(ordered_gene_matrix)
set.seed(314)
train_indices <- sample(1:n_samples, size = n_samples)
test_indices <- setdiff(1:n_samples, train_indices)

train_data_exp <- ordered_gene_matrix[train_indices,]
test_data_exp <- ordered_gene_matrix[test_indices,]

train_data_clinial <- case_df_stepwise[train_indices,]
test_data_clinial <- case_df_stepwise[test_indices,]

train_data <- cbind(train_data_clinial,train_data_exp)
test_data <- cbind(test_data_clinial,test_data_exp)

# ------------------------------
# 4. Log-rank Test
# ------------------------------
results_logrank <- data.frame(GeneID = character(),
    ChiSq = numeric(),
    pValue = numeric(),
    stringsAsFactors = FALSE)

for (gene_id in colnames(train_data_exp)) {
    gene_expression <- as.numeric(train_data_exp[, gene_id])
    median_expression <- median(gene_expression, na.rm = TRUE)
    train_data_clinial$gene_group <- ifelse(gene_expression > median_expression, "High", "Low")
    
    # Perform log-rank test
    surv_fit <- survdiff(Surv(interval_time, as.logical(is_dead)) ~ gene_group, data = train_data_clinial)
    chi_sq <- surv_fit$chisq
    p_value <- 1 - pchisq(chi_sq, df = 1)  
    results_logrank <- rbind(results_logrank, data.frame(
        GeneID = gene_id,
        ChiSq = chi_sq,
        pValue = p_value,
        stringsAsFactors = FALSE
    ))
}

significant_results <- results_logrank %>%
    filter(pValue <= 0.01) #575 genes 

# ------------------------------
# 5. WGCNA 
# ------------------------------
allowWGCNAThreads()
gene_names <- significant_results[,1]
wgcna_expr <- ordered_gene_matrix[, colnames(ordered_gene_matrix) %in% gene_names, drop = FALSE]

## detect the outlier samples
# pdf(file = "cluster_before_cutting_cox_cClinVar.pdf", wi = 10, he = 7)
sampleTree = hclust(dist(wgcna_expr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex = 0.2)
# dev.off()

## set cutHeight to remove outliers
clust <- cutreeStatic(sampleTree, cutHeight = 30, minSize=10)
keepSample <- (clust == 1) # no outlier detected for cox cClinVar
wgcna_expr_left <- wgcna_expr[keepSample, ]
# dim(wgcna_expr_left) # no outlier detected and removed

## Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_expr_left, powerVector = powers, verbose = 5)
# pdf(file = "softPowerSelection_rmOutliers_cox_cClinVar.pdf", wi = 10, he = 7)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()

sft[["fitIndices"]]
best_beta = 6; 
cor <- WGCNA::cor
## Co-expression networks constrction and Module dectection
net = blockwiseModules(wgcna_expr_left,
                       power = best_beta,
                       maxBlockSize = ncol(wgcna_expr_left),
                       TOMType = "unsigned",
                       minModuleSize = 5,
                       mergeCutHeight = 0.01,
                       deepSplit = 4, 
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "case_rm-TOM",
                       verbose = 3,
                       corType = "pearson",
                       networkType = "unsigned")
cor<-stats::cor

moduleColors = labels2colors(net$colors)
# pdf(file = "gene_clusters_allsamples_UM100_cClinVar.pdf", wi = 10, he = 7)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05)
table(net$colors)
# dev.off()

# Calculate eigengenes for all modules
moduleEigengenes <- moduleEigengenes(wgcna_expr_left, net$colors)

## Based on eigengene to identify hub genes in each module 
MEs <- moduleEigengenes(wgcna_expr_left, colors = net$colors)$eigengenes
# Calculate module membership (kME) for each gene in each module
moduleMembership <- as.data.frame(cor(wgcna_expr_left, MEs, use = "p"))
rownames(moduleMembership) <- colnames(wgcna_expr_left)
moduleMembership$module <- net$colors

hubGenes <- list()
# Set the top 10% genes in each module as hub genes
topPercent <- 0.1
moduleLabels <- unique(net$colors)
for (module in moduleLabels) {
  eigengene_name <- paste0("ME", module)
  moduleGenes <- names(net$colors[net$colors == module])
  kME <- moduleMembership[moduleGenes, eigengene_name, drop = FALSE]
  
  # Sort genes by their kME in descending order
  kME_sorted <- kME[order(-kME[,1]), , drop = FALSE]

  # Special case for module 0: output all genes
  if (module == 0) {
    hubGenes[[eigengene_name]] <- rownames(kME_sorted)
  } else {
    # For other modules, determine the number of top genes to consider as hub genes
    numTopGenes <- ceiling(nrow(kME_sorted) * topPercent)
    
    # Ensure at least one gene is selected as a hub gene
    if (numTopGenes < 1) {
      numTopGenes <- 1
    }
    hubGenes[[eigengene_name]] <- rownames(kME_sorted)[1:numTopGenes]
  }
}

# Output the hub genes for all modules
# for (module in names(hubGenes)) {
#   cat("Module:", module, "\n")
#   cat("Hub Genes:", hubGenes[[module]], "\n\n")
# }

# Save hub genes to a file (all genes for module 0, top hub genes for others)
hubGenesDF <- do.call(rbind, lapply(names(hubGenes), function(module) {
  data.frame(module = module, gene = hubGenes[[module]])
}))

wgcna_df <- hubGenesDF

# ------------------------------
# 6. Training the Stepwise model with MASS
# ------------------------------
wgcna_gene <- as.character(wgcna_df[[2]])

# Ensure only existing gene IDs are used
existing_wgcna_gene <- intersect(wgcna_gene, colnames(train_data_exp))
existing_wgcna_gene_df<-data.frame(geneID=existing_wgcna_gene)
wgcna_gene_expr <- train_data_exp[,existing_wgcna_gene,drop = FALSE]

## combine expr data with clinical data
wgcna_gene_expr_order <- wgcna_gene_expr[rownames(train_data_clinial), ]
stepwise_combined_data <- cbind(train_data_clinial,wgcna_gene_expr_order)

## try both forward and backward stepwise selection with MASS library
base_formula <- as.formula("Surv(interval_time, as.logical(is_dead)) ~ onset_segment + 
                           as.numeric(symptom_onset_date_age_at) + sex")
initial_model <- coxph(base_formula, data = stepwise_combined_data)

# Define the full model formula with all potential covariates (genes)
full_formula <- as.formula(paste("Surv(interval_time, as.logical(is_dead)) ~ onset_segment + 
                                  as.numeric(symptom_onset_date_age_at) + sex +", 
                                  paste(existing_wgcna_gene, collapse = " + ")))

# Perform forward stepwise selection
mass_model <- stepAIC(initial_model, scope = list(lower = base_formula, upper = full_formula), trace = FALSE, direction = "both")
selected_genes <- names(coef(mass_model)) #21

# Filter out the non-gene variables (onset_segment, symptom_onset_date_age_at, sex)
selected_gene_names <- selected_genes[selected_genes %in% existing_wgcna_gene] # 18 genes
model_summary <- capture.output(summary(mass_model))
feature_names <- names(coef(mass_model))
gene_symbols <- rpkm_annotation$gene_symbol[match(feature_names, rpkm_annotation$geneid)]
hazard_ratios <- exp(coef(mass_model))
p_values <- summary(mass_model)$coefficients[, "Pr(>|z|)"]
log_p_values <- -10 * log10(p_values)
result_df <- data.frame(
  Ensembl_ID = feature_names,
  Gene_Symbol = gene_symbols,
  Hazard_Ratio = hazard_ratios,
  p_Value = p_values,
  neg10Log10_p_Value = log_p_values
)

result_df$Gene_Symbol[1:3] <- c("onset segment", "onset age","sex")
order_column <- "neg10Log10_p_Value"
stepwise_features <- result_df %>%
    mutate(Gene_Symbol = fct_reorder(Gene_Symbol, !!sym(order_column), .desc = FALSE))

# ------------------------------
# 7. Figure 4b & c
# ------------------------------
# Create the orange barplot in Fig 4b
# pdf("Stepwise.CoxPH.barplot.orange.pdf", height = 7, width = 6)
ggplot(stepwise_features, aes(x = Gene_Symbol, y = !!sym(order_column), fill = Hazard_Ratio)) + 
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = c("#FCBE80", "#BC3C29")) +  
    scale_y_continuous(limits = c(0, 60)) +
    labs(x = "", y = "", fill = "Hazard ratio") +
    theme_classic() +
    coord_flip() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
# dev.off()

# Create the blue barplot in Fig 4b
# pdf("Stepwise.CoxPH.barplot.blue.pdf", height = 7, width = 6)
ggplot(stepwise_features, aes(x = Gene_Symbol, y = !!sym(order_column), fill = Hazard_Ratio)) + 
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = c("#0072B5", "white")) +  
    scale_y_continuous(limits = c(0, 60)) +
    labs(x = "", y = "", fill = "Hazard ratio") +
    theme_classic() +
    coord_flip() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
# dev.off()

# ------------------------------
# 8. External Validation
# ------------------------------
stepwise_gene_list <- selected_gene_names
# stepwise_gene_list <- c(
# "ENSG00000118564", "ENSG00000173914", "ENSG00000153551", "ENSG00000067182",
# "ENSG00000126903", "ENSG00000057704", "ENSG00000153310", "ENSG00000104960",
# "ENSG00000106603", "ENSG00000179933", "ENSG00000104915", "ENSG00000159128",
# "ENSG00000215114", "ENSG00000008838", "ENSG00000169180", "ENSG00000169891",
# "ENSG00000116906", "ENSG00000108848")

xgboost_gene_list <- c("ENSG00000169994", "ENSG00000187837", "ENSG00000082014", "ENSG00000087088", "ENSG00000102710", "ENSG00000103154", "ENSG00000150990", "ENSG00000150527")

# Baseline model
formula_baseline <- as.formula(
  "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment  + as.numeric(symptom_onset_date_age_at) + sex"
)
model_baseline <- coxph(formula_baseline, data = train_data, method = "breslow")

# Stepwise model
formula_stepwise <- as.formula(
  paste(
    "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment  + as.numeric(symptom_onset_date_age_at) + sex +",
    paste(stepwise_gene_list, collapse = " + ")
  )
)
model_stepwise <- coxph(formula_stepwise, data = train_data, method = "breslow")

# Xgboost selected features testing on CoxPH
formula_xgboost <- as.formula(
    paste(
        "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment  + as.numeric(symptom_onset_date_age_at) + sex +",
        paste(xgboost_gene_list, collapse = " + ")
    )
)
model_xgboost <- coxph(formula_xgboost, data = train_data, method = "breslow")

# Load external dataset1 (Grima et al.) and evaluate risk scores
load('external_dataset1.rda')
external_dataset1$is_dead <- as.logical(external_dataset1$is_dead)
risk_scores_baseline_test <- predict(model_baseline, newdata = external_dataset1, type = "risk")
risk_scores_stepwise_test <- predict(model_stepwise, newdata = external_dataset1, type = "risk")
risk_scores_xgboostCox_test <- predict(model_xgboost, newdata = external_dataset1, type = "risk")

# Time-dependent ROC curves
tROC_baseline_test <- timeROC(
    T = external_dataset1$interval_time,
    delta = as.logical(external_dataset1$is_dead),
    marker = risk_scores_baseline_test,
    cause = 1,
    weighting = "marginal",
    ROC = TRUE,
    times = c(1, 2, 3, 4, 5, 6, 7, 8),
    iid = TRUE
)

tROC_stepwise_test <- timeROC(
  T = external_dataset1$interval_time,
  delta = as.logical(external_dataset1$is_dead),
  marker = risk_scores_stepwise_test,
  cause = 1,
  weighting = "marginal",
  ROC = TRUE,
  times = c(1, 2, 3, 4, 5, 6, 7, 8),
  iid = TRUE
)

tROC_xgboostCox_test <- timeROC(
    T = external_dataset1$interval_time,
    delta = as.logical(external_dataset1$is_dead),
    marker = risk_scores_xgboostCox_test,
    cause = 1,
    weighting = "marginal",
    ROC = TRUE,
    times = c(1, 2, 3, 4, 5, 6, 7, 8),
    iid = TRUE
)

# ------------------------------
# 9. tROC curves, Fig 4e
# ------------------------------
# pdf("ThreeModels.4yrs.external.pdf")
plot(
  tROC_baseline_test, time = 4, add = FALSE, col = "grey50",
  title = FALSE, lwd = 3, cex.lab = 1.5, cex.axis = 1.2
)
plot(tROC_stepwise_test, time = 4, add = TRUE, col = "#548235", lwd = 3)
plot(tROC_xgboostCox_test, time = 4, add = TRUE, col = "#B21239", lwd = 3)
legend(
  "bottomright",
  legend = c(
    paste("4-year Baseline AUC =", round(tROC_baseline_test$AUC[4],3)),
    paste("4-year Stepwise AUC =", round(tROC_stepwise_test$AUC[4],3)),
    paste("4-year XGBoost AUC =", round(tROC_xgboostCox_test$AUC[4],3))
  ),
  col = c("grey50", "#548235","#B21239"), lwd = 3, cex = 1.5
)
# dev.off()

tROC_AUC_baseline_df<-as.data.frame(tROC_baseline_test[["AUC"]])
tROC_AUC_stepwise_df<-as.data.frame(tROC_stepwise_test[["AUC"]])
tROC_AUC_xgboostCox_df<-as.data.frame(tROC_xgboostCox_test[["AUC"]])
tROC_AUC_df <- cbind.data.frame(tROC_AUC_baseline_df,tROC_AUC_stepwise_df,tROC_AUC_xgboostCox_df)
# tROC_AUC_df

external_dataset1 <- external_dataset1[, colnames(external_dataset1) != ""]

# ------------------------------------------------------------------------------------------
# 10. KM Plot, log-rank test, Fig 4f, g, h; Pairwise log-rank test, supplementary Fig S6c
# ------------------------------------------------------------------------------------------
KM3_plot <- function(data, multicox_model){
    # Predict risk scores
    risk_scores <- predict(multicox_model, newdata = data, type = "risk")
    
    # Assign risk scores to data
    data$risk_score <- risk_scores
    
    # Define quantiles for risk stratification
    quantile_25 <- quantile(risk_scores, 0.25)
    quantile_75 <- quantile(risk_scores, 0.75)
    
    # Assign risk groups based on quantiles
    data$risk_group <- ifelse(data$risk_score > quantile_75, "High", 
        ifelse(data$risk_score < quantile_25, "Low", "Middle"))
    data$risk_group <- factor(data$risk_group, levels = c("High", "Middle", "Low"))
    
    # Fit Kaplan-Meier survival curves
    km_fit <- survfit(Surv(interval_time, is_dead) ~ risk_group, data = data)
    
    # Pairwise log-rank tests
    pw <- pairwise_survdiff(
        Surv(interval_time, is_dead) ~ risk_group,
        data             = data,
        p.adjust.method  = "bonferroni"
    )
    print(pw)
    
    # Calculate median survival times
    median_surv_time <- summary(km_fit)$table[, "median"]
    diff_high_low <- median_surv_time[3] - median_surv_time[1]
    
    message1 <- sprintf("Diff in median survival between high and low risk groups is %f", diff_high_low)
    print(message1)
    print(median_surv_time)
    print(quantile_25)
    print(quantile_75)
    
    # Generate Kaplan-Meier plot
    ggsurvplot(km_fit, data = data, pval = TRUE, conf.int = FALSE, risk.table = TRUE, surv.median.line = "hv",
        palette = c("#E41A1C", "#377EB8", "#4DAF4A"), # Red, Green, Blue for High, Middle, Low
        xlim = c(0, 8), break.time.by = 1, main = "Kaplan-Meier Curve", 
        xlab = "Time from Diagnosis (Years)", ylab = "Survival Probability",
        font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), font.tickslab = c(20, "black"),
        font.legend = c(20, "bold", "black"), tables.col = "strata", fontsize = 5, risk.table.fontsize = 5)   
}

# Fig 4f, g, h
pdf("01_stepwise/baseline.KM3plot.pdf")
KM3_plot(external_dataset1, model_baseline)
dev.off()
pdf("01_stepwise/stepwise.KM3plot.pdf")
KM3_plot(external_dataset1, model_stepwise)
dev.off()
pdf("03_xgboost/xgboostCox.KM3plot.pdf")
KM3_plot(external_dataset1, model_xgboost)
dev.off()

# ------------------------------
# 11. Bar Plot, Fig 4c
# ------------------------------
# Barplot for XGBoost selected features, Fig 4c
feature_names <- names(coef(model_xgboost))
gene_symbols <- rpkm_annotation$gene_symbol[match(feature_names, rpkm_annotation$geneid)]
hazard_ratios <- exp(coef(model_xgboost))
p_values <- summary(model_xgboost)$coefficients[, "Pr(>|z|)"]
log_p_values <- -10 * log10(p_values)
result_df <- data.frame(
    Ensembl_ID = feature_names,          # Keep original Ensembl IDs
    Gene_Symbol = gene_symbols,          # Gene symbols (NA if not available)
    Hazard_Ratio = hazard_ratios,        # Hazard ratios
    p_Value = p_values,
    neg10Log10_p_Value = log_p_values               # -10 * log10(p-values)
)
result_df$Gene_Symbol[1:3] <- c("onset segment", "onset age","sex")
order_column <- "neg10Log10_p_Value"
xgboost_features <- result_df %>%
    mutate(Gene_Symbol = fct_reorder(Gene_Symbol, !!sym(order_column), .desc = FALSE))

# Create the orange barplot
pdf("03_xgboost/xgboost.CoxPH.barplot.orange.pdf", height = 5, width = 6)
ggplot(xgboost_features, aes(x = Gene_Symbol, y = !!sym(order_column), fill = Hazard_Ratio)) + 
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = c("#FCBE80", "#BC3C29")) +  # Orange gradient
    scale_y_continuous(limits = c(0, 60)) +
    labs(x = "", y = "", fill = "Hazard ratio") +
    theme_classic() +
    coord_flip() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()

# Create the blue barplot
pdf("03_xgboost/xgboost.CoxPH.barplot.blue.pdf", height = 5, width = 6)
ggplot(xgboost_features, aes(x = Gene_Symbol, y = !!sym(order_column), fill = Hazard_Ratio)) + 
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = c("#0072B5", "white")) +  # Blue gradient
    scale_y_continuous(limits = c(0, 60)) +
    labs(x = "", y = "", fill = "Hazard ratio") +
    theme_classic() +
    coord_flip() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()

# ------------------------------
# 12. Fig 4i, j
# ------------------------------
# Venn diagram for high and low risk groups, Fig 4i, 4j
models <- list(
    "Baseline" = model_baseline,
    "Stepwise" = model_stepwise,
    "XGBoost" = model_xgboost
)

high_risk_groups <- list()

for (model_name in names(models)) {
    risk_scores <- predict(models[[model_name]], newdata = external_dataset1, type = "risk")
    quantile_75 <- quantile(risk_scores, 0.75)
    high_risk <- rownames(external_dataset1)[risk_scores > quantile_75]
    high_risk_groups[[model_name]] <- high_risk
}

pdf("Fast_progressor_venn.pdf")
if (length(models) <= 3) {
    grid.newpage()  
    venn.plot <- venn.diagram(
        x = high_risk_groups,
        category.names = names(models),
        filename = NULL, 
        output = TRUE,
        #   fill = c("red", "green", "blue", "yellow"),
        fill = c("grey50","#548235","#B21239"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        main = "High-Risk Group Consistency Across Models"
    )
    grid.draw(venn.plot)
}
dev.off()

low_risk_groups <- list()
for (model_name in names(models)) {
    risk_scores <- predict(models[[model_name]], newdata = external_dataset1, type = "risk")
    quantile_25 <- quantile(risk_scores, 0.25)
    low_risk <- rownames(external_dataset1)[risk_scores < quantile_25]
    low_risk_groups[[model_name]] <- low_risk
}

pdf("Slow_progressor_venn.pdf")
if (length(models) <= 3) {
    grid.newpage()  
    venn.plot <- venn.diagram(
        x = low_risk_groups,
        category.names = names(models),
        filename = NULL, 
        output = TRUE,
        #   fill = c("red", "green", "blue", "yellow"),
        fill = c("grey50","#548235", "#B21239"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        main = "Low-Risk Group Consistency Across Models"
    )
    grid.draw(venn.plot)
}
dev.off()

# ------------------------------
# 13. C-index and Supp Fig. S6c
# ------------------------------
c_index_baseline <- concordance.index(
    x = risk_scores_baseline_test,
    surv.time = external_dataset1$interval_time,
    surv.event = external_dataset1$is_dead
)

c_index_stepwise <- concordance.index(
    x = risk_scores_stepwise_test,
    surv.time = external_dataset1$interval_time,
    surv.event = external_dataset1$is_dead
)

c_index_xgboostCox <- concordance.index(
    x = risk_scores_xgboostCox_test,
    surv.time = external_dataset1$interval_time,
    surv.event = external_dataset1$is_dead
)

c_index_baseline$c.index # 0.659854
c_index_stepwise$c.index # 0.6508168
c_index_xgboostCox$c.index # 0.6922404

result <- compareC(
    timeX   = external_dataset1$interval_time,
    statusX = external_dataset1$is_dead,
    scoreY = risk_scores_baseline_test,
    scoreZ = risk_scores_xgboostCox_test
)

result2 <- compareC(
    timeX   = external_dataset1$interval_time,
    statusX = external_dataset1$is_dead,
    scoreY = risk_scores_baseline_test,
    scoreZ = risk_scores_stepwise_test
)

result$pval # 0.05057299
result2$pval # 0.794869


# ------------------------------
# 14. Supp Fig. S7: ROC curves over multiple time points
# ------------------------------
# Function to plot ROC curves over multiple time points
plot_roc_curves <- function(roc_models, time_range, file_name, pdf_width=16, pdf_height = 9, roc_colors) {
    if (length(roc_colors) < length(roc_models)) {
        stop("Number of colors must match the number of models.")
    }
    
    pdf(file_name, width = pdf_width, height = pdf_height)
    par(mfrow = c(2, 4))
    
    for (t in time_range) {
        first <- TRUE  # Track first plot in each time iteration
        for (i in seq_along(roc_models)) {
            plot(
                roc_models[[i]], time = t, add = !first, col = roc_colors[i],
                lwd = 3, cex.lab = 1.5, cex.axis = 1.2, title = FALSE
            )
            first <- FALSE
        }
        
        legend("bottomright",
            legend = sapply(seq_along(roc_models), function(i)
                paste0(t, " Years ", names(roc_models)[i], " AUC = ", round(roc_models[[i]]$AUC[t], 3))
            ),
            col = roc_colors, lwd = 3, cex = 1.2
        )
    }
    dev.off()
}

roc_models <- list(
    "Baseline" = tROC_baseline_test,
    "Stepwise" = tROC_stepwise_test,
    "XGBoost" = tROC_xgboostCox_test
)

roc_colors <- c("grey50", "#548235", "#B21239")

plot_roc_curves(
    roc_models = roc_models,
    time_range = 2:8,
    file_name = "Models.2-8Years.External.pdf",
    roc_colors = roc_colors
)


# ------------------------------------------------------------
# 15. Fig. 4d: 30 random splits, bar plot
# ------------------------------------------------------------
seeds <- c(
    111, 222, 333, 444, 555, 666, 777, 888, 999,
    1825, 410, 4507, 4013, 3658, 2287, 1680,
    8936, 1425, 9675, 6913, 521, 489, 1536,
    3583, 3812, 8280, 9864, 435, 9196, 3258
)

roc_times <- c(1, 2, 3, 4, 5, 6, 7, 8)

auc_results <- data.frame()

for (i in 1:30) {
    set.seed(seeds[i])
    
    n_samples <- nrow(ordered_gene_matrix)
    train_indices <- sample(1:n_samples, size = n_samples * 0.8)
    test_indices <- setdiff(1:n_samples, train_indices)
    
    train_data_exp <- ordered_gene_matrix[train_indices,]
    test_data_exp <- ordered_gene_matrix[test_indices,]
    
    train_data_clinical <- case_df_stepwise[train_indices,]
    test_data_clinical <- case_df_stepwise[test_indices,]
    
    train_data <- cbind(train_data_clinical, train_data_exp)
    test_data <- cbind(test_data_clinical, test_data_exp)
    
    # Baseline Model
    formula_baseline <- as.formula(
        "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment + as.numeric(symptom_onset_date_age_at) + sex"
    )
    model_baseline <- coxph(formula_baseline, data = train_data, method = "breslow")
    
    # Clinical + Stepwise Genes
    formula_stepwise <- as.formula(
        paste(
            "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment + as.numeric(symptom_onset_date_age_at) + sex +",
            paste(stepwise_gene_list, collapse = " + ")
        )
    )
    model_stepwise <- coxph(formula_stepwise, data = train_data, method = "breslow")
    
    # Clinical + XGBoost Genes
    formula_xgboostCox <- as.formula(
        paste(
            "Surv(as.numeric(interval_time), as.logical(is_dead)) ~ onset_segment + as.numeric(symptom_onset_date_age_at) + sex +",
            paste(xgboost_gene_list, collapse = " + ")
        )
    )
    model_xgboostCox <- coxph(formula_xgboostCox, data = train_data, method = "breslow")
    
    # Genes-only Stepwise
    formula_genes_only_stepwise <- as.formula(
        paste(
            "Surv(as.numeric(interval_time), as.logical(is_dead)) ~",
            paste(stepwise_gene_list, collapse = " + ")
        )
    )
    model_genes_only_stepwise <- coxph(formula_genes_only_stepwise, data = train_data, method = "breslow")
    
    # Genes-only XGBoost
    formula_genes_only_xgboost <- as.formula(
        paste(
            "Surv(as.numeric(interval_time), as.logical(is_dead)) ~",
            paste(xgboost_gene_list, collapse = " + ")
        )
    )
    model_genes_only_xgboost <- coxph(formula_genes_only_xgboost, data = train_data, method = "breslow")
    
    # Risk Scores
    rs_baseline <- predict(model_baseline, newdata = test_data, type = "risk")
    rs_stepwise <- predict(model_stepwise, newdata = test_data, type = "risk")
    rs_xgboost <- predict(model_xgboostCox, newdata = test_data, type = "risk")
    rs_genes_only_stepwise <- predict(model_genes_only_stepwise, newdata = test_data, type = "risk")
    rs_genes_only_xgboost <- predict(model_genes_only_xgboost, newdata = test_data, type = "risk")
    
    # Time-dependent ROC
    tROC_baseline <- timeROC(test_data$interval_time, test_data$is_dead, rs_baseline, cause=1, times=roc_times)
    tROC_stepwise <- timeROC(test_data$interval_time, test_data$is_dead, rs_stepwise, cause=1, times=roc_times)
    tROC_xgboost <- timeROC(test_data$interval_time, test_data$is_dead, rs_xgboost, cause=1, times=roc_times)
    tROC_genes_only_stepwise <- timeROC(test_data$interval_time, test_data$is_dead, rs_genes_only_stepwise, cause=1, times=roc_times)
    tROC_genes_only_xgboost <- timeROC(test_data$interval_time, test_data$is_dead, rs_genes_only_xgboost, cause=1, times=roc_times)
    
    auc_results <- rbind(
        auc_results,
        data.frame(
            Seed = seeds[i],
            Time = roc_times,
            Baseline = tROC_baseline$AUC,
            Clinical_Stepwise = tROC_stepwise$AUC,
            Clinical_XGBoost = tROC_xgboost$AUC,
            Genes_Only_Stepwise = tROC_genes_only_stepwise$AUC,
            Genes_Only_XGBoost = tROC_genes_only_xgboost$AUC
        )
    )
}

auc_summary <- auc_results %>%
    group_by(Time) %>%
    summarise(across(Baseline:Genes_Only_XGBoost, list(mean=mean, sd=sd)))

# write.table(auc_results, "01_stepwise/AUC_results_5models_30seeds.txt", sep="\t", quote=FALSE, row.names=FALSE)
# write.table(auc_summary, "01_stepwise/AUC_summary_5models_30seeds.txt", sep="\t", quote=FALSE, row.names=FALSE)

df <- auc_summary %>%
    filter(!is.na(Baseline_mean))

df_long <- df %>%
    pivot_longer(
        cols = -Time,
        names_to = c("Method", "Measure"),
        names_pattern = "(.*)_(mean|sd)",
        values_to = "Value"
    ) %>%
    pivot_wider(
        names_from = Measure,
        values_from = Value
    ) %>%
    dplyr::rename(
        MeanAUC = mean,
        SDAUC    = sd
    )

df_long <- df_long %>%
    mutate(
        Method = recode(
            Method,
            "Baseline"               = "Baseline",
            "Clinical_Stepwise"      = "Clinical+Stepwise",
            "Clinical_XGBoost"       = "Clinical+XGBoost",
            "Genes_Only_Stepwise"    = "Genes-only Stepwise",
            "Genes_Only_XGBoost"     = "Genes-only XGBoost"
        )
    )

method_colors <- c(
    "Baseline"               = "gray50",
    "Clinical+Stepwise"      = "#823C96",
    "Clinical+XGBoost"       = "#D78296",
    "Genes-only Stepwise"    = "#4C772B",
    "Genes-only XGBoost"     = "#A7D486"
)

# Bar plot for 5 models
pdf("internal_5Random_bar.pdf", height = 3.5, width = 8)

ggplot(df_long, aes(x = factor(Time), y = MeanAUC, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
        aes(ymin = MeanAUC - SDAUC, ymax = MeanAUC + SDAUC),
        width = 0.2,
        position = position_dodge(width = 0.8)
    ) +
    scale_fill_manual(values = method_colors) +  
    labs(
        title = "Mean AUC Over Time (with SD Error Bars)",
        x = "Year",
        y = "Mean AUC",
        fill = "Model"
    ) +
    scale_x_discrete(labels = sort(unique(df_bar$Time))) +
    theme_minimal() +
    theme(
        text = element_text(size = 14),
        legend.position = "right"
    )

dev.off()


# -------------------------------------------------------------------------
# 16. Supp Fig. S6a: Compare difference between three models in UM data
# -------------------------------------------------------------------------
auc_long <- auc_results %>%  
    dplyr::select(Seed, Time, Baseline, Clinical_Stepwise, Clinical_XGBoost) %>%
    pivot_longer(cols = c(Baseline, Clinical_Stepwise, Clinical_XGBoost),
        names_to  = "Model",
        values_to = "AUC") %>% 
    filter(!is.na(AUC))                 

auc_long$Model <- factor(auc_long$Model, levels=c("Clinical_XGBoost","Clinical_Stepwise","Baseline"))

per_time_results <- auc_long %>% 
    group_by(Time) %>% 
    pairwise_wilcox_test(
        AUC ~ Model,
        paired  = TRUE,
        p.adjust.method = "bonferroni",
        detailed = TRUE
    )

per_time_wide <- per_time_results %>%
    mutate(
        Comparison = case_when(
            (group1 == "Baseline" & group2 == "Clinical_Stepwise") |
                (group2 == "Baseline" & group1 == "Clinical_Stepwise") ~ "Stepwise vs Clinical",
            (group1 == "Baseline" & group2 == "Clinical_XGBoost") |
                (group2 == "Baseline" & group1 == "Clinical_XGBoost") ~ "XGBoost vs Clinical",
            (group1 == "Clinical_Stepwise" & group2 == "Clinical_XGBoost") |
                (group2 == "Clinical_Stepwise" & group1 == "Clinical_XGBoost") ~ "XGBoost vs Stepwise",
            TRUE ~ NA_character_
        ),
        `ΔAUC` = round(estimate, 3),
        `95% CI` = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
        `adj. p` = signif(p.adj, 2)
    ) %>%
    select(Time, Comparison, `ΔAUC`, `95% CI`, `adj. p`) %>%
    pivot_wider(
        names_from = Comparison,
        values_from = c(`ΔAUC`, `95% CI`, `adj. p`)
    ) %>%
    arrange(Time)

per_time_wide <- per_time_wide %>%
    select(
        Time,
        `ΔAUC_Stepwise vs Clinical`,
        `95% CI_Stepwise vs Clinical`,
        `adj. p_Stepwise vs Clinical`,
        `ΔAUC_XGBoost vs Clinical`,
        `95% CI_XGBoost vs Clinical`,
        `adj. p_XGBoost vs Clinical`,
        `ΔAUC_XGBoost vs Stepwise`,
        `95% CI_XGBoost vs Stepwise`,
        `adj. p_XGBoost vs Stepwise`
    )

write.table(per_time_wide, "AUC_p_values_3models_30seeds.txt", sep="\t", quote=FALSE, row.names=FALSE)