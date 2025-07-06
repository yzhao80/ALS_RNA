###############################################################################
# TDP-43 and postmortem spinal cord (PM) and blood core genes drug connectivity analysis 
# Author: Kai Guo
# Date: 2024-12-19
#
# Key analyses include:
# - Differential expression filtering
# - Drug connectivity mapping (CMAP/LINCS)
# - Venn diagram comparisons
# - Network visualization of drug-target interactions
# - Mechanism of action (MOA) analysis
#
# Input files required:
# - neuron_spinal_cord.rda (DEG data)
# - Various .gct files from CMAP/LINCS
# - log2fc_df.rda (fold change data)
###############################################################################

# --------------------------
# 1. INITIAL SETUP
# --------------------------

# Set working directory to project folder
setwd("Desktop/Project/ALS_transcriptomics/new/")

# Load required libraries
library(tidyverse)    # Data manipulation and visualization
library(cmapR)        # For parsing GCT files (CMAP/LINCS data)
library(igraph)       # Network analysis and visualization
library(GGally)       # Network visualization
library(venndetail)   # Advanced Venn diagrams

# --------------------------
# 2. LOAD AND FILTER DIFFERENTIALLY EXPRESSED GENES (DEGs)
# --------------------------

# Load pre-processed DEG data
load('neuron_spinal_cord.rda')

# Filter significant DEGs for TDP-43 condition:
# - Adjusted p-value < 0.01
# - Absolute log2 fold change > 0.1
tdp <- subset(TDP43_deg_sig_strong, padj < 0.01 & abs(log2FoldChange) > 0.1)

# Filter significant DEGs for PM condition:
# - FDR < 0.01
# - Absolute log2 fold change > 0.1
pm <- subset(PM_DEG_filter_strong, FDR < 0.01 & abs(log2FoldChange) > 0.1)

# Save filtered DEG lists sorted by fold change
tdp %>% arrange(desc(log2FoldChange)) %>% write.csv(file = "tdp_deg.csv")
pm %>% arrange(desc(log2FoldChange)) %>% write.csv(file = "pm_deg.csv")

# Load and save log2FC data for adjusted and non-adjusted analyses
load('log2fc_df.rda')
write.csv(log2fc_df_yesAdj, file = "log2fc_df_yesAdj.csv")
write.csv(log2fc_df_noAdj, file = "log2fc_df_noAdj.csv")

# --------------------------
# 3. DRUG CONNECTIVITY ANALYSIS - TDP-43 CONDITION
# --------------------------

# Parse drug connectivity data from CMAP/LINCS
drug <- parse_gctx("TDP/arfs/TAG/query_result.gct")

# Prepare drug data matrix:
# - drug_m: Connectivity scores matrix
# - drug_d: Drug metadata (row descriptors)
drug_m <- as.data.frame(drug@mat)
drug_d <- drug@rdesc
dd <- cbind(drug_d, drug_m[rownames(drug_d),])
colnames(dd)[ncol(dd)] <- "norm_cs"  # Normalized connectivity score

# Filter drug data with strict criteria:
# - Moderate connectivity (-2 < norm_cs < 2)
# - Treatment type is compound (not genetic/other)
# - Significant FDR (fdr_q_nlog10 > 2 => q < 0.01)
# - Quality control passed
# - Mechanism of action (MOA) available
drug_ds <- dd %>%
    filter(abs(norm_cs) < 2,
        pert_type == "trt_cp",
        fdr_q_nlog10 > 2,
        qc_pass == 1,
        moa != -666)

# Convert dose to numeric (removing 'uM' unit)
drug_ds$dose <- as.numeric(sub(' uM', '', drug_ds$pert_idose))

# Generate diagnostic plots:

# 1. Treatment time distribution
barplot(table(drug_ds$pert_itime), las = 2, col = distcolor[8],
    main = "Treatment Time Distribution",
    xlab = "Time", ylab = "Number of Treatments")
dev.print(pdf, file = "TDP_select_all_time_bar.pdf")

# 2. Signature strength distribution
hist(drug_ds$ss_ngene, breaks = 100, col = distcolor[9],
    border = "white", xlab = "Signature Strength",
    main = "Drug Signature Strength Distribution")
dev.print(pdf, file = "TDP_Signature_Strength_all_hist.pdf")

# 3. Top 50 MOA categories
drug_ds %>%
    select(pert_iname, moa) %>%
    distinct() %>%
    group_by(moa) %>%
    summarise(num = n()) %>%
    arrange(desc(num)) %>%
    head(n = 50) %>%
    mutate(moa = ft_sort(moa, num)) %>%
    ggplot(aes(moa, num)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_classic() +
    labs(title = "Top 50 Mechanism of Action Categories",
        x = "MOA", y = "Number of Drugs")
dev.print(pdf, file = "TDP_moa_all_number.pdf")

# Save filtered TDP drug data
tdpd <- drug_ds

# --------------------------
# 4. DRUG CONNECTIVITY ANALYSIS - PM CONDITION
# --------------------------

# Repeat similar analysis for PM condition
drug <- parse_gctx("PM/arfs/TAG/query_result.gct")
drug_m <- as.data.frame(drug@mat)
drug_d <- drug@rdesc
dd <- cbind(drug_d, drug_m[rownames(drug_d),])
colnames(dd)[ncol(dd)] <- "norm_cs"

# Filter with same criteria
drug_ds <- dd %>%
    filter(abs(norm_cs) < 2,
        pert_type == "trt_cp",
        fdr_q_nlog10 > 2,
        qc_pass == 1,
        moa != -666)
drug_ds$dose <- as.numeric(sub(' uM', '', drug_ds$pert_idose))

# Generate PM-specific plots
barplot(table(drug_ds$pert_itime), las = 2, col = distcolor[8])
dev.print(pdf, file = "PM_select_all_time_bar.pdf")

hist(drug_ds$ss_ngene, breaks = 100, col = distcolor[9],
    border = "white", xlab = "Signature Strength",
    main = "PM Signature Strength Distribution")
dev.print(pdf, file = "PM_Signature_Strength_all_hist.pdf")

drug_ds %>%
    select(pert_iname, moa) %>%
    distinct() %>%
    group_by(moa) %>%
    summarise(num = n()) %>%
    arrange(desc(num)) %>%
    head(n = 50) %>%
    mutate(moa = ft_sort(moa, num)) %>%
    ggplot(aes(moa, num)) +
    geom_col(fill = "darkred") +
    coord_flip() +
    theme_classic() +
    labs(title = "PM: Top 50 MOA Categories",
        x = "MOA", y = "Number of Drugs")
dev.print(pdf, file = "PM_moa_all_number.pdf")

pmd <- drug_ds

# --------------------------
# 5. COMMON ADJUSTED ANALYSIS
# --------------------------

drug <- parse_gctx("CommonAdj/arfs/TAG/query_result.gct")
drug_m <- as.data.frame(drug@mat)
drug_d <- drug@rdesc
dd <- cbind(drug_d, drug_m[rownames(drug_d),])
colnames(dd)[ncol(dd)] <- "norm_cs"

drug_ds <- dd %>%
    filter(abs(norm_cs) < 2,
        pert_type == "trt_cp",
        fdr_q_nlog10 > 2,
        qc_pass == 1,
        moa != -666)
drug_ds$dose <- as.numeric(sub(' uM', '', drug_ds$pert_idose))

# Generate CommonAdj plots
barplot(table(drug_ds$pert_itime), las = 2, col = distcolor[8])
dev.print(pdf, file = "CommonADJ_select_all_time_bar.pdf")

hist(drug_ds$ss_ngene, breaks = 100, col = distcolor[9],
    border = "white", xlab = "Signature Strength",
    main = "CommonAdj Signature Strength")
dev.print(pdf, file = "CommonADJ_Signature_Strength_all_hist.pdf")

drug_ds %>%
    select(pert_iname, moa) %>%
    distinct() %>%
    group_by(moa) %>%
    summarise(num = n()) %>%
    arrange(desc(num)) %>%
    head(n = 50) %>%
    mutate(moa = ft_sort(moa, num)) %>%
    ggplot(aes(moa, num)) +
    geom_col(fill = "darkgreen") +
    coord_flip() +
    theme_classic() +
    labs(title = "CommonAdj: Top 50 MOA Categories",
        x = "MOA", y = "Number of Drugs")
dev.print(pdf, file = "CommonADJ_moa_all_number.pdf")

adjd <- drug_ds

# --------------------------
# 6. COMMON NON-ADJUSTED ANALYSIS
# --------------------------

drug <- parse_gctx("CommonNoAdj/arfs/TAG/query_result.gct")
drug_m <- as.data.frame(drug@mat)
drug_d <- drug@rdesc
dd <- cbind(drug_d, drug_m[rownames(drug_d),])
colnames(dd)[ncol(dd)] <- "norm_cs"

drug_ds <- dd %>%
    filter(abs(norm_cs) < 2,
        pert_type == "trt_cp",
        fdr_q_nlog10 > 2,
        qc_pass == 1,
        moa != -666)
drug_ds$dose <- as.numeric(sub(' uM', '', drug_ds$pert_idose))

# Generate CommonNoAdj plots
barplot(table(drug_ds$pert_itime), las = 2, col = distcolor[8])
dev.print(pdf, file = "CommonNoADJ_select_all_time_bar.pdf")

hist(drug_ds$ss_ngene, breaks = 100, col = distcolor[9],
    border = "white", xlab = "Signature Strength",
    main = "CommonNoAdj Signature Strength")
dev.print(pdf, file = "CommonNoADJ_Signature_Strength_all_hist.pdf")

drug_ds %>%
    select(pert_iname, moa) %>%
    distinct() %>%
    group_by(moa) %>%
    summarise(num = n()) %>%
    arrange(desc(num)) %>%
    head(n = 50) %>%
    mutate(moa = ft_sort(moa, num)) %>%
    ggplot(aes(moa, num)) +
    geom_col(fill = "purple") +
    coord_flip() +
    theme_classic() +
    labs(title = "CommonNoAdj: Top 50 MOA Categories",
        x = "MOA", y = "Number of Drugs")
dev.print(pdf, file = "CommonNoADJ_moa_all_number.pdf")

adjnd <- drug_ds

# --------------------------
# 7. COMPARATIVE ANALYSES
# --------------------------

# Create subsets for 6h treatment at 10uM dose
tdpds <- subset(tdpd, dose == 10 & pert_itime == "6 h")
pmds <- subset(pmd, dose == 10 & pert_itime == "6 h")
adjds <- subset(adjd, dose == 10 & pert_itime == "6 h")
adjnds <- subset(adjnd, dose == 10 & pert_itime == "6 h")

# Generate Venn diagrams comparing drug sets:

# 1. Selected conditions (6h, 10uM)
venas <- venndetail(list(TDP = tdpds$pert_iname,
    PM = pmds$pert_iname,
    CommonAdj = adjds$pert_iname,
    CommonNoAdj = adjnds$pert_iname))
plot(venas, main = "Drug Overlap (6h, 10uM)")
dev.print(pdf, file = "Drug_venn_select.pdf")

# 2. All conditions
plot(venndetail(list(TDP = tdpd$pert_iname,
    PM = pmd$pert_iname,
    CommonAdj = adjd$pert_iname,
    CommonNoAdj = adjnd$pert_iname)),
    main = "Drug Overlap (All Conditions)")
dev.print(pdf, file = "Drug_venn.pdf")

# Save detailed Venn results
vena <- venndetail(list(TDP = tdpd$pert_iname,
    PM = pmd$pert_iname,
    CommonAdj = adjd$pert_iname,
    CommonNoAdj = adjnd$pert_iname))
write.csv(getFeature(vena,
    rlist = list(TDP = tdpd, PM = pmd,
        CommonAdj = adjd, CommonNoAdj = adjnd),
    userowname = F,
    gind = c(2, 2, 2, 2)),
    file = "All_Drug_venndetail.csv")

# --------------------------
# 8. NETWORK ANALYSIS FUNCTIONS
# --------------------------

#' Create drug-target interaction networks
#'
#' @param dss Dataframe containing drug information (with target_name column)
#' @param rna Dataframe containing RNA expression data
#' @return List containing two igraph objects:
#'         - g1: Full drug-target network
#'         - g2: DEG-filtered drug-target network
get_network <- function(dss, rna) {
    # Split target names and create edge list
    target <- lapply(strsplit(dss$target_name, "\\|"), function(x) gsub(' ', '', x))
    dx <- cbind(rep(dss$pert_iname, times = lengths(target)),
        unlist(target))
    dx <- as.data.frame(dx)
    
    # Remove missing targets
    dx <- subset(dx, V2 != -666)
    
    # Create graph from edge list
    g <- graph_from_data_frame(dx)
    g <- simplify(g)  # Remove duplicate edges
    
    # Prepare RNA data
    rna$gene_symbol <- rna$gene_name
    rna$logFC <- rna$log2FoldChange
    
    # Identify significant DEGs
    deg <- subset(rna, padj < 0.01 & abs(logFC) > 0.1)
    gened <- subset(deg, gene_expresion == "downregulated")$gene_symbol
    geneu <- subset(deg, gene_expresion == "upregulated")$gene_symbol
    
    # Color nodes:
    # - Drugs: lightblue
    # - Upregulated targets: pink
    # - Downregulated targets: cyan4
    # - Other genes: grey
    V(g)$color <- ifelse(V(g)$name %in% dx$V1, "lightblue", "grey")
    V(g)$color <- ifelse(V(g)$name %in% geneu, "pink", V(g)$color)
    V(g)$color <- ifelse(V(g)$name %in% gened, "cyan4", V(g)$color)
    
    # Create DEG-filtered network
    dxx <- subset(dx, V2 %in% deg$gene_symbol)
    g2 <- graph_from_data_frame(dxx)
    g2 <- simplify(g2)
    
    # Apply same coloring scheme
    V(g2)$color <- ifelse(V(g2)$name %in% dxx$V1, "lightblue", "grey")
    V(g2)$color <- ifelse(V(g2)$name %in% geneu, "pink", V(g2)$color)
    V(g2)$color <- ifelse(V(g2)$name %in% gened, "cyan4", V(g2)$color)
    
    return(list(g1 = g, g2 = g2))
}

# --------------------------
# 9. NETWORK VISUALIZATION
# --------------------------

# TDP-43 Network Analysis
gg <- get_network(tdpd, TDP43_deg_sig_strong)

# Plot full network
ggnet2(gg$g1,
    label = V(gg$g1)$name,
    label.size = 3,
    node.color = V(gg$g1)$color,
    node.alpha = 0.8,
    node.size = 5,
    edge.alpha = 0.2) +
    ggtitle("TDP-43: Full Drug-Target Network")
dev.print(pdf, file = "TDP_drug_target_net_color_new.pdf")

# Plot DEG-filtered network
ggnet2(gg$g2,
    label = V(gg$g2)$name,
    label.size = 3,
    node.color = V(gg$g2)$color,
    node.alpha = 0.8,
    node.size = 5,
    edge.alpha = 0.3) +
    ggtitle("TDP-43: DEG-Filtered Drug-Target Network")
dev.print(pdf, file = "TDP_drug_target_net_color_new_deg.pdf")

# PM Network Analysis (similar process)
PM_DEG_filter_strong$padj <- PM_DEG_filter_strong$FDR
PM_DEG_filter_strong$gene_expresion <- ifelse(PM_DEG_filter_strong$log2FoldChange > 0,
    "upregulated", "downregulated")

gg <- get_network(pmd, PM_DEG_filter_strong)

# [Additional network visualizations for PM, CommonAdj, CommonNoAdj...]

# --------------------------
# 10. SPECIAL CASE ANALYSES
# --------------------------

# Analyze specific drug subsets from Venn results
d9 <- venad %>% filter(Subset == "Shared") %>%
    select(Detail, contains("norm_cs"), contains("target_name"))

# [Additional special case network analyses...]

# --------------------------
# 11. SAVE KEY RESULTS
# --------------------------

# Save data for known drugs of interest
write.csv(adjd %>% filter(pert_iname %in% c("tofacitinib", "tamoxifen")),
    file = "ADJ_known_drug.csv")

write.csv(adjnd %>% filter(pert_iname %in% c("tofacitinib", "tamoxifen")),
    file = "ADJ_NO_known_drug.csv")

write.csv(pmd %>% filter(pert_iname %in% c("tofacitinib", "tamoxifen")),
    file = "PM_known_drug.csv")

write.csv(tdpd %>% filter(pert_iname %in% c("tofacitinib", "tamoxifen")),
    file = "TDP_known_drug.csv")

# --------------------------
# 12. SESSION INFORMATION
# --------------------------

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")