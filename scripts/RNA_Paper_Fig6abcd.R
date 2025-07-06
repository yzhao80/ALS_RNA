# Overlapping genes between blood DEGs and a model of early ALS, iPSC-neurons with TDP-43 knockdown, and a model of late ALS, postmortem spinal cord tissue
# Author: Yue Zhao
# Date: 2024-11-15

library(dplyr)
library(gplots)

load("../data/Fig6abcd_input.rda") #all_deg_sig contains primary blood DEGs, DEG_filtered contains adjusted blood DEGs.

# Fig 6a: Overlap with iPSC-derived neuron DEG between TDP43 knock down vs. controls
TDP43KDvsCtr_DEG <- read.csv("../data/UNC13A paper_Differentially Expressed Genes.csv")
TDP43_deg_sig_strong <- filter(TDP43KDvsCtr_DEG, padj<0.01 & abs(log2FoldChange)>0.1) #2742
# with cell type adjustment
DEG=DEG_filtered
gene_up= unlist(unique(DEG[DEG$adj.P.Val<0.01&DEG$logFC>0.1,'gene_symbol'] ))
gene_down= unlist(unique(DEG[DEG$adj.P.Val<0.01&DEG$logFC< -0.1,'gene_symbol']))
compare_list8 <- list(our_study_noCellAdjust=unique(all_deg_sig$gene_symbol),our_study_withCellAdjust=unique(c(gene_up,gene_down)),iPCS_neuron=unique(TDP43_deg_sig_strong$gene_name))
venn(compare_list8)

# Fig 6b: Overlap with Post-mortem spinal cord
PM_DEG <- read.csv("../data/Post_mortem_spinal_cord_DEG.csv")
PM_DEG$FDR <- p.adjust(PM_DEG$pvalue, method="fdr")
PM_DEG_filter <- filter(PM_DEG, FDR<0.05)
PM_DEG_filter_strong <- filter(PM_DEG, FDR<0.01 & abs(log2FoldChange)>0.1)
compare_list9 <- list(our_study_noCellAdjust=unique(all_deg_sig$gene_symbol),our_study_withCellAdjust=unique(c(gene_up,gene_down)),spinal_cord=unique(PM_DEG_filter_strong$gene_name))
venn(compare_list9)

# Fig 6c: Overlap DEGs between our study (adjusted cell types), neurons, and spinal cord
commonThree <- list(withCellAdjust=unique(c(gene_up,gene_down)),iPCS_neuron=unique(TDP43_deg_sig_strong$gene_name),spinal_cord=unique(PM_DEG_filter_strong$gene_name))
venn(commonThree)

# Fig 6d: Overlap DEGs between our study (primary analysis), neurons, and spinal cord
commonThreeNoA <- list(noCellAdjust=unique(all_deg_sig$gene_symbol),iPCS_neuron=unique(TDP43_deg_sig_strong$gene_name),spinal_cord=unique(PM_DEG_filter_strong$gene_name))
venn(commonThreeNoA)
