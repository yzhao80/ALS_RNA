# Generate enriched pathway results (Hallmark + KEGG) for DEGs
# Author: Yue Zhao
# Date: 2025-06-15

library(dplyr)
library(purrr)
library(stringr)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
setwd("C:/Users/yuedz/OneDrive - Michigan Medicine/Documents/Project/ALS/RNA-seq/")

# 1. Load DEGs (without and with cell-type adjustment) ---------------------------------
load("UM_deg.rda") # deg_all,deg_male,deg_female,PC8_all_df,PC8_m_df,PC8_f_df

# 2. Helper functions for ENTREZ conversion, hallmark geneset clean up, and KEGG/Hallmark enrichment ---------------------------------------
get_gene_list <- function(df, pval_thresh = 0.01, lfc_thresh = 0.1) {
    # Keep highest AveExpr per gene_symbol
    df_unique <- df %>% group_by(gene_symbol) %>% slice_max(AveExpr, n = 1, with_ties = FALSE) %>% ungroup()
    
    # Map symbols to ENTREZID
    mapping <- bitr(unique(df_unique$gene_symbol),fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    # Join and split up/down
    df_mapped <- df_unique %>% inner_join(mapping, by = c("gene_symbol" = "SYMBOL"))
    
    list(up = df_mapped %>% filter(adj.P.Val < pval_thresh, logFC >  lfc_thresh) %>% pull(ENTREZID) %>% unique(),
        down = df_mapped %>% filter(adj.P.Val < pval_thresh, logFC < -lfc_thresh) %>% pull(ENTREZID) %>% unique())
}
genes_no <- list(Overall = get_gene_list(deg_all),Males = get_gene_list(deg_male), Females = get_gene_list(deg_female))

clean_description <- function(descs) {
    cleaned <- descs %>%
        tolower() %>%
        str_remove("hallmark_") %>%
        str_replace_all("_", " ") %>%
        str_to_sentence()
    cleaned
}

enrich_pathways <- function(ids, type) {
    if (type == "H") {
        term2gene <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
        enr <- enricher(ids, TERM2GENE = term2gene)
    } else {
        enr <- enrichKEGG(gene = ids, organism = 'hsa', pvalueCutoff = 0.99, qvalueCutoff = 0.99)
    }
    setReadable(enr, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')@result %>% as_tibble()
}

# 3. Run enrichment & combine results for no cell adjustment (primary analysis) --------------------------------------
enrich_no <- lapply(names(genes_no), function(cmp) {
    gl <- genes_no[[cmp]]
    
    # Hallmark up/down
    hm_up   <- enrich_pathways(gl$up,   "H") %>% filter(p.adjust < 0.05)
    hm_down <- enrich_pathways(gl$down, "H") %>% filter(p.adjust < 0.05)
    
    hallmark <- bind_rows(hm_up, hm_down) %>%
        dplyr::mutate(
            Description = clean_description(Description),
            source      = "Hallmark",
            comparison  = cmp,
            method   = if_else(row_number() <= nrow(hm_up), "Up", "Down"),
            category    = "Hallmark",
            subcategory = "Hallmark"
        )
    
    # KEGG up/down
    kk_up   <- enrich_pathways(gl$up,   "KEGG") %>% filter(p.adjust < 0.05)
    kk_down <- enrich_pathways(gl$down, "KEGG") %>% filter(p.adjust < 0.05)
    
    kegg <- bind_rows(
        kk_up   %>% mutate(method = "Up"),
        kk_down %>% mutate(method = "Down")
    ) %>%
        mutate(
            source     = "KEGG",
            comparison = cmp
        )
    
    bind_rows(hallmark, kegg)
})

plot_no <- bind_rows(enrich_no) %>% filter(!is.na(Description), Description != "") %>%
    dplyr::rename(FDR = p.adjust) %>%
    mutate(
        comparison = factor(comparison, levels = c("Overall", "Males", "Females")),
        method  = factor(method,  levels = c("Up", "Down")),
        Description = factor(Description, levels = unique(Description))
    )

# 4. Repeat enrichment for cell-type–adjusted DEGs ------------------------
# DEGs with cell-type adjustment already loaded as PC8_all_df, PC8_m_df, PC8_f_df
# Build genes_yes
genes_yes <- list(Overall = get_gene_list(PC8_all_df), Males   = get_gene_list(PC8_m_df), Females = get_gene_list(PC8_f_df))

# Run enrichment for yes
enrich_yes <- lapply(names(genes_yes), function(cmp) {
    gl <- genes_yes[[cmp]]
    
    # Hallmark
    hm_up   <- enrich_pathways(gl$up,   "H") %>% filter(p.adjust < 0.05)
    hm_down <- enrich_pathways(gl$down, "H") %>% filter(p.adjust < 0.05)
    
    hallmark <- bind_rows(hm_up, hm_down) %>%
        mutate(
            Description = clean_description(Description),
            source      = "Hallmark",
            comparison  = cmp,
            method      = rep(c("Up","Down"), c(nrow(hm_up), nrow(hm_down))),
            category    = "Hallmark",
            subcategory = "Hallmark"
        )
    
    # KEGG
    kk_up   <- enrich_pathways(gl$up,   "KEGG") %>% filter(p.adjust < 0.05)
    kk_down <- enrich_pathways(gl$down, "KEGG") %>% filter(p.adjust < 0.05)
    
    kegg <- bind_rows(
        kk_up   %>% mutate(method = "Up"),
        kk_down %>% mutate(method = "Down")
    ) %>%
        mutate(
            source     = "KEGG",
            comparison = cmp
        )
    
    bind_rows(hallmark, kegg)
})

plot_yes <- bind_rows(enrich_yes) %>% filter(!is.na(Description), Description != "") %>%
    dplyr::rename(FDR = p.adjust) %>%
    mutate(
        comparison = factor(comparison, levels = c("Overall", "Males", "Females")),
        method     = factor(method,     levels = c("Up", "Down")),
        Description = factor(Description, levels = unique(Description))
    )

plot_no$method2 <- paste(plot_no$comparison, plot_no$method, sep = "_")
plot_no$comparison2<-"Not adjust celltype"
plot_yes$method2 <- paste(plot_yes$comparison, plot_yes$method, sep = "_")
plot_yes$comparison2<-"Adjust celltype"
plot_df3<-rbind.data.frame(plot_no,plot_yes)

save(plot_no,plot_yes,plot_df3, file="plot_no_yes_df3.rda")
