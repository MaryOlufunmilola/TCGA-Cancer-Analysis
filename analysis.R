# TCGA Cancer Analysis: Differential Expression, Pathway Enrichment,
# Immune Signature Scoring, and Survival Analysis
#
# Downloads RNA-seq and clinical data for one TCGA cohort via TCGAbiolinks,
# runs DESeq2 differential expression, follows up with fgsea pathway
# enrichment (Hallmark gene sets) and marker-gene-based immune signature
# scoring, then runs Kaplan-Meier survival analysis. Designed to be
# reproducible on a laptop for one cohort at a time.
#
# Default cohort: TCGA-UCEC (Uterine Corpus Endometrial Carcinoma)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)
  library(msigdbr)
  library(tibble)
  library(tidyr)
  library(pheatmap)
})

source("R/functions.R")

# Config
project     <- "TCGA-UCEC"        # Uterine Corpus Endometrial Carcinoma; change to any valid TCGA project, e.g. "TCGA-BRCA"
results_dir <- "results"
data_dir    <- "data"
top_n_genes <- 20
alpha       <- 0.05
lfc_threshold <- log2(1.5)  # 1.5-fold change cutoff, used alongside padj for "significant" calls

dir.create(results_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

# 1. Query and download RNA-seq counts
message("Querying GDC for ", project, " RNA-seq data...")

query <- GDCquery(
  project       = project,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, directory = data_dir)
data <- GDCprepare(query, directory = data_dir)

# 2. Define comparison groups
sample_types <- substr(colnames(data), 14, 15)
n_normal <- sum(sample_types == "11")
message(
  if (n_normal >= 5) paste0("Using tumor vs. normal comparison (", n_normal, " normal samples available).")
  else paste0("Too few normal samples (", n_normal, "); falling back to expression-based split.")
)
colData(data)$group <- determine_group(sample_types, assay(data, "unstranded"))

# 3. Differential expression with DESeq2
message("Running DESeq2...")

keep <- rowSums(assay(data, "unstranded") >= 10) >= 10
data <- data[keep, ]

assay(data, "counts") <- assay(data, "unstranded")
assays(data) <- assays(data)[c("counts", setdiff(names(assays(data)), "counts"))]

dds <- DESeqDataSet(data, design = ~group)
dds <- DESeq(dds)
res <- results(dds, alpha = alpha)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

# Join gene symbols in early so every downstream plot and CSV has gene_name
gene_id_to_symbol <- data.frame(
  gene_id = rownames(data),
  gene_name = rowData(data)$gene_name,
  stringsAsFactors = FALSE
)
res_df <- res_df %>% left_join(gene_id_to_symbol, by = "gene_id")

res_df$change <- dplyr::case_when(
  !is.na(res_df$padj) & res_df$padj < alpha & res_df$log2FoldChange > lfc_threshold  ~ "Up",
  !is.na(res_df$padj) & res_df$padj < alpha & res_df$log2FoldChange < -lfc_threshold ~ "Down",
  TRUE ~ "Not significant"
)
res_df$change <- factor(res_df$change, levels = c("Down", "Not significant", "Up"))

write.csv(res_df, file.path(results_dir, "differential_expression_results.csv"), row.names = FALSE)

top_genes <- head(res_df$gene_id[!is.na(res_df$padj) & res_df$padj < alpha], top_n_genes)
message("Top DE genes (padj < ", alpha, "): ", length(top_genes), " found.")

# 4. Volcano plot
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Down" = "steelblue", "Not significant" = "grey70", "Up" = "firebrick")) +
  geom_text_repel(
    data = res_df %>% filter(change != "Not significant") %>% arrange(padj) %>% head(20),
    aes(label = gene_name),
    size = 3, color = "black", max.overlaps = 20,
    segment.size = 0.3, segment.color = "grey50"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(project, ": Differential Expression (", levels(colData(data)$group)[1],
                    " vs. ", levels(colData(data)$group)[2], ")"),
    x = "log2 Fold Change", y = "-log10(p-value)", color = NULL
  )

ggsave(file.path(results_dir, "volcano_plot.png"), volcano, width = 7, height = 5, dpi = 150)

# 5. HSPA gene family (Hsp70) expression check
message("Checking HSPA (Hsp70) gene family expression...")

HSPA_GENES <- c(
  "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5",
  "HSPA6", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14"
)

hspa_results <- res_df %>%
  filter(gene_name %in% HSPA_GENES) %>%
  arrange(padj)

write.csv(hspa_results, file.path(results_dir, "hspa_family_de_results.csv"), row.names = FALSE)

n_hspa_found <- nrow(hspa_results)
n_hspa_sig <- sum(
  !is.na(hspa_results$padj) & hspa_results$padj < alpha &
    abs(hspa_results$log2FoldChange) > lfc_threshold,
  na.rm = TRUE
)
message(
  "Found ", n_hspa_found, " HSPA family genes in this dataset; ",
  n_hspa_sig, " reached significance (padj < ", alpha, " and |log2FC| > ", round(lfc_threshold, 3), ")."
)

if (n_hspa_found > 0) {
  hspa_results$direction <- ifelse(hspa_results$log2FoldChange > 0, "Up", "Down")
  hspa_results$is_sig <- !is.na(hspa_results$padj) & hspa_results$padj < alpha &
    abs(hspa_results$log2FoldChange) > lfc_threshold

  hspa_plot <- ggplot(
    hspa_results,
    aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange,
        fill = direction, alpha = is_sig)
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue"), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0(project, ": HSPA (Hsp70) Family Differential Expression"),
      subtitle = "Red = up-regulated, blue = down-regulated; faded = not significant (padj \u2265 0.05)",
      x = NULL, y = "log2 Fold Change"
    )
  ggsave(file.path(results_dir, "hspa_family_plot.png"), hspa_plot, width = 7, height = 5, dpi = 150)
}

# 6. Pathway enrichment with fgsea (Hallmark gene sets)
message("Running fgsea pathway enrichment...")
pathway_list <- get_hallmark_sets(local_path = "genesets/h.all.v2023.2.Hs.symbols.gmt")

# Rank genes by DESeq2 test statistic for preranked GSEA
ranked_genes <- res_df %>%
  filter(!is.na(stat)) %>%
  distinct(gene_id, .keep_all = TRUE)
ranks <- setNames(ranked_genes$stat, ranked_genes$gene_id)
ranks <- sort(ranks, decreasing = TRUE)

fgsea_res <- fgsea(pathways = pathway_list, stats = ranks, minSize = 15, maxSize = 500)
fgsea_res <- fgsea_res %>% arrange(padj)

fgsea_out <- fgsea_res %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(fgsea_out, file.path(results_dir, "gsea_hallmark_results.csv"), row.names = FALSE)

top_pathways <- fgsea_res %>% filter(padj < alpha) %>% arrange(padj) %>% head(15)

if (nrow(top_pathways) > 0) {
  gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("steelblue", "firebrick"), guide = "none") +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0(project, ": Top Enriched Hallmark Pathways"),
      x = NULL, y = "Normalized Enrichment Score (NES)"
    )
  ggsave(file.path(results_dir, "gsea_hallmark_plot.png"), gsea_plot, width = 9, height = 6, dpi = 150)
  message("Found ", nrow(top_pathways), " significantly enriched pathways (padj < ", alpha, ").")
} else {
  message("No significantly enriched pathways found at padj < ", alpha, ".")
}

# 7. Immune and tumor compartment signature scoring
message("Running immune and tumor compartment signature scoring...")

immune_marker_sets <- list(
  "T cell"               = c("CD3D", "CD3E", "CD3G", "CD2"),
  "CD8 T cell"           = c("CD8A", "CD8B"),
  "B cell"               = c("CD19", "MS4A1", "CD79A"),
  "NK cell"              = c("NCAM1", "KLRD1", "NKG7", "GNLY"),
  "Monocyte/Macrophage"  = c("CD14", "CD68", "CSF1R", "ITGAM"),
  "Dendritic cell"       = c("ITGAX", "CD1C", "CLEC9A"),
  "Neutrophil"           = c("FCGR3B", "CSF3R", "S100A8"),
  "Endothelial"          = c("PECAM1", "VWF"),
  "Fibroblast"           = c("COL1A1", "DCN"),
  "Epithelial"           = c("EPCAM", "KRT8", "KRT18", "PAX8")
)

norm_counts <- counts(dds, normalized = TRUE)
rownames(norm_counts) <- rowData(data)$gene_name
norm_counts <- norm_counts[!is.na(rownames(norm_counts)) & !duplicated(rownames(norm_counts)), ]
log_counts <- log2(norm_counts + 1)

immune_scores <- score_immune_signatures(log_counts, immune_marker_sets)

IMMUNE_CELL_TYPES <- c(
  "T cell", "CD8 T cell", "B cell", "NK cell",
  "Monocyte/Macrophage", "Dendritic cell", "Neutrophil"
)
immune_scores_display <- immune_scores[, IMMUNE_CELL_TYPES, drop = FALSE]

write.csv(
  as.data.frame(immune_scores_display) %>% tibble::rownames_to_column("sample"),
  file.path(results_dir, "immune_signature_scores.csv"), row.names = FALSE
)

sample_annotation <- data.frame(group = colData(data)$group)
rownames(sample_annotation) <- rownames(immune_scores_display)

png(file.path(results_dir, "immune_signature_heatmap.png"), width = 1400, height = 1000, res = 150)
pheatmap(
  t(immune_scores_display),
  annotation_col = sample_annotation,
  show_colnames = FALSE,
  main = paste0(project, ": Immune Signature Scores (mean marker-gene z-score)")
)
dev.off()

message("Immune signature scoring complete. Results in immune_signature_scores.csv and immune_signature_heatmap.png")

# 7b. Immune-hot vs. immune-cold tumor classification
# samples with high T cell signature relative to their Epithelial 
# signature are classified "immune-hot", the reverse "immune-cold" 
message("Classifying samples as immune-hot vs. immune-cold...")

immune_hot_cold <- as.data.frame(immune_scores) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(
    hot_cold = ifelse(`T cell` > `Epithelial`, "Immune-hot", "Immune-cold")
  )

write.csv(immune_hot_cold, file.path(results_dir, "immune_hot_cold_classification.csv"), row.names = FALSE)

n_hot <- sum(immune_hot_cold$hot_cold == "Immune-hot", na.rm = TRUE)
n_cold <- sum(immune_hot_cold$hot_cold == "Immune-cold", na.rm = TRUE)
message(project, ": ", n_hot, " samples classified immune-hot, ", n_cold, " immune-cold.")

hot_cold_plot <- ggplot(immune_hot_cold, aes(x = `Epithelial`, y = `T cell`, color = hot_cold)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Immune-hot" = "firebrick", "Immune-cold" = "steelblue")) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(project, ": Immune-hot vs. Immune-cold Classification"),
    subtitle = "T cell signature score vs. Epithelial signature score (dashed line = equal)",
    x = "Epithelial signature score", y = "T cell signature score", color = NULL
  )
ggsave(file.path(results_dir, "immune_hot_cold_plot.png"), hot_cold_plot, width = 7, height = 6, dpi = 150)

# 8. Survival analysis on top gene
message("Running survival analysis...")

clinical <- colData(data) %>% as.data.frame()
top_gene_id <- res_df$gene_id[1]
top_gene_name <- res_df$gene_name[1]
top_gene_label <- if (!is.na(top_gene_name) && nzchar(top_gene_name)) top_gene_name else top_gene_id

if (!is.null(top_gene_id) && top_gene_id %in% rownames(assay(data, "unstranded"))) {
  expr_top <- assay(data, "unstranded")[top_gene_id, ]
  clinical$expr_group <- ifelse(expr_top > median(expr_top), "High", "Low")

  clinical$time  <- clinical$days_to_death
  clinical$time[is.na(clinical$time)] <- clinical$days_to_last_follow_up[is.na(clinical$time)]
  clinical$event <- ifelse(clinical$vital_status == "Dead", 1, 0)

  clinical_surv <- clinical %>% filter(!is.na(time) & !is.na(event))

  surv_fit <- survfit(Surv(time, event) ~ expr_group, data = clinical_surv)

  km_plot <- ggsurvplot(
    surv_fit, data = clinical_surv, pval = TRUE, risk.table = TRUE,
    title = paste0(project, ": Survival by ", top_gene_label, " Expression"),
    xlab = "Days", legend.title = top_gene_label
  )

  ggsave(file.path(results_dir, "survival_km_plot.png"), km_plot$plot, width = 7, height = 5, dpi = 150)
  message("Survival plot saved for gene: ", top_gene_label, " (", top_gene_id, ")")
} else {
  message("No significant top gene found for survival analysis; skipping.")
}

message("Done. All results written to '", results_dir, "/'.")

# 9. HSPA gene family survival analysis
message("Running survival analysis for significant HSPA genes...")

run_km_survival <- function(gene_id, gene_label, expr_matrix, clinical_data, results_dir, project) {
  if (!(gene_id %in% rownames(expr_matrix))) return(NULL)

  expr_gene <- expr_matrix[gene_id, ]
  clin <- clinical_data
  clin$expr_group <- ifelse(expr_gene > median(expr_gene), "High", "Low")
  clin$time  <- clin$days_to_death
  clin$time[is.na(clin$time)] <- clin$days_to_last_follow_up[is.na(clin$time)]
  clin$event <- ifelse(clin$vital_status == "Dead", 1, 0)
  clin_surv <- clin %>% filter(!is.na(time) & !is.na(event))

  if (nrow(clin_surv) < 10 || length(unique(clin_surv$expr_group)) < 2) return(NULL)

  surv_fit <- survfit(Surv(time, event) ~ expr_group, data = clin_surv)
  surv_diff <- survdiff(Surv(time, event) ~ expr_group, data = clin_surv)
  logrank_p <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

  km_plot <- ggsurvplot(
    surv_fit, data = clin_surv, pval = TRUE, risk.table = TRUE,
    title = paste0(project, ": Survival by ", gene_label, " Expression"),
    xlab = "Days", legend.title = gene_label
  )

  safe_name <- gsub("[^A-Za-z0-9]", "_", gene_label)
  ggsave(file.path(results_dir, paste0("hspa_survival_", safe_name, ".png")),
         km_plot$plot, width = 7, height = 5, dpi = 150)

  data.frame(gene_name = gene_label, gene_id = gene_id, logrank_p = logrank_p)
}

significant_hspa <- hspa_results %>% filter(is_sig)

if (nrow(significant_hspa) > 0) {
  clinical_base <- colData(data) %>% as.data.frame()
  expr_matrix <- assay(data, "unstranded")

  hspa_survival_list <- lapply(seq_len(nrow(significant_hspa)), function(i) {
    run_km_survival(
      gene_id = significant_hspa$gene_id[i],
      gene_label = significant_hspa$gene_name[i],
      expr_matrix = expr_matrix,
      clinical_data = clinical_base,
      results_dir = results_dir,
      project = project
    )
  })
  hspa_survival_df <- bind_rows(hspa_survival_list)

  if (nrow(hspa_survival_df) > 0) {
    write.csv(hspa_survival_df, file.path(results_dir, "hspa_survival_results.csv"), row.names = FALSE)
    n_sig_survival <- sum(hspa_survival_df$logrank_p < alpha, na.rm = TRUE)
    message(
      nrow(hspa_survival_df), " significant HSPA gene(s) tested for survival association; ",
      n_sig_survival, " showed a significant survival difference (log-rank p < ", alpha, ")."
    )
  } else {
    message("No HSPA genes had sufficient clinical data for survival analysis.")
  }
} else {
  message("No significant HSPA genes found; skipping HSPA-specific survival analysis.")
}

# 10. HSPA expression by TCGA molecular subtype
message("Retrieving TCGA molecular subtype information...")

tumor_code <- sub("^TCGA-", "", project)

subtype_info <- tryCatch(
  TCGAquery_subtype(tumor = tumor_code),
  error = function(e) {
    message("Could not retrieve molecular subtype info for ", tumor_code, ": ", conditionMessage(e))
    NULL
  }
)

if (!is.null(subtype_info) && nrow(significant_hspa) > 0) {
  message("Subtype table columns: ", paste(colnames(subtype_info), collapse = ", "))

  patient_col <- intersect(c("patient", "Patient", "bcr_patient_barcode"), colnames(subtype_info))
  subtype_col_candidates <- grep("subtype|cluster|integrat", colnames(subtype_info), ignore.case = TRUE, value = TRUE)

  if (length(patient_col) == 0 || length(subtype_col_candidates) == 0) {
    message(
      "Could not auto-detect patient/subtype columns from the table above. ",
      "Inspect colnames(subtype_info) directly and adjust this section to match."
    )
  } else {
    patient_col <- patient_col[1]
    # Prefer a column with "subtype" in the name over a more generic "cluster" match
    subtype_col <- if (any(grepl("subtype", subtype_col_candidates, ignore.case = TRUE))) {
      grep("subtype", subtype_col_candidates, ignore.case = TRUE, value = TRUE)[1]
    } else {
      subtype_col_candidates[1]
    }
    message("Using patient column '", patient_col, "' and subtype column '", subtype_col, "'.")

    subtype_lookup <- subtype_info[, c(patient_col, subtype_col)]
    colnames(subtype_lookup) <- c("patient_barcode", "molecular_subtype")

    expr_matrix <- assay(data, "unstranded")
    sample_barcodes <- colnames(expr_matrix)
    patient_barcodes <- substr(sample_barcodes, 1, 12)

    subtype_by_sample <- data.frame(
      sample = sample_barcodes,
      patient_barcode = patient_barcodes,
      stringsAsFactors = FALSE
    ) %>%
      left_join(subtype_lookup, by = "patient_barcode") %>%
      filter(!is.na(molecular_subtype))

    if (nrow(subtype_by_sample) > 0) {
      hspa_subtype_rows <- lapply(seq_len(nrow(significant_hspa)), function(i) {
        gid <- significant_hspa$gene_id[i]
        gname <- significant_hspa$gene_name[i]
        if (!(gid %in% rownames(expr_matrix))) return(NULL)

        expr_vals <- expr_matrix[gid, subtype_by_sample$sample]
        data.frame(
          gene_name = gname,
          sample = subtype_by_sample$sample,
          molecular_subtype = subtype_by_sample$molecular_subtype,
          expression = as.numeric(expr_vals)
        )
      })
      hspa_subtype_df <- bind_rows(hspa_subtype_rows)

      write.csv(hspa_subtype_df, file.path(results_dir, "hspa_by_subtype.csv"), row.names = FALSE)

      subtype_plot <- ggplot(hspa_subtype_df, aes(x = molecular_subtype, y = log2(expression + 1), fill = molecular_subtype)) +
        geom_boxplot(outlier.size = 0.5) +
        facet_wrap(~gene_name, scales = "free_y") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        labs(
          title = paste0(project, ": Significant HSPA Genes by Molecular Subtype"),
          x = NULL, y = "log2(count + 1)"
        )
      ggsave(file.path(results_dir, "hspa_by_subtype_plot.png"), subtype_plot, width = 9, height = 6, dpi = 150)

      message(
        "HSPA-by-subtype analysis complete: ", nrow(subtype_by_sample), " samples matched to a molecular subtype. ",
        "Results in hspa_by_subtype.csv and hspa_by_subtype_plot.png"
      )
    } else {
      message("No samples could be matched to a molecular subtype; skipping HSPA-by-subtype plot.")
    }
  }
} else if (nrow(significant_hspa) == 0) {
  message("No significant HSPA genes found; skipping HSPA-by-subtype analysis.")
}
