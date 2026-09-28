# TCGA Cancer Analysis: Differential Expression, Pathway Enrichment,
# Immune Signature Scoring, and Survival Analysis
#
# Downloads RNA-seq and clinical data for one TCGA cohort via TCGAbiolinks,
# runs DESeq2 differential expression, follows up with fgsea pathway
# enrichment (Hallmark gene sets) and marker-gene-based immune signature
# scoring, then runs survival analysis. Designed to be reproducible on a
# laptop for one cohort at a time.
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
project       <- "TCGA-UCEC"   # Uterine Corpus Endometrial Carcinoma; change to any valid TCGA project, e.g. "TCGA-BRCA"
results_dir   <- "results"
data_dir      <- "data"
top_n_genes   <- 20
alpha         <- 0.05
lfc_threshold <- log2(1.5)     # 1.5-fold change cutoff, used alongside padj for "significant" calls
seed          <- 42

set.seed(seed)
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

# 2. Keep one primary tumor and at most one solid normal per patient
# TCGA has multiple aliquots/portions for some patients; keeping all of them
# double-counts patients. Recurrent/metastatic samples are excluded.
n_before <- ncol(data)
data <- data[, select_one_sample_per_patient(colnames(data), c("01", "11"))]
message("Kept ", ncol(data), " of ", n_before, " samples (one primary tumor / solid normal per patient).")

# 3. Define comparison groups
sample_types <- sample_type_code(colnames(data))

colData(data)$group <- determine_group(sample_types, min_normal = 5)

data <- data[, !is.na(colData(data)$group)]
if (ncol(data) == 0) {
  stop(
    "determine_group() returned no usable samples -- likely too few normal ",
    "samples for a tumor-vs-normal comparison (see the warning above). This ",
    "no longer silently falls back to a different analysis; if you want the ",
    "MKI67-based proliferation split, call determine_proliferation_group() ",
    "explicitly as its own, separately-labeled analysis instead."
  )
}
ref_level  <- levels(colData(data)$group)[1]
test_level <- levels(colData(data)$group)[2]
comparison_label <- paste0(test_level, " vs. ", ref_level)
message(
  "Using tumor vs. normal comparison (", sum(data$group == "normal"), " normal, ",
  sum(data$group == "tumor"), " tumor samples)."
)

# 4. Differential expression with DESeq2
message("Running DESeq2...")

keep <- rowSums(assay(data, "unstranded") >= 10) >= 10
data <- data[keep, ]

assay(data, "counts") <- assay(data, "unstranded")
assays(data) <- assays(data)[c("counts", setdiff(names(assays(data)), "counts"))]

dds <- DESeqDataSet(data, design = ~group)
dds <- DESeq(dds)
# Explicit contrast: positive log2FC = higher in test_level (tumor, or MKI67-high).
# Significance tests H0: |log2FC| <= log2(1.5) directly
contrast_vec <- c("group", test_level, ref_level)
res <- results(dds, contrast = contrast_vec, alpha = alpha,
               lfcThreshold = lfc_threshold, altHypothesis = "greaterAbs")

# Standard Wald test (H0: LFC = 0): its statistic is used to rank genes for GSEA,
# which needs a signed score for every gene rather than a thresholded call.
res_wald <- results(dds, contrast = contrast_vec, alpha = alpha)

# Shrunken fold changes (apeglm) for display: stabilizes noisy LFCs of
# low-count genes in the volcano and HSPA plots. Falls back to "normal".
coef_name <- paste0("group_", test_level, "_vs_", ref_level)
stopifnot(coef_name %in% resultsNames(dds))
res_shrunk <- tryCatch(
  lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
  error = function(e) {
    message("apeglm shrinkage unavailable (", conditionMessage(e), "); using type = 'normal'.")
    lfcShrink(dds, coef = coef_name, type = "normal", quiet = TRUE)
  }
)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  mutate(
    log2FoldChange_shrunk = res_shrunk$log2FoldChange[match(gene_id, rownames(res_shrunk))],
    wald_stat = res_wald$stat[match(gene_id, rownames(res_wald))],
    wald_padj = res_wald$padj[match(gene_id, rownames(res_wald))]
  ) %>%
  left_join(
    data.frame(gene_id = rownames(dds), gene_name = rowData(dds)$gene_name, stringsAsFactors = FALSE),
    by = "gene_id"
  ) %>%
  arrange(padj)

res_df$change <- dplyr::case_when(
  !is.na(res_df$padj) & res_df$padj < alpha & res_df$log2FoldChange > 0 ~ "Up",
  !is.na(res_df$padj) & res_df$padj < alpha & res_df$log2FoldChange < 0 ~ "Down",
  TRUE ~ "Not significant"
)
res_df$change <- factor(res_df$change, levels = c("Down", "Not significant", "Up"))

write.csv(res_df, file.path(results_dir, "differential_expression_results.csv"), row.names = FALSE)

message(
  "DE genes with |log2FC| significantly > ", round(lfc_threshold, 3), " (padj < ", alpha, "): ",
  sum(res_df$change != "Not significant"), " (",
  sum(res_df$change == "Up"), " up, ", sum(res_df$change == "Down"), " down). ",
  "For comparison, the LFC = 0 test calls ", sum(!is.na(res_df$wald_padj) & res_df$wald_padj < alpha), "."
)

# 4b. Matched-pair sensitivity analysis (~patient + group)
# this fits a SEPARATE, secondary model on just the patients with both a
# tumor and a normal sample, using design = ~patient + group which
# properly blocks on patient identity. If the primary analysis's top genes
# are not artifacts of ignoring pairing, they should look similar here too.
if (ref_level == "normal") {
  message("Running matched-pair sensitivity analysis...")

  patient_ids <- patient_barcode(colnames(data))
  has_tumor  <- patient_ids[data$group == "tumor"]
  has_normal <- patient_ids[data$group == "normal"]
  matched_patients <- intersect(has_tumor, has_normal)

  MIN_MATCHED_PATIENTS <- 10  # below this, a paired model is underpowered to be informative
  if (length(matched_patients) >= MIN_MATCHED_PATIENTS) {
    paired_idx <- which(patient_ids %in% matched_patients)
    paired_data <- data[, paired_idx]
    colData(paired_data)$patient <- factor(patient_barcode(colnames(paired_data)))

    dds_paired <- DESeqDataSet(paired_data, design = ~patient + group)
    dds_paired <- DESeq(dds_paired)
    res_paired <- results(dds_paired, contrast = contrast_vec, alpha = alpha,
                          lfcThreshold = lfc_threshold, altHypothesis = "greaterAbs")

    paired_df <- as.data.frame(res_paired) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(
        data.frame(gene_id = rownames(dds_paired), gene_name = rowData(dds_paired)$gene_name, stringsAsFactors = FALSE),
        by = "gene_id"
      )
    write.csv(paired_df, file.path(results_dir, "paired_sensitivity_results.csv"), row.names = FALSE)

    # Compare against the primary (unpaired) result on the genes both models tested
    comparison_df <- res_df %>%
      select(gene_id, gene_name, primary_lfc = log2FoldChange, primary_padj = padj) %>%
      inner_join(
        paired_df %>% select(gene_id, paired_lfc = log2FoldChange, paired_padj = padj),
        by = "gene_id"
      ) %>%
      filter(!is.na(primary_lfc) & !is.na(paired_lfc))

    lfc_cor <- stats::cor(comparison_df$primary_lfc, comparison_df$paired_lfc, method = "spearman")
    primary_sig_genes <- comparison_df %>% filter(!is.na(primary_padj) & primary_padj < alpha)
    pct_still_sig <- if (nrow(primary_sig_genes) > 0) {
      100 * mean(!is.na(primary_sig_genes$paired_padj) & primary_sig_genes$paired_padj < alpha)
    } else NA_real_

    write.csv(comparison_df, file.path(results_dir, "paired_sensitivity_comparison.csv"), row.names = FALSE)
    message(
      "Matched-pair sensitivity: ", length(matched_patients), " patients had both a tumor and normal sample. ",
      "log2FC Spearman correlation (primary vs. paired) = ", round(lfc_cor, 3), "; ",
      round(pct_still_sig, 1), "% of the primary analysis's significant genes remain significant ",
      "(padj < ", alpha, ") in the paired model."
    )
  } else {
    message(
      "Only ", length(matched_patients), " patient(s) have both a tumor and normal sample ",
      "(need >= ", MIN_MATCHED_PATIENTS, "); skipping the matched-pair sensitivity analysis."
    )
  }
}

# Variance-stabilized expression, used for every downstream per-sample analysis
# (immune scoring, survival, subtype plots) instead of raw counts.
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)                      # rownames: gene_id
vst_sym <- vst_mat
rownames(vst_sym) <- rowData(dds)$gene_name
vst_sym <- vst_sym[!is.na(rownames(vst_sym)), ]
# When multiple Ensembl IDs map to the same gene symbol, keep the one with
# the HIGHEST mean expression, not just whichever row happened to come
# first. 
if (any(duplicated(rownames(vst_sym)))) {
  gene_mean_expr <- rowMeans(vst_sym)
  keep_idx <- tapply(seq_len(nrow(vst_sym)), rownames(vst_sym), function(idx) {
    idx[which.max(gene_mean_expr[idx])]
  })
  n_collapsed <- sum(duplicated(rownames(vst_sym)))
  vst_sym <- vst_sym[sort(unlist(keep_idx)), ]
  message(
    n_collapsed, " duplicate gene symbol row(s) collapsed to the highest-",
    "mean-expression representative per symbol."
  )
}

tumor_idx <- which(sample_type_code(colnames(dds)) == "01")

# 5. Volcano plot
res_df$plot_p <- pmax(res_df$pvalue, .Machine$double.xmin)

# lfcShrink works from the original counts 
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(plot_p), color = change)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Down" = "steelblue", "Not significant" = "grey70", "Up" = "firebrick")) +
  geom_text_repel(
    data = res_df %>% filter(change != "Not significant") %>% arrange(padj) %>% head(top_n_genes),
    aes(label = gene_name),
    size = 3, color = "black", max.overlaps = 20,
    segment.size = 0.3, segment.color = "grey50"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(project, ": Differential Expression (", comparison_label, ")"),
    subtitle = paste0("Positive log2FC = higher in ", test_level,
                      "\nColored = |log2FC| significantly > log2(1.5) (dashed lines), padj < ", alpha),
    x = "log2 Fold Change", y = "-log10(p-value)", color = NULL
  )

ggsave(file.path(results_dir, "volcano_plot.png"), volcano, width = 7.5, height = 5.5, dpi = 150)
res_df$plot_p <- NULL

# 6. HSPA gene family (Hsp70) expression check
message("Checking HSPA (Hsp70) gene family expression...")

HSPA_GENES <- c(
  "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5",
  "HSPA6", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14"
)

hspa_results <- res_df %>%
  filter(gene_name %in% HSPA_GENES) %>%
  arrange(padj)

# padj comes from the threshold test, so padj < alpha already means |log2FC| > log2(1.5)
hspa_results$is_sig <- !is.na(hspa_results$padj) & hspa_results$padj < alpha

write.csv(hspa_results, file.path(results_dir, "hspa_family_de_results.csv"), row.names = FALSE)

message(
  "Found ", nrow(hspa_results), " HSPA family genes in this dataset; ",
  sum(hspa_results$is_sig), " significantly exceed |log2FC| > ", round(lfc_threshold, 3),
  " (padj < ", alpha, ")."
)

if (nrow(hspa_results) > 0) {
  hspa_results$direction <- ifelse(hspa_results$log2FoldChange > 0, "Up", "Down")

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
      subtitle = paste0(comparison_label, ". Faded = |log2FC| not significantly > log2(1.5) (padj \u2265 ", alpha, ")"),
      x = NULL, y = "log2 Fold Change"
    )
  ggsave(file.path(results_dir, "hspa_family_plot.png"), hspa_plot, width = 7, height = 5, dpi = 150)
}

# 7. Pathway enrichment with fgsea (Hallmark gene sets)
message("Running fgsea pathway enrichment...")
pathway_list <- get_hallmark_sets(local_path = "genesets/h.all.v2023.2.Hs.symbols.gmt")

# Hallmark sets use gene symbols, so ranks must be keyed by symbol too.
# Rank on the LFC = 0 Wald statistic (a signed score for every gene).
ranks <- ranks_by_symbol(res_df, stat_col = "wald_stat")
n_overlap <- length(intersect(names(ranks), unique(unlist(pathway_list))))
message("Ranked ", length(ranks), " genes; ", n_overlap, " overlap the Hallmark gene sets.")
if (n_overlap == 0) stop("No overlap between ranked genes and gene sets; check gene identifiers.")

set.seed(seed)  # fgsea p-values are permutation-based
# eps = 0: estimate very small p-values instead of flooring them at 1e-50
fgsea_res <- fgsea(pathways = pathway_list, stats = ranks, minSize = 15, maxSize = 500, eps = 0) %>%
  arrange(padj)

fgsea_out <- fgsea_res %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(fgsea_out, file.path(results_dir, "gsea_hallmark_results.csv"), row.names = FALSE)

top_pathways <- fgsea_res %>% filter(padj < alpha) %>% arrange(padj) %>% head(15)

if (nrow(top_pathways) > 0) {
  gsea_plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = "steelblue", `TRUE` = "firebrick"), guide = "none") +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0(project, ": Top Enriched Hallmark Pathways"),
      subtitle = paste0(comparison_label, "; positive NES = enriched in ", test_level),
      x = NULL, y = "Normalized Enrichment Score (NES)"
    )
  ggsave(file.path(results_dir, "gsea_hallmark_plot.png"), gsea_plot, width = 9, height = 6, dpi = 150)
  message("Found ", sum(fgsea_res$padj < alpha, na.rm = TRUE), " significantly enriched pathways (padj < ", alpha, ").")
} else {
  message("No significantly enriched pathways found at padj < ", alpha, ".")
}

# 8. Immune and tumor compartment signature scoring (primary tumors only)
# Z-scores are computed across tumors only, so scores are relative to other
# tumors rather than to a tumor/normal mixture.
message("Running immune and tumor compartment signature scoring on ", length(tumor_idx), " primary tumors...")

immune_marker_sets <- list(
  "T cell"               = c("CD3D", "CD3E", "CD3G", "CD2"),
  "CD8 T cell"           = c("CD8A", "CD8B"),
  "B cell"               = c("CD19", "MS4A1", "CD79A"),
  "NK cell"              = c("NCAM1", "KLRD1", "NKG7", "GNLY"),
  "Monocyte/Macrophage"  = c("CD14", "CD68", "CSF1R", "ITGAM"),
  "Dendritic cell"       = c("ITGAX", "CD1C", "CLEC9A"),
  "Neutrophil"           = c("FCGR3B", "CSF3R", "S100A8"),
  # Ayers et al. 2017 (J Clin Invest) 6-gene IFN-gamma signature
  "IFN-gamma"            = c("IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"),
  "Endothelial"          = c("PECAM1", "VWF"),
  "Fibroblast"           = c("COL1A1", "DCN"),
  "Epithelial"           = c("EPCAM", "KRT8", "KRT18", "PAX8")
)

immune_scores <- score_immune_signatures(vst_sym[, tumor_idx, drop = FALSE], immune_marker_sets)

# 8b. Immune-hot / intermediate / immune-cold (IFN-gamma signature tertiles)
hot_cold <- classify_immune_hot_cold(immune_scores[, "IFN-gamma"])

IMMUNE_CELL_TYPES <- c(
  "T cell", "CD8 T cell", "B cell", "NK cell",
  "Monocyte/Macrophage", "Dendritic cell", "Neutrophil", "IFN-gamma"
)

immune_df <- as.data.frame(immune_scores[, IMMUNE_CELL_TYPES, drop = FALSE]) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(hot_cold = hot_cold)
write.csv(immune_df, file.path(results_dir, "immune_signature_scores.csv"), row.names = FALSE)


sample_annotation <- data.frame(immune_class = hot_cold, row.names = rownames(immune_scores))
if (ref_level != "normal") {
  sample_annotation$group <- colData(dds)$group[tumor_idx]
}

png(file.path(results_dir, "immune_signature_heatmap.png"), width = 1400, height = 1000, res = 150)
pheatmap(
  t(immune_scores[, IMMUNE_CELL_TYPES, drop = FALSE]),
  annotation_col = sample_annotation,
  annotation_colors = list(immune_class = HOT_COLD_COLORS),
  show_colnames = FALSE,
  main = paste0(project, ": Immune Signature Scores, primary tumors (mean marker z-score)")
)
invisible(dev.off())

message(project, ": ", paste(names(table(hot_cold)), table(hot_cold), sep = " = ", collapse = ", "),
        " (IFN-gamma signature tertiles).")

hot_cold_plot <- ggplot(immune_df, aes(x = `IFN-gamma`, y = `CD8 T cell`, color = hot_cold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = HOT_COLD_COLORS) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(project, ": Immune-hot vs. Immune-cold Tumors"),
    subtitle = "Classes = tertiles of the Ayers IFN-gamma signature",
    x = "IFN-gamma signature score", y = "CD8 T cell signature score", color = NULL
  )
ggsave(file.path(results_dir, "immune_hot_cold_plot.png"), hot_cold_plot, width = 7, height = 6, dpi = 150)

# 9. Survival analysis (primary tumors, one per patient)
message("Running survival analysis...")

clinical_tumor <- as.data.frame(colData(dds))[tumor_idx, , drop = FALSE]
km_max_days <- 3650  # KM display cap (days); tests still use all follow-up
stage <- pick_stage(clinical_tumor)
stage_col <- NULL
if (!is.null(stage)) {
  clinical_tumor$stage_collapsed <- stage
  stage_col <- "stage_collapsed"
  message("Adjusting Cox models for stage (column '", attr(stage, "source_column"), "').")
} else {
  message("No sufficiently complete stage column found; Cox models are unadjusted.")
}

run_gene_survival <- function(gene_id, gene_label, file_prefix) {
  if (!(gene_id %in% rownames(vst_mat))) return(NULL)
  surv_df <- build_survival_data(clinical_tumor, vst_mat[gene_id, tumor_idx])
  if (nrow(surv_df) < 20 || sum(surv_df$event) < 5 || nlevels(droplevels(surv_df$expr_group)) < 2) {
    message("  Skipping ", gene_label, ": too few patients/events, or median split produced one group.")
    return(NULL)
  }

  stats <- fit_survival_models(surv_df, stage_col = stage_col)

  surv_fit <- survfit(Surv(time, event) ~ expr_group, data = surv_df)
  # Cap at 10 years: few patients remain at risk beyond that, so the tail is noise
  km_plot <- ggsurvplot(
    surv_fit, data = surv_df, pval = TRUE, risk.table = TRUE,
    xlim = c(0, km_max_days), break.time.by = 1000,
    title = paste0(project, ": Survival by ", gene_label, " Expression (primary tumors)"),
    xlab = "Days", legend.title = gene_label, legend.labs = c("Low", "High")
  )
  # ggsave(km_plot$plot) would drop the risk table, so print the full object
  safe_name <- gsub("[^A-Za-z0-9]", "_", gene_label)
  png(file.path(results_dir, paste0(file_prefix, safe_name, ".png")), width = 1050, height = 1050, res = 150)
  print(km_plot, newpage = FALSE)
  invisible(dev.off())

  cbind(data.frame(gene_name = gene_label, gene_id = gene_id), stats)
}

# 9a. Exploratory: top significant DE gene. Differential expression vs. the
# reference group does not imply prognostic value; treat as hypothesis-generating.
top_row <- res_df %>% filter(!is.na(padj), padj < alpha, !is.na(gene_name), nzchar(gene_name)) %>% head(1)
if (nrow(top_row) == 1) {
  top_surv <- run_gene_survival(top_row$gene_id, top_row$gene_name, "survival_km_")
  if (!is.null(top_surv)) {
    write.csv(top_surv, file.path(results_dir, "survival_top_gene_results.csv"), row.names = FALSE)
    file.copy(file.path(results_dir, paste0("survival_km_", gsub("[^A-Za-z0-9]", "_", top_row$gene_name), ".png")),
              file.path(results_dir, "survival_km_plot.png"), overwrite = TRUE)
    message("Survival analysis done for top DE gene ", top_row$gene_name,
            " (log-rank p = ", signif(top_surv$logrank_p, 3), ", Cox p = ", signif(top_surv$cox_p, 3), ").")
  }
} else {
  message("No significant DE gene with a symbol; skipping top-gene survival analysis.")
}

# 9b. Significant HSPA genes, with BH correction across the genes tested
significant_hspa <- hspa_results %>% filter(is_sig)

if (nrow(significant_hspa) > 0) {
  hspa_survival_df <- bind_rows(lapply(seq_len(nrow(significant_hspa)), function(i) {
    run_gene_survival(significant_hspa$gene_id[i], significant_hspa$gene_name[i], "hspa_survival_")
  }))

  if (nrow(hspa_survival_df) > 0) {
    hspa_survival_df$logrank_padj <- p.adjust(hspa_survival_df$logrank_p, method = "BH")
    hspa_survival_df$cox_padj     <- p.adjust(hspa_survival_df$cox_p, method = "BH")
    hspa_survival_df$adj_cox_padj <- p.adjust(hspa_survival_df$adj_cox_p, method = "BH")
    write.csv(hspa_survival_df, file.path(results_dir, "hspa_survival_results.csv"), row.names = FALSE)
    message(
      nrow(hspa_survival_df), " significant HSPA gene(s) tested for survival association; ",
      sum(hspa_survival_df$cox_padj < alpha, na.rm = TRUE),
      " significant after BH correction (Cox, continuous expression)."
    )
  } else {
    message("No HSPA genes had sufficient clinical data for survival analysis.")
  }
} else {
  message("No significant HSPA genes found; skipping HSPA-specific survival analysis.")
}

# 10. HSPA expression by TCGA molecular subtype (primary tumors only)
message("Retrieving TCGA molecular subtype information...")
subtype_lookup <- load_subtype_lookup(sub("^TCGA-", "", project))

if (!is.null(subtype_lookup) && nrow(significant_hspa) > 0) {
  # Subtypes describe tumors: join primary tumor samples only, never normals
  tumor_barcodes <- colnames(dds)[tumor_idx]
  subtype_by_sample <- data.frame(
    sample = tumor_barcodes,
    patient_barcode = patient_barcode(tumor_barcodes),
    stringsAsFactors = FALSE
  ) %>%
    inner_join(subtype_lookup, by = "patient_barcode")

  # Subtype labels come from the 2013 marker paper and cover only part of the
  # cohort; "Not assigned" is a mixed group, not a subtype.
  subtype_coverage <- data.frame(
    tumors_analyzed = length(tumor_barcodes),
    matched_to_subtype_table = nrow(subtype_by_sample),
    not_assigned = sum(subtype_by_sample$molecular_subtype == "Not assigned"),
    assigned_subtype = sum(subtype_by_sample$molecular_subtype != "Not assigned")
  )
  write.csv(subtype_coverage, file.path(results_dir, "subtype_coverage.csv"), row.names = FALSE)
  message(
    "Subtype coverage: ", subtype_coverage$assigned_subtype, " of ", subtype_coverage$tumors_analyzed,
    " primary tumors have an assigned subtype (", subtype_coverage$not_assigned, " 'Not assigned', ",
    subtype_coverage$tumors_analyzed - subtype_coverage$matched_to_subtype_table, " absent from the subtype table)."
  )

  if (nrow(subtype_by_sample) > 0) {
    hspa_subtype_df <- bind_rows(lapply(seq_len(nrow(significant_hspa)), function(i) {
      gid <- significant_hspa$gene_id[i]
      if (!(gid %in% rownames(vst_mat))) return(NULL)
      data.frame(
        gene_name = significant_hspa$gene_name[i],
        sample = subtype_by_sample$sample,
        molecular_subtype = subtype_by_sample$molecular_subtype,
        vst_expression = as.numeric(vst_mat[gid, subtype_by_sample$sample])
      )
    }))
    write.csv(hspa_subtype_df, file.path(results_dir, "hspa_by_subtype.csv"), row.names = FALSE)

    # Kruskal-Wallis across assigned subtypes, BH-corrected across genes
    kw <- hspa_subtype_df %>%
      group_by(gene_name) %>%
      summarise(kruskal_p = kruskal_p(vst_expression, molecular_subtype), .groups = "drop") %>%
      mutate(kruskal_padj = p.adjust(kruskal_p, method = "BH"))
    write.csv(kw, file.path(results_dir, "hspa_subtype_kruskal.csv"), row.names = FALSE)

    hspa_subtype_df$molecular_subtype <- order_subtypes(hspa_subtype_df$molecular_subtype)

    facet_labels <- setNames(
      paste0(kw$gene_name, "  (KW padj = ", formatC(kw$kruskal_padj, format = "g", digits = 2), ")"),
      kw$gene_name
    )

    subtype_plot <- ggplot(hspa_subtype_df, aes(x = molecular_subtype, y = vst_expression, fill = molecular_subtype)) +
      geom_boxplot(outlier.size = 0.5) +
      facet_wrap(~gene_name, scales = "free_y", labeller = as_labeller(facet_labels)) +
      scale_fill_manual(values = subtype_palette(levels(hspa_subtype_df$molecular_subtype))) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(
        title = paste0(project, ": Significant HSPA Genes by Molecular Subtype (primary tumors)"),
        subtitle = "Kruskal-Wallis across assigned subtypes ('Not assigned', grey, excluded), BH-corrected across genes",
        x = NULL, y = "VST expression"
      )
    ggsave(file.path(results_dir, "hspa_by_subtype_plot.png"), subtype_plot, width = 10, height = 6.5, dpi = 150)

    message(
      "HSPA-by-subtype analysis complete: ", nrow(subtype_by_sample),
      " primary tumors matched to a molecular subtype; ",
      sum(kw$kruskal_padj < alpha, na.rm = TRUE), " gene(s) differ by subtype (BH padj < ", alpha, ")."
    )
  } else {
    message("No samples could be matched to a molecular subtype; skipping HSPA-by-subtype plot.")
  }
} else if (nrow(significant_hspa) == 0) {
  message("No significant HSPA genes found; skipping HSPA-by-subtype analysis.")
}

writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
message("Done. All results written to '", results_dir, "/'.")
