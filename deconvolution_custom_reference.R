# TCGA Custom-Reference Deconvolution
#
# Prerequisite: run scRNA-seq-analysis/scripts/export_signature_matrix.py
# and copy its output to data/cell_type_signature_matrix.csv here.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(nnls)
  library(pheatmap)
  library(tibble)
})

project        <- "TCGA-UCEC"
results_dir    <- "results"
data_dir       <- "data"
signature_path <- file.path(data_dir, "cell_type_signature_matrix.csv")

dir.create(results_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

if (!file.exists(signature_path)) {
  stop(
    "Signature matrix not found at ", signature_path, ". Run ",
    "scRNA-seq-analysis/scripts/export_signature_matrix.py and copy its ",
    "output (results/cell_type_signature_matrix.csv there) to ", signature_path, " here."
  )
}

message("Loading custom single-cell-derived signature matrix...")
signature_matrix <- read.csv(signature_path, row.names = 1, check.names = FALSE)
message(
  "Signature matrix: ", nrow(signature_matrix), " genes x ",
  ncol(signature_matrix), " cell types."
)
message("Cell types: ", paste(colnames(signature_matrix), collapse = ", "))

# 1. Download bulk RNA-seq data 
message("Querying GDC for ", project, " RNA-seq data...")

query <- GDCquery(
  project       = project,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, directory = data_dir)
data <- GDCprepare(query, directory = data_dir)

# 2. Use TPM, not raw/DESeq2-normalized counts 
tpm_expr <- assay(data, "tpm_unstrand")
rownames(tpm_expr) <- rowData(data)$gene_name
tpm_expr <- tpm_expr[!is.na(rownames(tpm_expr)) & !duplicated(rownames(tpm_expr)), ]

# 3. Restrict both matrices to their shared gene set
common_genes <- intersect(rownames(signature_matrix), rownames(tpm_expr))
message("Using ", length(common_genes), " genes common to both the signature matrix and bulk data.")

if (length(common_genes) < 50) {
  stop(
    "Fewer than 50 genes overlap between the signature matrix and bulk data."
  )
}

sig_mat  <- as.matrix(signature_matrix[common_genes, , drop = FALSE])
bulk_mat <- as.matrix(tpm_expr[common_genes, , drop = FALSE])

# 4. NNLS deconvolution, one bulk sample at a time
message("Running NNLS deconvolution per sample (this may take a minute for ", ncol(bulk_mat), " samples)...")

deconv_list <- apply(bulk_mat, 2, function(sample_expr) {
  fit <- nnls::nnls(sig_mat, sample_expr)
  proportions <- fit$x / sum(fit$x)  # normalize to sum to 1 per sample
  names(proportions) <- colnames(sig_mat)
  proportions
})
deconv_df <- as.data.frame(t(deconv_list))
deconv_df <- tibble::rownames_to_column(deconv_df, "sample")

write.csv(deconv_df, file.path(results_dir, "custom_reference_deconvolution.csv"), row.names = FALSE)

# 5. Heatmap of estimated cell-type proportions across samples
deconv_mat <- as.matrix(deconv_df[, -1, drop = FALSE])
rownames(deconv_mat) <- deconv_df$sample

png(file.path(results_dir, "custom_reference_deconvolution_heatmap.png"), width = 1400, height = 1000, res = 150)
pheatmap(
  t(deconv_mat),
  show_colnames = FALSE,
  main = paste0(project, ": Cell-Type Proportions (Custom scRNA-seq Reference, NNLS)")
)
dev.off()

# 6. Sanity check: compare against the simpler marker-gene immune
# signature scores from analysis.R, if available 
immune_scores_path <- file.path(results_dir, "immune_signature_scores.csv")
if (file.exists(immune_scores_path) && "T cell" %in% colnames(deconv_mat)) {
  immune_scores <- read.csv(immune_scores_path)
  immune_scores$patient_barcode <- substr(immune_scores$sample, 1, 12)
  deconv_df$patient_barcode <- substr(deconv_df$sample, 1, 12)

  merged_check <- inner_join(
    deconv_df %>% select(patient_barcode, nnls_t_cell = `T cell`),
    immune_scores %>% select(patient_barcode, marker_t_cell = `T.cell`) %>% distinct(patient_barcode, .keep_all = TRUE),
    by = "patient_barcode"
  )

  if (nrow(merged_check) >= 3) {
    cor_check <- cor.test(merged_check$nnls_t_cell, merged_check$marker_t_cell, method = "spearman", exact = FALSE)
    message(
      "\nCross-check: NNLS T-cell proportion vs. marker-gene T-cell signature score, ",
      "Spearman rho = ", round(cor_check$estimate, 3), ", p = ", format.pval(cor_check$p.value, digits = 3)
    )
  }
} else {
  message("\nNo immune_signature_scores.csv found (or no 'T cell' column in the deconvolution) -- skipping cross-check.")
}

message("\nDone. Custom-reference deconvolution results written to '", results_dir, "/'.")
