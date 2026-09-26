# TCGA Mutation Analysis: MAF loading, tumor mutational burden (TMB),
# top mutated genes, HSPA (Hsp70) family mutation frequency, and
# correlation with immune signature scores.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(maftools)
  library(dplyr)
  library(ggplot2)
})

source("R/functions.R")

# Config
project     <- "TCGA-UCEC"   # keep in sync with analysis.R
results_dir <- "results"
data_dir    <- "data"

dir.create(results_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

# 1. Download and load MAF data
message("Querying GDC for ", project, " somatic mutation data...")

maf_query <- GDCquery(
  project       = project,
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  access        = "open",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(maf_query, directory = data_dir)
maf_df <- GDCprepare(maf_query, directory = data_dir)
maf <- read.maf(maf = maf_df)

message("Loaded MAF: ", nrow(maf_df), " mutation records across ", length(getSampleSummary(maf)$Tumor_Sample_Barcode), " samples.")

# 2. Cohort mutation summary
png(file.path(results_dir, "maf_summary.png"), width = 1600, height = 1200, res = 150)
plotmafSummary(maf, addStat = "median", dashboard = TRUE)
dev.off()

# 3. Top mutated genes (oncoplot)
png(file.path(results_dir, "oncoplot.png"), width = 1600, height = 1000, res = 150)
oncoplot(maf = maf, top = 20)
dev.off()

# 4. Tumor mutational burden
message("Calculating tumor mutational burden...")

# captureSize is the approximate exome capture size in Mb, used to normalize
# mutation counts into mutations/Mb (the standard TMB unit). 
tmb_result <- tmb(maf = maf, captureSize = 38, logScale = FALSE)

write.csv(tmb_result, file.path(results_dir, "tmb_summary.csv"), row.names = FALSE)
message("Median TMB: ", round(median(tmb_result$total_perMB, na.rm = TRUE), 2), " mutations/Mb")

# 5. Correlate TMB with immune signature scores
immune_path <- file.path(results_dir, "immune_signature_scores.csv")

if (file.exists(immune_path)) {
  message("Correlating TMB with immune signature scores from analysis.R...")

  immune_scores <- read.csv(immune_path)

  # TCGA barcodes: MAF uses full sample barcodes, immune scores used the
  # colnames from the expression data
  tmb_result$patient_barcode <- substr(tmb_result$Tumor_Sample_Barcode, 1, 12)

  # immune_scores includes both tumor AND normal samples 
  immune_sample_type <- substr(immune_scores$sample, 14, 15)
  immune_scores_tumor_only <- immune_scores[immune_sample_type == "01", ]
  immune_scores_tumor_only$patient_barcode <- substr(immune_scores_tumor_only$sample, 1, 12)

  # Explicitly deduplicate both sides to exactly one row per patient
  # before joining, rather than relying on the sample-type filter alone 
  tmb_result_dedup <- tmb_result %>% distinct(patient_barcode, .keep_all = TRUE)
  immune_scores_dedup <- immune_scores_tumor_only %>% distinct(patient_barcode, .keep_all = TRUE)

  n_tmb_dropped <- nrow(tmb_result) - nrow(tmb_result_dedup)
  n_immune_dropped <- nrow(immune_scores_tumor_only) - nrow(immune_scores_dedup)
  if (n_tmb_dropped > 0 || n_immune_dropped > 0) {
    message(
      "Deduplicated to one sample per patient before TMB-immune correlation: ",
      "dropped ", n_tmb_dropped, " duplicate TMB row(s), ",
      n_immune_dropped, " duplicate immune-score row(s)."
    )
  }

  merged <- inner_join(tmb_result_dedup, immune_scores_dedup, by = "patient_barcode")

  if (nrow(merged) >= 3) {
    # T cell infiltration is the most commonly examined axis for a
    # TMB-immune relationship.
    cor_test <- cor.test(merged$total_perMB, merged$`T.cell`, method = "spearman", exact = FALSE)

    corr_plot <- ggplot(merged, aes(x = total_perMB, y = `T.cell`)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
      theme_minimal(base_size = 13) +
      labs(
        title = paste0(project, ": TMB vs. T Cell Signature Score"),
        subtitle = paste0(
          "Spearman rho = ", round(cor_test$estimate, 3),
          ", p = ", format.pval(cor_test$p.value, digits = 3)
        ),
        x = "Tumor Mutational Burden (mutations/Mb)",
        y = "T cell signature score (mean marker z-score)"
      )

    ggsave(file.path(results_dir, "tmb_immune_correlation.png"), corr_plot, width = 7, height = 5, dpi = 150)
    message(
      "TMB vs. T cell correlation: rho = ", round(cor_test$estimate, 3),
      ", p = ", format.pval(cor_test$p.value, digits = 3)
    )
  } else {
    message("Fewer than 3 matched samples between TMB and immune scores; skipping correlation plot.")
  }
} else {
  message("No immune_signature_scores.csv found -- run analysis.R first to enable TMB-immune correlation.")
}

# 6. HSPA gene family (Hsp70) mutation frequency
message("Checking HSPA (Hsp70) gene family mutation frequency...")

HSPA_GENES <- c(
  "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5",
  "HSPA6", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14"
)

gene_summary <- getGeneSummary(maf)
hspa_mutation_freq <- gene_summary %>%
  filter(Hugo_Symbol %in% HSPA_GENES) %>%
  arrange(desc(MutatedSamples))

write.csv(hspa_mutation_freq, file.path(results_dir, "hspa_mutation_frequency.csv"), row.names = FALSE)

hspa_genes_present <- intersect(
  HSPA_GENES,
  gene_summary$Hugo_Symbol[gene_summary$MutatedSamples > 0]
)

message(
  length(hspa_genes_present), " of ", length(HSPA_GENES),
  " HSPA family genes carry at least one mutation in this cohort",
  if (length(hspa_genes_present) > 0) paste0(": ", paste(hspa_genes_present, collapse = ", ")) else "."
)

if (length(hspa_genes_present) > 0) {
  png(file.path(results_dir, "hspa_oncoplot.png"), width = 1400, height = 800, res = 150)
  oncoplot(maf = maf, genes = hspa_genes_present)
  dev.off()
  message("Saved HSPA-specific oncoplot to hspa_oncoplot.png")
} else {
  message("No HSPA family mutations found in this cohort; skipping HSPA oncoplot.")
}

message("Done. Mutation analysis results written to '", results_dir, "/'.")
