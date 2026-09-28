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
invisible(dev.off())

# 3. Top mutated genes (oncoplot)
png(file.path(results_dir, "oncoplot.png"), width = 1600, height = 1000, res = 150)
oncoplot(maf = maf, top = 20)
invisible(dev.off())

# 4. Tumor mutational burden
message("Calculating tumor mutational burden...")

# captureSize is the approximate exome capture size in Mb, used to normalize
# mutation counts into mutations/Mb (the standard TMB unit). 38 Mb is a common
# approximation for TCGA exome capture; treat absolute values as approximate.
tmb_result <- tmb(maf = maf, captureSize = 38, logScale = FALSE)

write.csv(tmb_result, file.path(results_dir, "tmb_summary.csv"), row.names = FALSE)
message(
  "Median TMB: ", round(median(tmb_result$total_perMB, na.rm = TRUE), 2), " mutations/Mb ",
  "(approximate -- assumes a 38 Mb exome capture size; absolute values would shift under a ",
  "different assumed capture size, though relative comparisons across samples are unaffected)."
)

# 5. Correlate TMB with immune signature scores
immune_path <- file.path(results_dir, "immune_signature_scores.csv")

if (file.exists(immune_path)) {
  message("Correlating TMB with immune signature scores from analysis.R...")

  immune_scores <- read.csv(immune_path)
  immune_scores$hot_cold <- as_hot_cold_factor(immune_scores$hot_cold)  # CSV loses factor order

  # analysis.R already writes primary tumors only, one per patient. Apply the
  # same rule to the MAF side so both are one primary tumor per patient.
  tmb_tumor <- tmb_result[sample_type_code(tmb_result$Tumor_Sample_Barcode) == "01", ]
  tmb_result_dedup <- tmb_tumor[select_one_sample_per_patient(tmb_tumor$Tumor_Sample_Barcode, "01"), ]
  tmb_result_dedup$patient_barcode <- patient_barcode(tmb_result_dedup$Tumor_Sample_Barcode)

  immune_scores <- immune_scores[sample_type_code(immune_scores$sample) == "01", ]
  immune_scores_dedup <- immune_scores[select_one_sample_per_patient(immune_scores$sample, "01"), ]
  immune_scores_dedup$patient_barcode <- patient_barcode(immune_scores_dedup$sample)

  n_tmb_dropped <- nrow(tmb_result) - nrow(tmb_result_dedup)
  if (n_tmb_dropped > 0) {
    message("Dropped ", n_tmb_dropped, " non-primary or duplicate TMB row(s) before correlation.")
  }

  # Defensive checks before correlating -- barcode format mismatches or
  # upstream column changes should fail loudly and clearly here, not
  # produce a cryptic error (or worse, a silently wrong result) downstream.
  stopifnot(
    "Duplicate patient_barcode in tmb_result_dedup after dedup -- barcode format may differ from expected TCGA length" =
      !any(duplicated(tmb_result_dedup$patient_barcode)),
    "Duplicate patient_barcode in immune_scores_dedup after dedup" =
      !any(duplicated(immune_scores_dedup$patient_barcode))
  )

  # as.data.frame: maftools returns data.tables; plain data frames keep dplyr/ggplot behaviour predictable
  merged <- as.data.frame(inner_join(as.data.frame(tmb_result_dedup), immune_scores_dedup, by = "patient_barcode"))

  if (!("T.cell" %in% colnames(merged))) {
    stop(
      "Expected column 'T.cell' not found in immune_signature_scores.csv (columns present: ",
      paste(colnames(merged), collapse = ", "), "). Check that analysis.R's immune_marker_sets ",
      "still includes a 'T cell' entry with that exact name."
    )
  }

  # Drop rows with missing TMB or T cell score explicitly, rather than relying
  # on cor.test()'s own (less transparent) handling of NA values.
  n_before_na_filter <- nrow(merged)
  merged <- merged[!is.na(merged$total_perMB) & !is.na(merged$`T.cell`), , drop = FALSE]
  n_dropped_na <- n_before_na_filter - nrow(merged)
  if (n_dropped_na > 0) {
    message("Dropped ", n_dropped_na, " sample(s) with missing TMB or T cell score before correlation.")
  }

  if (nrow(merged) >= 3) {
    # T cell infiltration is the most commonly examined axis for a
    # TMB-immune relationship.
    cor_test <- cor.test(merged$total_perMB, merged$`T.cell`, method = "spearman", exact = FALSE)

    # TMB is heavily right-skewed (POLE-ultramutated and MSI-H tumors), so plot
    # it on a log axis. 
    corr_plot <- ggplot(merged, aes(x = total_perMB, y = `T.cell`)) +
      geom_point(aes(color = hot_cold), alpha = 0.6) +
      scale_x_log10() +
      scale_color_manual(values = HOT_COLD_COLORS) +
      theme_minimal(base_size = 13) +
      labs(
        title = paste0(project, ": TMB vs. T Cell Signature Score"),
        subtitle = paste0(
          "Spearman rho = ", round(cor_test$estimate, 3),
          ", p = ", format.pval(cor_test$p.value, digits = 3),
          "\nColor = cohort-relative immune class (IFN-gamma signature tertiles), not an absolute phenotype"
        ),
        x = "Tumor Mutational Burden (mutations/Mb, log scale)",
        color = "IFN-gamma class",
        y = "T cell signature score (mean marker z-score)"
      )

    ggsave(file.path(results_dir, "tmb_immune_correlation.png"), corr_plot, width = 7, height = 5, dpi = 150)
    message(
      "TMB vs. T cell correlation (", nrow(merged), " primary tumors): rho = ", round(cor_test$estimate, 3),
      ", p = ", format.pval(cor_test$p.value, digits = 3)
    )

    # 5b. Same relationship split by molecular subtype. In UCEC, POLE-ultramutated
    # and MSI-H tumors drive most of the high-TMB tail, so the pooled correlation
    # can mix a between-subtype effect with any within-subtype effect.
    subtype_lookup <- load_subtype_lookup(sub("^TCGA-", "", project))

    if (!is.null(subtype_lookup)) {
      merged_sub <- merged %>%
        inner_join(subtype_lookup, by = "patient_barcode")
      write.csv(
        merged_sub %>% select(patient_barcode, molecular_subtype, total_perMB, T.cell, IFN.gamma, hot_cold),
        file.path(results_dir, "tmb_immune_by_subtype.csv"), row.names = FALSE
      )

      # Per-subtype Spearman correlation (assigned subtypes with >= 10 tumors).
      # "Not assigned" is a mixture that includes unlabeled POLE/MSI tumors, so it
      # is excluded from the statistics rather than treated as a subtype.
      merged_assigned <- merged_sub %>% filter(molecular_subtype != "Not assigned")

      subtype_cor <- merged_assigned %>%
        group_by(molecular_subtype) %>%
        filter(n() >= 10) %>%
        summarise(
          n = n(),
          median_tmb_per_Mb = median(total_perMB),
          spearman_rho = suppressWarnings(cor(total_perMB, T.cell, method = "spearman")),
          spearman_p = suppressWarnings(cor.test(total_perMB, T.cell, method = "spearman", exact = FALSE)$p.value),
          .groups = "drop"
        ) %>%
        mutate(
          spearman_padj = p.adjust(spearman_p, method = "BH"),
          # n >= 10 is the minimum to compute a correlation at all here, not
          # a threshold for it being reliable -- flag anything still small
          # (arbitrary but reasonable cutoff at 20) so a reader doesn't read
          # a striking rho from 11 tumors with the same confidence as one
          # from 80.
          small_n_caution = n < 20
        ) %>%
        arrange(desc(median_tmb_per_Mb))
      write.csv(subtype_cor, file.path(results_dir, "tmb_immune_subtype_correlation.csv"), row.names = FALSE)

      # Immune class composition per subtype ("Not assigned" kept here, listed last)
      class_by_subtype <- merged_sub %>%
        mutate(molecular_subtype = order_subtypes(molecular_subtype)) %>%
        count(molecular_subtype, hot_cold, .drop = FALSE) %>%
        group_by(molecular_subtype) %>%
        mutate(fraction = n / sum(n)) %>%
        ungroup()
      write.csv(class_by_subtype, file.path(results_dir, "immune_class_by_subtype.csv"), row.names = FALSE)

      # One panel per assigned subtype, ordered by median TMB, with rho in the
      # panel title. 
      plot_df <- merged_assigned %>% filter(molecular_subtype %in% subtype_cor$molecular_subtype)
      panel_labels <- setNames(
        paste0(subtype_cor$molecular_subtype, ifelse(subtype_cor$small_n_caution, "*", ""),
               " (n = ", subtype_cor$n,
               ", \u03c1 = ", sprintf("%.2f", subtype_cor$spearman_rho),
               ", padj = ", formatC(subtype_cor$spearman_padj, format = "g", digits = 2), ")"),
        subtype_cor$molecular_subtype
      )
      plot_df$molecular_subtype <- factor(plot_df$molecular_subtype, levels = subtype_cor$molecular_subtype)

      subtype_tmb_plot <- ggplot(plot_df, aes(x = total_perMB, y = T.cell, color = hot_cold)) +
        geom_point(alpha = 0.7) +
        facet_wrap(~molecular_subtype, labeller = as_labeller(panel_labels)) +
        scale_x_log10() +
        scale_color_manual(values = HOT_COLD_COLORS, drop = FALSE) +
        theme_minimal(base_size = 12) +
        labs(
          title = paste0(project, ": TMB vs. T Cell Signature within Molecular Subtypes"),
          subtitle = paste0(
            "Spearman correlation per assigned subtype (>= 10 tumors), BH-corrected; 'Not assigned' excluded.\n",
            "* = fewer than 20 tumors, interpret cautiously. Color = cohort-relative immune class, not an absolute phenotype."
          ),
          x = "Tumor Mutational Burden (mutations/Mb, log scale)",
          y = "T cell signature score (mean marker z-score)",
          color = "IFN-gamma class"
        )
      ggsave(file.path(results_dir, "tmb_immune_by_subtype.png"), subtype_tmb_plot, width = 10, height = 6.5, dpi = 150)

      message(
        "TMB-immune by subtype: ", nrow(merged_assigned), " tumors with an assigned subtype (",
        sum(merged_sub$molecular_subtype == "Not assigned"), " 'Not assigned' excluded); within-subtype rho: ",
        paste0(subtype_cor$molecular_subtype, " = ", round(subtype_cor$spearman_rho, 2), collapse = ", ")
      )
    }
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
n_tumors_assessed <- length(getSampleSummary(maf)$Tumor_Sample_Barcode)

hspa_mutation_freq <- gene_summary %>%
  filter(Hugo_Symbol %in% HSPA_GENES) %>%
  mutate(
    tumors_assessed = n_tumors_assessed,
    pct_tumors_mutated = round(100 * MutatedSamples / n_tumors_assessed, 1)
  ) %>%
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
  invisible(dev.off())
  message("Saved HSPA-specific oncoplot to hspa_oncoplot.png")
} else {
  message("No HSPA family mutations found in this cohort; skipping HSPA oncoplot.")
}

message("Done. Mutation analysis results written to '", results_dir, "/'.")
