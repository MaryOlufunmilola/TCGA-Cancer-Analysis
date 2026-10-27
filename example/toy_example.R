#!/usr/bin/env Rscript
# Toy example: run the core logic in R/functions.R on a small synthetic cohort
# with known answers. Needs only base R and the survival package (no GDC
# download, no Bioconductor); takes a few seconds.
#
# Usage (from the repo root):  Rscript example/toy_example.R
# Outputs: example/output/*.csv. Exits with an error if any check fails.
#
# What is planted, and what each check expects:
#   - 60 patients: 60 primary tumors (01), 12 solid normals (11), plus 5
#     duplicate tumor aliquots and 3 metastatic samples (06) that sample
#     selection must drop.
#   - The first 20 tumors have high IFN-gamma signature genes -> "Immune-hot".
#   - GENE_RISK expression raises the death hazard -> Cox HR per SD > 1.

suppressPackageStartupMessages(library(survival))
source("R/functions.R")

set.seed(7)
out_dir <- "example/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Synthetic barcodes in TCGA format ----------------------------------------
n_patients <- 60
patients <- sprintf("TCGA-XX-%04d", seq_len(n_patients))
barcode <- function(patient, type, vial = "A", portion = "01") {
  paste0(patient, "-", type, vial, "-", portion, "R-A000-07")
}
barcodes <- c(
  barcode(patients, "01"),                              # one primary tumor each
  barcode(patients[1:12], "11"),                        # 12 solid normals
  barcode(patients[1:5], "01", vial = "B"),             # duplicate tumor aliquots
  barcode(patients[50:52], "06")                        # metastatic, must be excluded
)

keep <- select_one_sample_per_patient(barcodes, c("01", "11"))
kept <- barcodes[keep]
types <- sample_type_code(kept)
group <- determine_group(types, min_normal = 5)

# 2. Expression with a planted group of immune hot tumors ----------------------------
ifng_genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")  # Ayers 2017 IFN-gamma core genes
other_genes <- sprintf("GENE%03d", 1:100)
genes <- c(ifng_genes, "GENE_RISK", other_genes)
expr <- matrix(rnorm(length(genes) * length(kept), mean = 5, sd = 1),
               nrow = length(genes), dimnames = list(genes, kept))
tumor_idx <- which(types == "01")
hot_tumors <- kept[tumor_idx[1:20]]
expr[ifng_genes, hot_tumors] <- expr[ifng_genes, hot_tumors] + 3

tumor_expr <- expr[, tumor_idx]
scores <- score_immune_signatures(tumor_expr, list(IFNG_Ayers = ifng_genes))
immune_class <- classify_immune_hot_cold(scores[, "IFNG_Ayers"])

# 3. Survival with a planted risk gene ------------------------------------
risk <- tumor_expr["GENE_RISK", ]
hazard <- exp(1.2 * as.numeric(scale(risk)))   # true HR per SD = 3.3
death_time <- round(rexp(length(risk), rate = hazard / 1500)) + 1
censor_time <- round(runif(length(risk), 300, 3000))
clinical <- data.frame(
  patient = patient_barcode(colnames(tumor_expr)),
  vital_status = ifelse(death_time <= censor_time, "Dead", "Alive"),
  days_to_death = ifelse(death_time <= censor_time, death_time, NA),
  days_to_last_follow_up = ifelse(death_time <= censor_time, NA, censor_time),
  ajcc_pathologic_stage = sample(c("Stage I", "Stage IA", "Stage II", "Stage IIIC", "Stage IV"),
                                 length(risk), replace = TRUE),
  stringsAsFactors = FALSE
)
surv_df <- build_survival_data(clinical, risk)
surv_df$stage <- pick_stage(surv_df)
surv <- fit_survival_models(surv_df, stage_col = "stage")

# 4. Write outputs --------------------------------------------------------
sample_table <- data.frame(barcode = kept, sample_type = types, group = group)
write.csv(sample_table, file.path(out_dir, "selected_samples.csv"), row.names = FALSE)
write.csv(data.frame(barcode = rownames(scores), ifng_score = scores[, 1], immune_class = immune_class,
                     planted_hot = rownames(scores) %in% hot_tumors),
          file.path(out_dir, "immune_classes.csv"), row.names = FALSE)
write.csv(surv, file.path(out_dir, "survival_GENE_RISK.csv"), row.names = FALSE)

# 5. Checks against the planted answers -----------------------------------
check <- function(ok, what) {
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", what))
  ok
}
hot_recovered <- mean(immune_class[rownames(scores) %in% hot_tumors] == "Immune-hot")
results <- c(
  check(length(kept) == 72, sprintf("sample selection kept %d of %d (60 tumors + 12 normals)", length(kept), length(barcodes))),
  check(!any(types == "06"), "metastatic samples excluded"),
  check(identical(as.vector(table(group)), c(12L, 60L)), "12 normal vs. 60 tumor groups"),
  check(hot_recovered >= 0.9, sprintf("%.0f%% of planted hot tumors classed Immune-hot", 100 * hot_recovered)),
  check(surv$cox_hr_per_sd > 1 && surv$cox_p < 0.05,
        sprintf("GENE_RISK Cox HR per SD = %.2f (p = %.2g), expected > 1", surv$cox_hr_per_sd, surv$cox_p))
)
cat("\nSurvival summary:\n")
print(surv[, c("n", "n_events", "logrank_p", "cox_hr_per_sd", "cox_p", "adj_cox_hr_per_sd", "adj_cox_p")],
      row.names = FALSE)
cat("\nOutputs written to", out_dir, "\n")
if (!all(results)) stop("Toy example checks failed.")
cat("All toy example checks passed.\n")
