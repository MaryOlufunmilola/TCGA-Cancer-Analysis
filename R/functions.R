# Core, unit-tested logic for the TCGA analysis scripts.
# Everything here depends only on base R, stats, and survival (a recommended
# package shipped with R), so the test suite runs without Bioconductor.

# TCGA barcode helpers ----------------------------------------------------

#' Two-digit TCGA sample type code (e.g. "01" primary tumor, "11" solid normal)
sample_type_code <- function(barcodes) substr(barcodes, 14, 15)

#' Patient barcode (first 12 characters, e.g. "TCGA-AX-A1CR")
patient_barcode <- function(barcodes) substr(barcodes, 1, 12)

#' Pick one sample per patient for the given sample type code(s)
#'
#' TCGA often has several aliquots/portions per patient. Keeping all of them
#' double-counts patients in DE, survival, and correlation analyses. This keeps
#' the lexicographically first barcode per patient (vial A, lowest portion),
#' which is deterministic across runs.
#'
#' @param barcodes character vector of full TCGA sample barcodes
#' @param sample_codes sample type codes to keep, e.g. "01" or c("01", "11")
#' @return sorted integer indices into `barcodes`
select_one_sample_per_patient <- function(barcodes, sample_codes = "01") {
  barcodes <- as.character(barcodes)  # maftools may return factors
  idx <- which(sample_type_code(barcodes) %in% sample_codes)
  idx <- idx[order(barcodes[idx])]
  # Deduplicate per patient *and* sample type, so a patient keeps one tumor
  # and one normal when both codes are requested.
  key <- paste(patient_barcode(barcodes[idx]), sample_type_code(barcodes[idx]))
  sort(idx[!duplicated(key)])
}

# Comparison groups -------------------------------------------------------

#' Determine tumor vs. normal comparison groups
#'
#' Uses primary tumors ("01") and solid tissue normals ("11"); other sample
#' types (recurrent, metastatic, ...) get NA so callers can drop them.
#'
#' Factor levels are set explicitly so the reference level is predictable:
#' c("normal", "tumor").
#'
#' IMPORTANT: this function no longer silently substitutes a different
#' comparison (MKI67-high/low) when too few normals are available. An
#' earlier version did exactly that -- but "insufficient normal samples for
#' tumor-vs-normal" and "I want a proliferation-based tumor-only split" are
#' two independent, unrelated decisions, and conflating them meant a caller
#' could get back a completely different biological comparison (proliferation
#' status, not tumor-vs-normal) without any code-level signal that the
#' question being answered had silently changed. If too few normals exist,
#' this now returns all-NA with a clear warning, so callers see the failure
#' explicitly rather than an unannounced change of analysis. If you want the
#' proliferation-based split, call determine_proliferation_group() directly
#' -- it's now a separate, deliberately-invoked function (see below), not an
#' automatic fallback.
#'
#' @param sample_types character vector of TCGA sample type codes
#' @param min_normal minimum number of normal samples required
#' @return factor, same length as sample_types, possibly containing NA;
#'   entirely NA (with a warning) if fewer than min_normal normals exist
determine_group <- function(sample_types, min_normal = 5) {
  is_tumor  <- sample_types == "01"
  is_normal <- sample_types == "11"

  if (sum(is_normal) < min_normal) {
    warning(
      "Only ", sum(is_normal), " normal sample(s) found (need >= ", min_normal,
      "). Returning all-NA rather than silently substituting a different ",
      "comparison -- call determine_proliferation_group() directly if a ",
      "tumor-only proliferation split is what you actually want."
    )
    return(factor(rep(NA_character_, length(sample_types)), levels = c("normal", "tumor")))
  }

  group <- ifelse(is_normal, "normal", ifelse(is_tumor, "tumor", NA_character_))
  factor(group, levels = c("normal", "tumor"))
}

#' Determine MKI67-high/low proliferation groups among primary tumors only
#'
#' A DELIBERATE, standalone analysis (not a fallback for anything) -- splits
#' primary tumors by expression of a proliferation marker (default MKI67),
#' at the median among tumors. This answers "which tumors are highly
#' proliferative" -- a genuinely different question from tumor-vs-normal,
#' and DE results from this split should be described as reflecting
#' proliferation status, not general tumor biology.
#'
#' @param sample_types character vector of TCGA sample type codes
#' @param expr_mat genes x samples matrix of depth-normalized expression,
#'   gene symbols as rownames, columns aligned with sample_types
#' @param marker_gene marker gene (symbol) for the high/low split
#' @return factor, same length as sample_types; NA for non-tumor samples
determine_proliferation_group <- function(sample_types, expr_mat, marker_gene = "MKI67") {
  stopifnot(length(sample_types) == ncol(expr_mat))
  is_tumor <- sample_types == "01"

  if (marker_gene %in% rownames(expr_mat)) {
    marker_expr <- expr_mat[marker_gene, ]
  } else {
    warning("Marker gene '", marker_gene, "' not found; splitting on mean expression instead.")
    marker_expr <- colMeans(expr_mat)
  }
  cutoff <- stats::median(marker_expr[is_tumor])
  group <- ifelse(!is_tumor, NA_character_, ifelse(marker_expr > cutoff, "high", "low"))
  factor(group, levels = c("low", "high"))
}

# Pathway enrichment ------------------------------------------------------

#' Build a named, sorted ranking vector for preranked GSEA keyed by gene symbol
#'
#' Gene set collections such as MSigDB Hallmark use gene symbols, so ranks must
#' be named by symbol rather than Ensembl ID. When several Ensembl IDs map to
#' one symbol, the entry with the largest |stat| is kept.
#'
#' @param res_df data frame with a gene_name column and a ranking statistic
#' @param stat_col name of the ranking column (e.g. the Wald statistic)
#' @return named numeric vector sorted in decreasing order
ranks_by_symbol <- function(res_df, stat_col = "stat") {
  stat <- res_df[[stat_col]]
  ok <- !is.na(stat) & !is.na(res_df$gene_name) & nzchar(res_df$gene_name)
  df <- data.frame(gene_name = res_df$gene_name[ok], stat = stat[ok], stringsAsFactors = FALSE)
  df <- df[order(-abs(df$stat)), ]
  df <- df[!duplicated(df$gene_name), ]
  sort(stats::setNames(df$stat, df$gene_name), decreasing = TRUE)
}

#' Write gene sets to a GMT file (name, description, genes; tab-separated)
#'
#' @param gene_sets named list of character vectors
#' @param path output path; parent directories are created
#' @param description text for the second GMT column (e.g. source + version)
write_gmt <- function(gene_sets, path, description = "NA") {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  lines <- vapply(names(gene_sets), function(nm) {
    paste(c(nm, description, unique(gene_sets[[nm]])), collapse = "\t")
  }, character(1))
  writeLines(lines, path)
  invisible(path)
}

#' Get Hallmark gene sets, preferring a local GMT file over a live download
#'
#' @param local_path path to a local .gmt file
#' @param max_attempts retries for the live-download fallback
#' @param pause_seconds delay between retries
#' @param gmt_reader function used to read a local GMT file
#' @param remote_fetcher function used to fetch gene sets remotely
#' @param cache if TRUE, a successful remote fetch is written to local_path so
#'   later runs use the same gene set version without network access
#' @return named list of character vectors (gene sets)
get_hallmark_sets <- function(local_path,
                              max_attempts = 3,
                              pause_seconds = 10,
                              gmt_reader = fgsea::gmtPathways,
                              remote_fetcher = NULL,
                              cache = TRUE) {
  if (file.exists(local_path)) {
    message("Using local Hallmark gene set file: ", local_path)
    return(gmt_reader(local_path))
  }

  if (is.null(remote_fetcher)) {
    remote_fetcher <- function() {
      # msigdbr >= 10 renamed `category` to `collection`
      sets <- if (utils::packageVersion("msigdbr") >= "10.0.0") {
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
      } else {
        msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      }
      split(sets$gene_symbol, sets$gs_name)
    }
  }

  message(
    "No local gene set file found at '", local_path, "'. Falling back to ",
    "a live download (see README, 'Gene sets', to avoid this in future runs)."
  )
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(remote_fetcher(), error = function(e) {
      message("  Attempt ", attempt, " failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(result)) {
      if (cache) {
        desc <- if (requireNamespace("msigdbr", quietly = TRUE)) {
          paste0("msigdbr_", utils::packageVersion("msigdbr"))
        } else {
          "remote"
        }
        write_gmt(result, local_path, description = desc)
        message("Cached gene sets to '", local_path, "' (", desc, "). Commit this file to pin the version.")
      }
      return(result)
    }
    if (attempt < max_attempts) {
      message("  Retrying in ", pause_seconds, " seconds...")
      Sys.sleep(pause_seconds)
    }
  }
  stop(
    "Could not obtain MSigDB Hallmark gene sets: no local file found and ",
    "the remote fetch failed ", max_attempts, " times."
  )
}

# Immune signatures -------------------------------------------------------

#' Score samples against marker gene sets via mean per-gene z-score
#'
#' Z-scores are relative to the samples passed in, so callers should pass a
#' biologically coherent set (e.g. primary tumors only), not tumors + normals.
#'
#' A zero-variance gene (identical expression in every sample -- e.g. never
#' detected) produces NaN from scale() for every sample; na.rm = TRUE below
#' already prevents this from corrupting the per-sample mean (that one gene
#' is just excluded from the average for genes that have it), but previously
#' this happened silently with no record of it. Now reported via a message.
#'
#' @param log_counts genes x samples log-scale expression matrix (symbols as rownames)
#' @param marker_sets named list of character vectors (gene symbols per signature)
#' @return samples x signatures matrix
score_immune_signatures <- function(log_counts, marker_sets) {
  gene_zscores <- t(scale(t(log_counts)))

  zero_var_genes <- rownames(log_counts)[apply(log_counts, 1, stats::sd) == 0]
  if (length(zero_var_genes) > 0) {
    affected_sets <- names(marker_sets)[vapply(marker_sets, function(g) any(g %in% zero_var_genes), logical(1))]
    if (length(affected_sets) > 0) {
      message(
        length(zero_var_genes), " zero-variance gene(s) excluded from signature scoring ",
        "(affects: ", paste(affected_sets, collapse = ", "), "): ",
        paste(intersect(zero_var_genes, unique(unlist(marker_sets))), collapse = ", ")
      )
    }
  }

  scores <- sapply(marker_sets, function(genes) {
    genes_present <- intersect(genes, rownames(gene_zscores))
    if (length(genes_present) == 0) {
      return(rep(NA_real_, ncol(gene_zscores)))
    }
    colMeans(gene_zscores[genes_present, , drop = FALSE], na.rm = TRUE)
  })
  # sapply returns a vector, not a matrix, when there is a single sample
  scores <- matrix(scores, nrow = ncol(log_counts), dimnames = list(colnames(log_counts), names(marker_sets)))
  scores
}

#' Classify tumors as immune-hot / intermediate / immune-cold by score tertiles
#'
#' Uses a single inflammation score (e.g. the Ayers IFN-gamma signature) and
#' cohort tertiles, instead of comparing two differently-composed z-scores.
#'
#' @param score numeric vector of per-sample scores (NA allowed)
#' @param probs lower and upper quantile cutoffs
#' @return factor with levels Immune-cold, Intermediate, Immune-hot
classify_immune_hot_cold <- function(score, probs = c(1 / 3, 2 / 3)) {
  stopifnot(length(probs) == 2, probs[1] < probs[2])
  q <- stats::quantile(score, probs, na.rm = TRUE, names = FALSE)
  out <- ifelse(score >= q[2], "Immune-hot",
                ifelse(score <= q[1], "Immune-cold", "Intermediate"))
  factor(out, levels = c("Immune-cold", "Intermediate", "Immune-hot"))
}

# Survival ----------------------------------------------------------------

#' Build an analysis-ready survival data frame
#'
#' Dead patients use days_to_death; living patients use days_to_last_follow_up.
#' Rows with missing/non-positive time or unknown vital status are dropped, and
#' the High/Low split uses the median among the rows that remain.
#'
#' @param clinical data frame with vital_status, days_to_death, days_to_last_follow_up
#' @param expr_values numeric expression values aligned with clinical rows
#' @return clinical with added columns expr, expr_group, time, event (filtered)
build_survival_data <- function(clinical, expr_values) {
  stopifnot(nrow(clinical) == length(expr_values))

  status <- as.character(clinical$vital_status)
  event <- ifelse(status == "Dead", 1L, ifelse(status == "Alive", 0L, NA_integer_))
  death <- suppressWarnings(as.numeric(clinical$days_to_death))
  follow <- suppressWarnings(as.numeric(clinical$days_to_last_follow_up))
  time <- ifelse(event == 1L, death, follow)

  out <- clinical
  out$expr <- as.numeric(expr_values)
  out$time <- time
  out$event <- event
  out <- out[!is.na(out$time) & out$time > 0 & !is.na(out$event) & !is.na(out$expr), , drop = FALSE]

  cutoff <- stats::median(out$expr)
  out$expr_group <- factor(ifelse(out$expr > cutoff, "High", "Low"), levels = c("Low", "High"))
  out
}

#' Collapse detailed stage strings ("Stage IIIC", "Stage IA") to I-IV
collapse_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  # Order matters: "STAGE IV" and "STAGE III" also contain "STAGE I"
  out[grepl("STAGE IV", x)] <- "IV"
  out[is.na(out) & grepl("STAGE III", x)] <- "III"
  out[is.na(out) & grepl("STAGE II", x)] <- "II"
  out[is.na(out) & grepl("STAGE I", x)] <- "I"
  factor(out, levels = c("I", "II", "III", "IV"))
}

#' Find the best-populated stage column and collapse it to I-IV
#'
#' @param clinical data frame of clinical annotations
#' @param candidates column names to consider, in order of preference
#' @param min_fraction minimum fraction of non-missing collapsed values
#' @return factor aligned with clinical rows, or NULL if no usable column
pick_stage <- function(clinical,
                       candidates = c("ajcc_pathologic_stage", "figo_stage",
                                      "ajcc_clinical_stage", "clinical_stage"),
                       min_fraction = 0.5) {
  best <- NULL
  best_frac <- 0
  for (col in intersect(candidates, colnames(clinical))) {
    stage <- collapse_stage(clinical[[col]])
    frac <- mean(!is.na(stage))
    if (frac > best_frac) {
      best <- stage
      best_frac <- frac
      attr(best, "source_column") <- col
    }
  }
  if (is.null(best) || best_frac < min_fraction || nlevels(droplevels(best)) < 2) return(NULL)
  best
}

#' Log-rank test on the median split plus Cox models on continuous expression
#'
#' Also runs cox.zph() on each fitted Cox model to check the proportional
#' hazards assumption -- previously absent. A Cox model's hazard ratio and
#' p-value assume the covariate's effect is constant over time; if that
#' assumption is violated (cox.zph global p < 0.05), the reported HR is an
#' average over time that can be misleading, and this should be flagged
#' alongside the estimate, not silently omitted.
#'
#' @param surv_df output of build_survival_data()
#' @param stage_col optional name of a stage factor column for an adjusted Cox model
#' @return one-row data frame of test statistics, including cox.zph p-values
fit_survival_models <- function(surv_df, stage_col = NULL) {
  Surv <- survival::Surv  # local binding so formulas below resolve it
  surv_df$expr_z <- as.numeric(scale(surv_df$expr))

  lr <- survival::survdiff(Surv(time, event) ~ expr_group, data = surv_df)
  logrank_p <- stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

  cox_fit <- survival::coxph(Surv(time, event) ~ expr_z, data = surv_df)
  cox <- summary(cox_fit)
  # cox.zph() can itself fail on edge cases (e.g. too few events) -- fail
  # safe with NA rather than aborting the whole survival analysis for one gene.
  zph_p <- tryCatch(
    survival::cox.zph(cox_fit)$table["GLOBAL", "p"],
    error = function(e) NA_real_
  )

  out <- data.frame(
    n = nrow(surv_df),
    n_events = sum(surv_df$event),
    logrank_p = logrank_p,
    cox_hr_per_sd = cox$conf.int["expr_z", "exp(coef)"],
    cox_hr_lower95 = cox$conf.int["expr_z", "lower .95"],
    cox_hr_upper95 = cox$conf.int["expr_z", "upper .95"],
    cox_p = cox$coefficients["expr_z", "Pr(>|z|)"],
    cox_zph_p = zph_p,
    adj_cox_hr_per_sd = NA_real_,
    adj_cox_p = NA_real_,
    adj_cox_zph_p = NA_real_,
    adjusted_for = NA_character_
  )

  if (!is.null(stage_col) && stage_col %in% colnames(surv_df)) {
    adj_df <- surv_df[!is.na(surv_df[[stage_col]]), , drop = FALSE]
    adj_df[[stage_col]] <- droplevels(adj_df[[stage_col]])
    if (nlevels(adj_df[[stage_col]]) >= 2 && sum(adj_df$event) >= 10) {
      f <- stats::as.formula(paste("Surv(time, event) ~ expr_z +", stage_col))
      adj_fit <- survival::coxph(f, data = adj_df)
      adj <- summary(adj_fit)
      adj_zph_p <- tryCatch(
        survival::cox.zph(adj_fit)$table["GLOBAL", "p"],
        error = function(e) NA_real_
      )
      out$adj_cox_hr_per_sd <- adj$conf.int["expr_z", "exp(coef)"]
      out$adj_cox_p <- adj$coefficients["expr_z", "Pr(>|z|)"]
      out$adj_cox_zph_p <- adj_zph_p
      out$adjusted_for <- "stage"
    }
  }
  out
}

# Shared plotting constants ----------------------------------------------

HOT_COLD_LEVELS <- c("Immune-cold", "Intermediate", "Immune-hot")
HOT_COLD_COLORS <- c("Immune-cold" = "steelblue", "Intermediate" = "grey60", "Immune-hot" = "firebrick")

#' Restore the immune class factor order (lost when round-tripping through CSV)
as_hot_cold_factor <- function(x) factor(as.character(x), levels = HOT_COLD_LEVELS)

# Molecular subtypes -----------------------------------------------------

#' Pick patient and subtype columns from a TCGAquery_subtype() table
#'
#' Tries an explicit, known-good column name FIRST (if it's actually present
#' in this table), falling back to auto-detection only if it isn't. An
#' earlier version used auto-detection as the only mechanism -- convenient,
#' and it does log which column it picked, but for a portfolio repo an
#' explicit, reproducible default is preferable to "whichever column the
#' regex happened to match this time."
#'
#' @param cols column names of the subtype table
#' @param preferred_subtype_col the known column name to try first, e.g.
#'   "Subtype_Integrative" for the 2013 UCEC marker paper's table
#' @return list(patient = , subtype = ) or NULL if neither can be found
pick_subtype_columns <- function(cols, preferred_subtype_col = "Subtype_Integrative") {
  patient <- intersect(c("patient", "Patient", "bcr_patient_barcode"), cols)
  if (length(patient) == 0) return(NULL)

  if (preferred_subtype_col %in% cols) {
    return(list(patient = patient[1], subtype = preferred_subtype_col))
  }

  candidates <- grep("subtype|cluster|integrat", cols, ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0) return(NULL)
  # Prefer a column with "subtype" in the name over a generic "cluster" match
  preferred <- grep("subtype", candidates, ignore.case = TRUE, value = TRUE)
  list(patient = patient[1], subtype = if (length(preferred)) preferred[1] else candidates[1])
}

#' Tidy subtype labels: blanks and "NA" become NA; "Notassigned" gets a space
clean_subtype_labels <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA")] <- NA
  x[!is.na(x) & grepl("^not[ _]?assigned$", x, ignore.case = TRUE)] <- "Not assigned"
  x
}

#' One row per patient with a cleaned molecular subtype
#'
#' @param tumor_code TCGA code without prefix, e.g. "UCEC"
#' @param fetcher function(tumor) returning the subtype table (injectable for tests)
#' @param preferred_subtype_col explicit column name to try first -- see pick_subtype_columns()
#' @return data.frame(patient_barcode, molecular_subtype) or NULL on failure
load_subtype_lookup <- function(tumor_code, fetcher = NULL, preferred_subtype_col = "Subtype_Integrative") {
  if (is.null(fetcher)) fetcher <- function(tumor) TCGAbiolinks::TCGAquery_subtype(tumor = tumor)
  info <- tryCatch(fetcher(tumor_code), error = function(e) {
    message("Could not retrieve molecular subtype info for ", tumor_code, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(info)) return(NULL)
  info <- as.data.frame(info)

  cols <- pick_subtype_columns(colnames(info), preferred_subtype_col = preferred_subtype_col)
  if (is.null(cols)) {
    message("Could not auto-detect patient/subtype columns. Available: ", paste(colnames(info), collapse = ", "))
    return(NULL)
  }
  message("Using subtype column '", cols$subtype, "' (patient column '", cols$patient, "').")

  out <- data.frame(
    patient_barcode = as.character(info[[cols$patient]]),
    molecular_subtype = clean_subtype_labels(info[[cols$subtype]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$molecular_subtype), , drop = FALSE]
  out[!duplicated(out$patient_barcode), , drop = FALSE]
}

#' Kruskal-Wallis p-value across groups, ignoring excluded labels
#'
#' @param values numeric vector
#' @param groups group labels aligned with values
#' @param exclude labels left out of the test (e.g. "Not assigned")
#' @return p-value, or NA if fewer than two groups remain
kruskal_p <- function(values, groups, exclude = "Not assigned") {
  keep <- !is.na(values) & !is.na(groups) & !(groups %in% exclude)
  if (length(unique(groups[keep])) < 2) return(NA_real_)
  stats::kruskal.test(values[keep], factor(groups[keep]))$p.value
}

#' Order subtype labels alphabetically with "Not assigned" last
order_subtypes <- function(x, not_assigned = "Not assigned") {
  lv <- sort(unique(stats::na.omit(as.character(x))))
  factor(as.character(x), levels = c(setdiff(lv, not_assigned), intersect(not_assigned, lv)))
}

#' Named colors for subtype levels, with "Not assigned" in grey
subtype_palette <- function(levels, not_assigned = "Not assigned") {
  assigned <- setdiff(levels, not_assigned)
  cols <- stats::setNames(grDevices::hcl.colors(max(length(assigned), 1), "Dark 3")[seq_along(assigned)], assigned)
  if (not_assigned %in% levels) cols[not_assigned] <- "grey65"
  cols[levels]
}
