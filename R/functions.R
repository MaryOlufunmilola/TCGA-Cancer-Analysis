#' Determine comparison groups (tumor vs. normal, or expression-based fallback)
#'
#' @param sample_types character vector of TCGA sample type codes (e.g. "01", "11")
#' @param expr_mat genes x samples expression matrix, columns aligned with sample_types
#' @param min_normal minimum number of normal samples required to use tumor vs. normal
#' @param marker_gene fallback marker gene for expression-based high/low split
#' @return factor of group labels, same length as sample_types
determine_group <- function(sample_types, expr_mat, min_normal = 5, marker_gene = "MKI67") {
  stopifnot(length(sample_types) == ncol(expr_mat))

  n_normal <- sum(sample_types == "11")

  if (n_normal >= min_normal) {
    group <- ifelse(sample_types == "11", "normal", "tumor")
  } else {
    if (marker_gene %in% rownames(expr_mat)) {
      marker_expr <- expr_mat[marker_gene, ]
    } else {
      marker_expr <- colMeans(expr_mat)
    }
    group <- ifelse(marker_expr > median(marker_expr), "high", "low")
  }

  factor(group)
}

#' Score samples against marker gene sets via mean per-gene z-score
#'
#' @param log_counts genes x samples log-normalized expression matrix
#' @param marker_sets named list of character vectors (gene symbols per cell type)
#' @return samples x cell-type-scores matrix
score_immune_signatures <- function(log_counts, marker_sets) {
  gene_zscores <- t(scale(t(log_counts)))

  scores <- sapply(marker_sets, function(genes) {
    genes_present <- intersect(genes, rownames(gene_zscores))
    if (length(genes_present) == 0) {
      return(rep(NA_real_, ncol(gene_zscores)))
    }
    colMeans(gene_zscores[genes_present, , drop = FALSE], na.rm = TRUE)
  })

  rownames(scores) <- colnames(log_counts)
  scores
}

#' Get Hallmark gene sets, preferring a local GMT file over a live download
#'
#' @param local_path path to a local .gmt file
#' @param max_attempts retries for the live-download fallback
#' @param pause_seconds delay between retries
#' @param gmt_reader function used to read a local GMT file 
#' @param remote_fetcher function used to fetch gene sets remotely 
#' @return named list of character vectors (gene sets)
get_hallmark_sets <- function(local_path,
                               max_attempts = 3,
                               pause_seconds = 10,
                               gmt_reader = fgsea::gmtPathways,
                               remote_fetcher = NULL) {
  if (file.exists(local_path)) {
    message("Using local Hallmark gene set file: ", local_path)
    return(gmt_reader(local_path))
  }

  if (is.null(remote_fetcher)) {
    remote_fetcher <- function() {
      sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
      split(sets$gene_symbol, sets$gs_name)
    }
  }

  message(
    "No local gene set file found at '", local_path, "'. Falling back to ",
    "a live download (see README for how to avoid this in future runs)."
  )
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(remote_fetcher(), error = function(e) {
      message("  Attempt ", attempt, " failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(result)) return(result)
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
