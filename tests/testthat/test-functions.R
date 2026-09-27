source(testthat::test_path("..", "..", "R", "functions.R"))

# determine_group() and determine_proliferation_group()
# NOTE: an earlier version of determine_group() took an expr_mat argument
# and silently fell back to an MKI67-based split when too few normals were
# present. That fallback was split out into a separate function
# (determine_proliferation_group()) -- see review point 2 -- so the tests
# below are rewritten to match the current, narrower signature and behavior.

test_that("determine_group uses tumor vs normal when enough normals present", {
  sample_types <- c("01", "01", "01", "11", "11", "11")

  result <- determine_group(sample_types, min_normal = 3)

  expect_equal(as.character(result), c("tumor", "tumor", "tumor", "normal", "normal", "normal"))
  expect_s3_class(result, "factor")
})

test_that("determine_group returns all-NA with a warning when too few normals", {
  sample_types <- c("01", "01", "01", "01", "11")  # only 1 normal

  expect_warning(
    result <- determine_group(sample_types, min_normal = 3),
    "Only 1 normal sample"
  )
  expect_true(all(is.na(result)))
  expect_equal(length(result), length(sample_types))
})

test_that("determine_group excludes sample types other than 01/11", {
  sample_types <- c("01", "11", "06")  # "06" = metastatic, should become NA

  result <- determine_group(sample_types, min_normal = 1)

  expect_equal(as.character(result), c("tumor", "normal", NA))
})

test_that("determine_proliferation_group splits primary tumors on the marker gene median", {
  sample_types <- c("01", "01", "01", "01", "11")  # last one is a normal, should get NA
  expr_mat <- matrix(0, nrow = 1, ncol = 5, dimnames = list("MKI67", NULL))
  expr_mat["MKI67", ] <- c(10, 8, 2, 1, 999)  # median among tumors (10,8,2,1) = 5

  result <- determine_proliferation_group(sample_types, expr_mat)

  expect_equal(as.character(result), c("high", "high", "low", "low", NA))
})

test_that("determine_proliferation_group falls back to colMeans when marker gene absent", {
  sample_types <- c("01", "01", "11")
  expr_mat <- matrix(c(1, 1, 9, 9, 1, 1), nrow = 2, ncol = 3)
  rownames(expr_mat) <- c("GENE_A", "GENE_B")  # no MKI67 present

  expect_warning(
    result <- determine_proliferation_group(sample_types, expr_mat),
    "not found"
  )
  expect_true(all(as.character(result)[sample_types == "01"] %in% c("high", "low")))
})

test_that("determine_proliferation_group errors on mismatched lengths", {
  sample_types <- c("01", "01", "11")
  expr_mat <- matrix(0, nrow = 2, ncol = 5)  # 5 columns, only 3 sample types

  expect_error(determine_proliferation_group(sample_types, expr_mat))
})

# score_immune_signatures()

test_that("score_immune_signatures gives higher scores to samples with elevated markers", {
  genes <- c("CD3D", "CD3E", "OTHER1", "OTHER2")
  log_counts <- matrix(1, nrow = 4, ncol = 3, dimnames = list(genes, c("s1", "s2", "s3")))
  # s1 has strongly elevated T cell markers; s2, s3 do not
  log_counts["CD3D", "s1"] <- 20
  log_counts["CD3E", "s1"] <- 20

  marker_sets <- list("T cell" = c("CD3D", "CD3E"))
  scores <- score_immune_signatures(log_counts, marker_sets)

  expect_true(scores["s1", "T cell"] > scores["s2", "T cell"])
  expect_true(scores["s1", "T cell"] > scores["s3", "T cell"])
})

test_that("score_immune_signatures returns NA for a marker set with no matching genes", {
  log_counts <- matrix(1:12, nrow = 4, ncol = 3, dimnames = list(
    c("GENE_A", "GENE_B", "GENE_C", "GENE_D"), c("s1", "s2", "s3")
  ))
  marker_sets <- list("Nonexistent" = c("NOT_PRESENT_1", "NOT_PRESENT_2"))

  scores <- score_immune_signatures(log_counts, marker_sets)

  expect_true(all(is.na(scores[, "Nonexistent"])))
})

test_that("score_immune_signatures output has one row per sample", {
  log_counts <- matrix(rnorm(20), nrow = 5, ncol = 4)
  rownames(log_counts) <- paste0("gene", 1:5)
  colnames(log_counts) <- paste0("sample", 1:4)
  marker_sets <- list("A" = c("gene1", "gene2"), "B" = c("gene3"))

  scores <- score_immune_signatures(log_counts, marker_sets)

  expect_equal(nrow(scores), 4)
  expect_equal(ncol(scores), 2)
})

test_that("score_immune_signatures handles a zero-variance gene without producing NaN scores", {
  log_counts <- matrix(rnorm(20), nrow = 5, ncol = 4)
  rownames(log_counts) <- paste0("gene", 1:5)
  colnames(log_counts) <- paste0("sample", 1:4)
  log_counts["gene1", ] <- 5  # identical in every sample -> zero variance -> NaN from scale()

  marker_sets <- list("A" = c("gene1", "gene2"))  # gene1 (zero-var) + gene2 (real variance)

  expect_message(
    scores <- score_immune_signatures(log_counts, marker_sets),
    "zero-variance gene"
  )
  # na.rm = TRUE in the underlying colMeans should mean the score is still a
  # real number (from gene2 alone), not NaN just because gene1 is degenerate.
  expect_true(all(!is.nan(scores[, "A"])))
})

# fit_survival_models() -- proportional hazards check

test_that("fit_survival_models includes a cox.zph proportional-hazards p-value", {
  set.seed(1)
  n <- 60
  surv_df <- data.frame(
    time = rexp(n, rate = 0.01),
    event = rbinom(n, 1, 0.7),
    expr = rnorm(n),
    expr_group = factor(sample(c("Low", "High"), n, replace = TRUE), levels = c("Low", "High"))
  )

  result <- fit_survival_models(surv_df)

  expect_true("cox_zph_p" %in% colnames(result))
  expect_true(is.na(result$cox_zph_p) || (result$cox_zph_p >= 0 && result$cox_zph_p <= 1))
})

# pick_subtype_columns() -- explicit-first selection

test_that("pick_subtype_columns uses the explicit preferred column when present", {
  cols <- c("patient", "Subtype_Integrative", "some_other_cluster_col")

  result <- pick_subtype_columns(cols, preferred_subtype_col = "Subtype_Integrative")

  expect_equal(result$subtype, "Subtype_Integrative")
  expect_equal(result$patient, "patient")
})

test_that("pick_subtype_columns falls back to auto-detection when the preferred column is absent", {
  cols <- c("patient", "Some_Other_Subtype_Field")

  result <- pick_subtype_columns(cols, preferred_subtype_col = "Subtype_Integrative")

  expect_equal(result$subtype, "Some_Other_Subtype_Field")
})

# get_hallmark_sets()

test_that("get_hallmark_sets reads a local GMT file when present", {
  tmp_gmt <- tempfile(fileext = ".gmt")
  writeLines(c(
    "HALLMARK_TEST_SET\tdescription\tGENE1\tGENE2\tGENE3",
    "HALLMARK_OTHER_SET\tdescription\tGENE4\tGENE5"
  ), tmp_gmt)

  # Injected fake reader avoids requiring the real fgsea package for this test
  fake_reader <- function(path) {
    lines <- readLines(path)
    setNames(
      lapply(lines, function(l) strsplit(l, "\t")[[1]][-c(1, 2)]),
      sapply(lines, function(l) strsplit(l, "\t")[[1]][1])
    )
  }

  result <- get_hallmark_sets(tmp_gmt, gmt_reader = fake_reader)

  expect_true("HALLMARK_TEST_SET" %in% names(result))
  expect_equal(result[["HALLMARK_TEST_SET"]], c("GENE1", "GENE2", "GENE3"))

  unlink(tmp_gmt)
})

test_that("get_hallmark_sets falls back to remote fetcher when local file absent", {
  tmp_gmt <- tempfile(fileext = ".gmt")

  fake_remote <- function() {
    list("FAKE_SET" = c("GENEX", "GENEY"))
  }

  result <- get_hallmark_sets(
    tmp_gmt,
    remote_fetcher = fake_remote
  )

  expect_equal(result, list("FAKE_SET" = c("GENEX", "GENEY")))
  expect_true(file.exists(tmp_gmt))

  unlink(tmp_gmt)
})

test_that("get_hallmark_sets retries on remote failure and eventually errors", {
  tmp_gmt <- tempfile(fileext = ".gmt")
  attempt_count <- 0

  always_fails <- function() {
    attempt_count <<- attempt_count + 1
    stop("simulated network failure")
  }

  expect_error(
    get_hallmark_sets(
      tmp_gmt,
      max_attempts = 2,
      pause_seconds = 0,
      remote_fetcher = always_fails
    ),
    "Could not obtain MSigDB Hallmark gene sets"
  )

  expect_equal(attempt_count, 2)

  unlink(tmp_gmt)
})
