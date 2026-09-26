source(testthat::test_path("..", "..", "R", "functions.R"))

# determine_group()

test_that("determine_group uses tumor vs normal when enough normals present", {
  sample_types <- c("01", "01", "01", "11", "11", "11")
  expr_mat <- matrix(rnorm(6 * 10), nrow = 10, ncol = 6)
  colnames(expr_mat) <- paste0("sample", 1:6)

  result <- determine_group(sample_types, expr_mat, min_normal = 3)

  expect_equal(as.character(result), c("tumor", "tumor", "tumor", "normal", "normal", "normal"))
  expect_s3_class(result, "factor")
})

test_that("determine_group falls back to expression split when too few normals", {
  sample_types <- c("01", "01", "01", "01", "11")  # only 1 normal
  expr_mat <- matrix(0, nrow = 2, ncol = 5)
  rownames(expr_mat) <- c("MKI67", "OTHER")
  expr_mat["MKI67", ] <- c(10, 8, 2, 1, 5)  # median = 5

  result <- determine_group(sample_types, expr_mat, min_normal = 3)

  # marker_expr = c(10, 8, 2, 1, 5), median = 5; the function uses a strict
  # > comparison, so the value exactly equal to the median (index 5) is
  # correctly "low", not "high".
  expect_equal(as.character(result), c("high", "high", "low", "low", "low"))
})

test_that("determine_group falls back to colMeans when marker gene absent", {
  sample_types <- c("01", "01", "11")  # too few normals
  expr_mat <- matrix(c(1, 1, 9, 9, 1, 1), nrow = 2, ncol = 3)
  rownames(expr_mat) <- c("GENE_A", "GENE_B")  # no MKI67 present

  result <- determine_group(sample_types, expr_mat, min_normal = 5)

  expect_true(all(as.character(result) %in% c("high", "low")))
})

test_that("determine_group errors on mismatched lengths", {
  sample_types <- c("01", "01", "11")
  expr_mat <- matrix(0, nrow = 2, ncol = 5)  # 5 columns, only 3 sample types

  expect_error(determine_group(sample_types, expr_mat))
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
  fake_remote <- function() list("FAKE_SET" = c("GENEX", "GENEY"))

  result <- get_hallmark_sets(
    "/nonexistent/path/does_not_exist.gmt",
    remote_fetcher = fake_remote
  )

  expect_equal(result, list("FAKE_SET" = c("GENEX", "GENEY")))
})

test_that("get_hallmark_sets retries on remote failure and eventually errors", {
  attempt_count <- 0
  always_fails <- function() {
    attempt_count <<- attempt_count + 1
    stop("simulated network failure")
  }

  expect_error(
    get_hallmark_sets(
      "/nonexistent/path.gmt",
      max_attempts = 2,
      pause_seconds = 0,
      remote_fetcher = always_fails
    ),
    "Could not obtain MSigDB Hallmark gene sets"
  )
  expect_equal(attempt_count, 2)
})
