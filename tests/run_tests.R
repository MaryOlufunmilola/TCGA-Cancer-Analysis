#!/usr/bin/env Rscript
# Entry point for running the unit test suite.
# Usage: Rscript tests/run_tests.R

library(testthat)

results <- test_dir("tests/testthat", reporter = "summary")

# Exit non-zero on any failure so CI catches it
if (any(as.data.frame(results)$failed > 0)) {
  quit(status = 1)
}
