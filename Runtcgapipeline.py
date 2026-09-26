#!/usr/bin/env bash
# TCGA-Cancer-Analysis: full setup + run sequence
#
# MUST run AFTER run_scrna_pipeline.sh -- step 5 below depends on
# ~/Projects/scRNA-seq-analysis/results/cell_type_signature_matrix.csv,
# which that script's own step 5 produces.
#
# Apple Silicon note: needs ARM-native scientific packages, so Miniforge
# (not standard Anaconda) is required. One-time setup, not per-repo:
#   curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
#   bash Miniforge3-MacOSX-arm64.sh
#   (restart your terminal after)

set -e  # stop on first error, rather than continuing past a failure

# ---------------------------------------------------------------------------
# 1. Project folder + environment
# ---------------------------------------------------------------------------
mkdir -p ~/Projects/TCGA-Cancer-Analysis && cd ~/Projects/TCGA-Cancer-Analysis

conda env create -f environment.yml
conda activate tcga-cancer-analysis

# ---------------------------------------------------------------------------
# 2. maftools -- installed separately, deliberately NOT in environment.yml.
# bioconda's maftools build is pinned to old R versions incompatible with
# r-base=4.4 used here.
# ---------------------------------------------------------------------------
R -e 'BiocManager::install("maftools", update = FALSE, ask = FALSE)'

# ---------------------------------------------------------------------------
# 3. Core analyses: differential expression, immune signatures, mutations/TMB
# ---------------------------------------------------------------------------
Rscript analysis.R 2>&1 | tee run_log.txt
Rscript mutation_analysis.R 2>&1 | tee run_log_mut.txt

# ---------------------------------------------------------------------------
# 4. Update environment for the deconvolution step (picks up nnls, added
# to environment.yml after the core analyses above were first built)
# ---------------------------------------------------------------------------
conda env update -f environment.yml --prune

# ---------------------------------------------------------------------------
# 5. Custom-reference deconvolution -- CROSS-REPO DEPENDENCY.
# Requires results/cell_type_signature_matrix.csv from the scRNA repo,
# already copied into data/ if you ran run_scrna_pipeline.sh first. If
# that file isn't in data/, this step will fail -- go run that script's
# step 5 before continuing here.
# ---------------------------------------------------------------------------
if [ ! -f data/cell_type_signature_matrix.csv ]; then
  echo "ERROR: data/cell_type_signature_matrix.csv not found."
  echo "Run run_scrna_pipeline.sh first (its step 5 produces and copies this file)."
  exit 1
fi
Rscript deconvolution_custom_reference.R

# ---------------------------------------------------------------------------
# 6. Ligand-receptor survival analysis -- the TFF3-CXCR4 signaling axis
# (motivated by the scRNA repo's cell-communication findings) tested
# against real clinical survival outcomes in this bulk cohort.
# ---------------------------------------------------------------------------
Rscript ligand_receptor_survival.R

# ---------------------------------------------------------------------------
# 7. Tests
# ---------------------------------------------------------------------------
Rscript tests/run_tests.R

# ---------------------------------------------------------------------------
# 8. Render the final report -- knits everything above into one HTML doc
# ---------------------------------------------------------------------------
Rscript -e 'rmarkdown::render("report.Rmd", params = list(project = "TCGA-UCEC"))'

echo ""
echo "Done. Open report.html and check results/ for all outputs."

# ---------------------------------------------------------------------------
# 9. Push to GitHub (uncomment and edit once the repo exists on github.com)
# ---------------------------------------------------------------------------
# mv README_SLIM.md README.md   # if not already renamed
# git init
# git add .
# git status   # confirm data/ and results/ are NOT listed -- .gitignore should exclude them
# git commit -m "Initial commit: full TCGA analysis pipeline"
# git branch -M main
# git remote add origin https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis.git
# git push -u origin main