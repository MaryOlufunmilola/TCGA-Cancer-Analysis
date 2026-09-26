# TCGA Cancer Analysis

[![CI](https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis/actions/workflows/ci.yml)

Differential expression, pathway enrichment, immune signature scoring, mutation analysis, and survival analysis on public TCGA data. Reproducible end to end from a single script per module, with a combined HTML report.

## Pipeline

- **Differential expression** — DESeq2, tumor vs. normal (or expression-based split for cohorts with too few normal samples)
- **Pathway enrichment** — fgsea against MSigDB Hallmark gene sets
- **Immune signature scoring** — marker-gene-based scoring across 9 immune/stromal populations
- **Survival analysis** — Kaplan-Meier, stratified on the top DE gene
- **Mutation analysis** — TMB, oncoplot, and TMB vs. immune signature correlation (maftools)
- **Immune-hot/cold classification** — reuses the cell-type marker panel from my `scRNA-seq-analysis` repo to classify tumor immune context
- **HTML report** — knits all of the above into one document

## Quickstart

```bash
conda env create -f environment.yml
conda activate tcga-cancer-analysis
R -e 'BiocManager::install("maftools", update = FALSE, ask = FALSE)'

Rscript analysis.R
Rscript mutation_analysis.R
```

Render the report:

```r
rmarkdown::render("report.Rmd", params = list(project = "TCGA-UCEC"))
```

## Testing

```bash
Rscript tests/run_tests.R
```

## Structure

```
analysis.R              # DE, GSEA, immune scoring, survival
mutation_analysis.R      # TMB, oncoplot, TMB-immune correlation
report.Rmd               # combined HTML report
R/functions.R            # core logic, unit tested
tests/                   # testthat unit tests
```

## Notes

- Default cohort is `TCGA-UCEC`; change the `project` variable in `analysis.R` to analyze a different one.
- Immune signature scoring uses a lightweight marker-gene z-score approach rather than a dedicated deconvolution tool (xCell/ConsensusTME).

## Data source

Public TCGA RNA-seq and mutation data, queried live from the NCI Genomic Data Commons (GDC) via TCGAbiolinks.
