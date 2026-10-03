# TCGA Cancer Analysis

[![CI](https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis/actions/workflows/ci.yml)

Differential expression, pathway enrichment, immune signature scoring, mutation analysis, and survival analysis on public TCGA data. Reproducible end to end from a single script per module, with a combined HTML report.

## Pipeline

- **Sample selection** — one primary tumor (`01`) and at most one solid tissue normal (`11`) per patient; recurrent/metastatic samples and duplicate aliquots are excluded
- **Differential expression** — DESeq2, tumor vs. normal (or expression-based split for cohorts with too few normal samples)
- **Pathway enrichment** — fgsea against MSigDB Hallmark gene sets
- **Immune signature scoring** — marker-gene-based scoring across 9 immune/stromal populations
- **Survival analysis** — Kaplan-Meier and Cox models (stage-adjusted when available), stratified on the top DE gene, with a proportional-hazards check
- **Mutation analysis** — TMB, oncoplot, and TMB vs. T cell signature correlation (maftools), pooled and within each molecular subtype, plus immune-class composition per subtype
- **Immune-hot/cold classification** — tertiles (hot / intermediate / cold) of the Ayers 2017 IFN-γ signature
- **Subtype comparisons** — significant HSPA genes by TCGA molecular subtype with Kruskal-Wallis tests (BH-corrected)
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

## Gene sets

`analysis.R` looks for `genesets/h.all.v2023.2.Hs.symbols.gmt` first. If it is missing, Hallmark sets are downloaded via `msigdbr` (with retries) and written to that path, with the msigdbr version in the GMT description column. Commit the file to pin the gene set version for future runs. Alternatively, download the Hallmark GMT (gene symbols) from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) and save it there.

## Notes

- Default cohort is `TCGA-UCEC`; change the `project` variable in both `analysis.R` and `mutation_analysis.R` to analyze a different one.
- Immune signature scoring uses a lightweight marker-gene z-score approach rather than a dedicated deconvolution tool (xCell/ConsensusTME). Scores and hot/cold classes are relative to the cohort analyzed.
- Survival for the top DE gene is exploratory: being differentially expressed vs. normal tissue does not imply prognostic value. 
- Molecular subtypes come from the 2013 TCGA marker paper and cover only part of the cohort. "Not assigned" tumors are a mixed group; they are shown (in grey) but excluded from subtype-level statistics.

## Data source

Public TCGA RNA-seq and mutation data, queried live from the NCI Genomic Data Commons (GDC) via TCGAbiolinks.
