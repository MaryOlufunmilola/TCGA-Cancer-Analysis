# TCGA Cancer Analysis

[![CI](https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/MaryOlufunmilola/TCGA-Cancer-Analysis/actions/workflows/ci.yml)

Differential expression, pathway enrichment, immune signature scoring, mutation
analysis, and survival analysis on public TCGA data, run here on endometrial cancer
(TCGA-UCEC). Each module is one R script, the core logic lives in tested functions in
`R/functions.R`, and `report.Rmd` combines the results into one HTML report.

## Quick start

```bash
conda env create -f environment.yml
conda activate tcga-cancer-analysis
R -e 'BiocManager::install("maftools", update = FALSE, ask = FALSE)'

bash run_all.sh               # tests, toy example, analysis.R, mutation_analysis.R
bash run_all.sh --report      # also renders report.html
```

Or run the pieces yourself:

```bash
Rscript analysis.R
Rscript mutation_analysis.R
Rscript -e 'rmarkdown::render("report.Rmd", params = list(project = "TCGA-UCEC"))'
```

The scripts download everything they need from the NCI Genomic Data Commons (GDC) on
the first run. To analyze a different cohort, change the `project` variable in
both `analysis.R` and `mutation_analysis.R`.

To check the core logic without downloading anything, run the toy example. It needs
only base R and the `survival` package:

```bash
Rscript example/toy_example.R
```

It builds a small synthetic cohort with known answers (60 patients with tumors and
normals, plus duplicate aliquots and metastatic samples that should be dropped, 20
planted immune hot tumors, and a planted risk gene) and checks that sample selection,
hot and cold classification, and the Cox model recover them. It ends with "All toy
example checks passed." CI runs it with the unit tests.

## What it does

- **Differential expression:** DESeq2, primary tumor vs. solid normal, keeping one
  primary tumor and at most one solid normal per patient (recurrent and metastatic
  samples and duplicate aliquots are dropped). If a cohort has too few normals, the grouping function returns no groups and a
  warning rather than quietly switching to a different comparison; a proliferation
  split is available as a separate function you call on purpose.
- **Pathway enrichment:** fgsea against MSigDB Hallmark gene sets.
- **Immune scoring:** marker z-score signatures for immune cell populations, and immune
  hot, intermediate, or cold classes from tertiles of the Ayers 2017 IFN-gamma signature.
  Scores and classes are relative to the cohort analyzed.
- **Survival:** Kaplan-Meier and log-rank test on the top gene, Cox models on
  z-scored expression (hazard ratio per 1 SD), with and without stage, and a
  proportional hazards check (`cox.zph`).
- **Mutations:** tumor mutational burden, oncoplot, and TMB vs. the T cell signature
  (maftools), pooled and within each molecular subtype, plus the mix of immune classes
  in each subtype.
- **Subtype comparisons:** significant HSPA genes compared across TCGA molecular
  subtypes with Kruskal-Wallis tests (BH corrected).
- **Report:** `report.Rmd` knits all of the above into one HTML file.

## Results for TCGA-UCEC

After keeping one primary tumor or solid normal per patient, the comparison has 545
tumors and 35 normals. 9,795 genes differ at an adjusted p below 0.05. As a check, the
23 patients with both a tumor and a normal sample were refit as matched pairs: fold
changes correlate well with the main analysis (Spearman 0.86), and 58% of its
significant genes stay significant.

The strongest pathway signals in tumors are cell cycle programs: E2F targets and the
G2M checkpoint (normalized enrichment scores 3.4 and 3.0). CDC20, a cell cycle gene,
has the smallest adjusted p value (log2 fold change 4.3).

Immune context differs by molecular subtype. Among tumors with a subtype call, 47% of
POLE tumors (8 of 17) are immune hot, compared with 14% of copy number low tumors, which are
the most often immune cold (49%).

For survival, 497 primary tumors with follow up (90 deaths) were analyzed for CDC20:

| Model | Result |
|---|---|
| Kaplan-Meier, median split | log-rank p = 0.077 |
| Cox, CDC20 per 1 SD | HR 1.26 (95% CI 1.01 to 1.57), p = 0.040; proportional hazards check p = 0.16 |
| Cox, adjusted for stage | HR 1.16, p = 0.30 |

The unadjusted association does not hold after adjusting for stage. Proliferation
tracks with stage, and the median split was not prespecified, so I treat this as
exploratory rather than as a prognostic marker.

## Inputs and outputs

All inputs are downloaded from the GDC by TCGAbiolinks into `data/`: RNA-seq counts
with clinical data (STAR counts), the open access somatic mutation MAF, and the TCGA
molecular subtypes. Hallmark gene sets are read from
`genesets/h.all.v2023.2.Hs.symbols.gmt`; if that file is missing, they are fetched
with `msigdbr` and saved there, so committing the file pins the gene set version.

Main outputs in `results/`:

| File | What it is |
|---|---|
| `differential_expression_results.csv`, `volcano_plot.png` | DESeq2 tumor vs. normal, all genes |
| `paired_sensitivity_results.csv`, `paired_sensitivity_comparison.csv` | The matched pair check |
| `gsea_hallmark_results.csv`, `gsea_hallmark_plot.png` | Hallmark enrichment |
| `immune_signature_scores.csv`, `immune_signature_heatmap.png`, `immune_hot_cold_plot.png` | Immune scores and classes for each tumor |
| `survival_top_gene_results.csv`, `survival_km_plot.png` | Kaplan-Meier and Cox results |
| `tmb_summary.csv`, `oncoplot.png`, `tmb_immune_*` | Mutational burden and its link to immune signatures |
| `sessionInfo.txt` | Package versions for the run |

## Run time

`run_all.sh` records each step's run time and peak memory in `logs/timing.tsv`. The
first run of `analysis.R` is dominated by the GDC download. The toy example takes a
few seconds.

## Tests

```bash
Rscript tests/run_tests.R
```

Unit tests (testthat) cover the functions in `R/functions.R`, with external calls
mocked. CI runs them and the toy example on every push, installing only `testthat`
and `survival` so it stays fast.

## Structure

```
analysis.R              DE, GSEA, immune scoring, survival
mutation_analysis.R     TMB, oncoplot, TMB and immune correlation
report.Rmd              combined HTML report
run_all.sh              runs everything with timing
R/functions.R           core logic, unit tested
tests/                  testthat unit tests
example/toy_example.R   synthetic cohort with known answers
```

## Limitations

Immune scores use a simple marker z-score approach, not a dedicated deconvolution tool
such as xCell or ConsensusTME, and they are relative to the tumors in the cohort.
Molecular subtypes come from the 2013 TCGA marker paper and cover only part of the
cohort; "Not assigned" tumors are a mixed group, shown in plots but left out of
subtype level statistics. The survival
analysis looks at one gene chosen from the DE results, so it is a starting point, not
a validated model.

## Data source

Public TCGA RNA-seq, clinical, and mutation data, queried from the NCI Genomic Data
Commons through TCGAbiolinks.
