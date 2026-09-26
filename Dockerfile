# Reproducible environment for TCGA differential expression, pathway
# enrichment, immune deconvolution, and survival analysis.
#

FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL maintainer="Funmi Oyebamiji"
LABEL description="TCGA differential expression, GSEA, immune deconvolution, and survival analysis"

WORKDIR /analysis

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /usr/local/lib/R/etc && \
    echo "CXX11STD = -std=gnu++14" >> /usr/local/lib/R/etc/Makevars.site

RUN R -e " \
    pkgs <- c('TCGAbiolinks', 'SummarizedExperiment', 'DESeq2', 'fgsea'); \
    BiocManager::install(pkgs, update = FALSE, ask = FALSE); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Failed to install: ', paste(missing, collapse = ', ')) \
    "

RUN R -e " \
    pkgs <- c('survival', 'survminer', 'dplyr', 'tidyr', 'ggplot2', 'ggrepel', \
              'tibble', 'pheatmap', 'msigdbr', 'nnls', 'testthat', 'rmarkdown', 'knitr'); \
    install.packages(pkgs, repos = 'https://cloud.r-project.org'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Failed to install: ', paste(missing, collapse = ', ')) \
    "

RUN R -e " \
    BiocManager::install('maftools', update = FALSE, ask = FALSE); \
    if (!requireNamespace('maftools', quietly = TRUE)) stop('Failed to install maftools') \
    "

COPY analysis.R mutation_analysis.R deconvolution_custom_reference.R report.Rmd ./
COPY R/ ./R/
COPY tests/ ./tests/

# data/ and results/ are meant to be mounted from the host at run time 
RUN mkdir -p data results

CMD ["Rscript", "analysis.R"]
