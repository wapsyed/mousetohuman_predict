# R packages needed by RunYourAnalysis_MouseToHuman.Rmd
#
# This file is executed by Binder (and can be run locally with source("install.R")).

# --- CRAN packages -----------------------------------------------------------
install.packages(c("here", "tidyverse", "janitor", "tidymodels", "plotly", "ranger"))

# --- Bioconductor packages ---------------------------------------------------
# clusterProfiler (used for the GSEA) lives on Bioconductor.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Pin the Bioconductor release that matches R 4.5.3 (Bioc 3.22) and install
# clusterProfiler together with its dependencies.
BiocManager::install(
  "clusterProfiler",
  ask = FALSE,
  update = FALSE,
  version = "3.22"
)

# Fail the build early (instead of at library() time) if it did not install.
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("clusterProfiler failed to install - check the build log.")
}
