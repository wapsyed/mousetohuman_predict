# R packages needed by RunYourAnalysis_MouseToHuman.Rmd
install.packages(c("here", "tidyverse", "janitor", "tidymodels", "plotly"))

# GSEA (clusterProfiler) comes from Bioconductor
install.packages("BiocManager")
BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
