# =============================================================================
# required.R - packages + helpers for the mouse-to-human prediction notebook.
# Installs any missing package automatically (CRAN + Bioconductor), loads them,
# and defines theme_vaxgo() and autoGSEA(). Works locally and on Binder.
# (The full analysis pipeline uses scripts_notebooks/required.R in the main
# animals_vax_atlas repository.)
# =============================================================================

required_cran <- c("here", "tidyverse", "janitor", "tidymodels", "ranger", "glmnet", "plotly", "ggrepel", "patchwork")
required_bioc <- c("clusterProfiler")

# CRAN mirror fallback (avoids an interactive prompt when none is set)
if (is.null(getOption("repos")) || is.na(getOption("repos")["CRAN"]) ||
    getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# --- Install anything missing (CRAN) -----------------------------------------
missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

# --- Install anything missing (Bioconductor) ---------------------------------
# BiocManager picks the Bioconductor release that matches the running R version.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# --- Load --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(janitor)         # clean_names()
  library(tidymodels)      # predict() on the fitted workflows
  library(ranger)          # engine used by the saved random-forest models
  library(glmnet)          # engine used by the saved Lasso model (score_shared)
  library(patchwork)       # plot_layout() / wrap_plots()
  library(ggrepel)         # geom_text_repel()
  library(plotly)          # ggplotly()
  library(clusterProfiler) # GSEA()
})

# --- Plot theme --------------------------------------------------------------
theme_vaxgo <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      # Geral
      ggh4x.facet.nestline = ggplot2::element_line(colour = "black", linetype = 1),

      # Axis
      axis.text            = ggplot2::element_text(size = 10),
      axis.text.x          = ggplot2::element_text(size = 10, color = "black", angle = 0),
      axis.text.y          = ggplot2::element_text(size = 10, color = "black", angle = 0),
      axis.title.x         = ggplot2::element_text(color = "black", face = "bold"),
      axis.title.y         = ggplot2::element_text(color = "black", face = "bold"),
      axis.line.x          = ggplot2::element_line(linewidth = 0.5, colour = "black", linetype = 1),
      axis.line.y          = ggplot2::element_line(linewidth = 0.5, colour = "black", linetype = 1),
      axis.ticks.x         = ggplot2::element_line(linewidth = 0.5, color = "black"),
      axis.ticks.y         = ggplot2::element_line(linewidth = 0.5, color = "black"),

      # Legend
      legend.position      = "right",
      legend.location      = "plot",
      legend.text          = ggplot2::element_text(size = 10),
      legend.title         = ggplot2::element_text(size = 10),
      legend.key.width     = grid::unit(0.4, "cm"),
      legend.key.height    = grid::unit(0.4, "cm"),
      legend.margin        = margin(0, 0, 0, 0),
      legend.box.margin    = margin(0, 0, -10, 0),

      # Panel
      panel.border         = ggplot2::element_blank(),
      panel.grid.major.x   = ggplot2::element_blank(),
      panel.grid.major.y   = ggplot2::element_blank(),
      panel.grid.minor     = ggplot2::element_blank(),
      panel.spacing        = grid::unit(0.2, "cm"),

      # Strip
      strip.text           = ggplot2::element_text(size = 10, color = "black", margin = margin(1, 0, 1, 0)),

      # Plot
      plot.title           = ggplot2::element_text(size = 10),
      plot.subtitle        = ggplot2::element_text(size = 10),
      plot.caption         = ggplot2::element_text(hjust = 0, size = 5),
      plot.margin          = ggplot2::margin(0.1, 0.2, 0.1, 0.1, "cm")
    )
}

# --- GSEA for every condition in the table -----------------------------------
autoGSEA <- function(df, TERM2GENE, geneset_name) {

  gsea_results <- list()

  conditions <- df$condition %>% unique() %>% as.character()

  for (condition_i in conditions) {

    degs_condition <- df %>%
      filter(condition == condition_i) %>%
      select(genes, rank) %>%
      distinct() %>%
      arrange(desc(rank)) %>%
      deframe()

    auto_gsea <- tryCatch({

      GSEA(
        geneList = degs_condition,
        TERM2GENE = TERM2GENE,
        minGSSize = 1,
        maxGSSize = 1000,
        pvalueCutoff = 1,
        pAdjustMethod = "BH"
      ) %>%
        as.data.frame() %>%
        arrange(qvalue) %>%
        mutate(
          condition = condition_i,
          gsea_enrichment = geneset_name
        )

    }, error = function(e) NULL)

    if (!is.null(auto_gsea)) {
      key <- paste(condition_i, geneset_name, sep = "_")
      gsea_results[[key]] <- auto_gsea
    }
  }

  return(list(gsea = bind_rows(gsea_results)))
}
