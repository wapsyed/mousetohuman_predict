# =============================================================================
# required.R - minimal helpers for the mouse-to-human prediction notebook.
# Contains ONLY what the notebook needs: package loading, theme_vaxgo() and
# autoGSEA(). (The full analysis pipeline uses scripts_notebooks/required.R in
# the main animals_vax_atlas repository.)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(janitor)         # clean_names()
  library(tidymodels)      # predict() on the fitted workflows
  library(plotly)          # ggplotly()
})

# clusterProfiler (used for the GSEA) lives on Bioconductor. Give a clear
# message if it is missing, instead of an opaque library() error.
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("The 'clusterProfiler' package is not installed.\n",
       "Run source('install.R') first, or install it with:\n",
       "  BiocManager::install('clusterProfiler')")
}
suppressPackageStartupMessages(library(clusterProfiler)) # GSEA()

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
