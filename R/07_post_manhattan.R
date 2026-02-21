here::i_am("R/07_post_manhattan.R")

####################
# Run ORA on gene lists of genes with DMRs
####################

# Import packages
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(purrr)
library(tidyverse)
library(vctrs)

# Load data
results_mcsea <- readRDS(
  file = file.path(
    "output",
    "data_dmr",
    "results_mcsea_withmean.RDS"
  )
)

titles <- c(
  "P+13 vs. P+5 (Control D3) - Promoters",
  "P+13 vs. P+5 (Control D3) - Gene body",
  "Control vs. Hep @P+5 D3 - Promoters",
  "Control vs. Hep @P+5 D3 - Gene body",
  "Control vs. Hep @P+13 D3 - Promoters",
  "Control vs. Hep @P+13 D3 - Gene body"
)

goi <- list(
  timepoint_promoters = c(
    "CDK7",
    "FOS",
    "FOSL2",
    "NFKBIA",
    "IFNGR2",
    "MAPT"
  ),
  timepoint_genes = c(
    "GPC5",
    "HS3ST4",
    "BMP4"
  ),
  treatment_early_promoters = c(),
  treatment_early_genes = c(),
  treatment_late_promoters = c(),
  treatment_late_genes = c()
)

# Filter for necessary data
## Annotate with chromosome and location
## Select only relevant dataframes
results_mcsea <- results_mcsea %>%
  map(\(results) {
    # Get annotation (chromosome and loc for each gene)
    annotation <- AnnotationDbi::select(
      EnsDb.Hsapiens.v86,
      keys = rownames(results),
      columns = c(
        "SEQNAME",
        "GENESEQSTART"
      ),
      keytype = "GENENAME"
    ) %>%
      dplyr::filter(
        SEQNAME %in%
          c(
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X"
          )
      ) %>%
      mutate(
        SEQNAME = SEQNAME %>%
          factor(levels = c(1:22, "X"))
      ) %>%
      group_by(GENENAME) %>%
      slice(1)

    # Combine
    results_out <- results %>%
      rownames_to_column(var = "gene") %>%
      left_join(
        annotation,
        by = join_by("gene" == "GENENAME")
      ) %>%
      drop_na() %>%
      dplyr::relocate(
        gene,
        SEQNAME,
        GENESEQSTART,
        padj
      )

    # Return data
    return(results_out)
  })

# Create manhattan plots
plots_manhattan <- pmap(
  list(
    results_mcsea,
    goi,
    titles
  ),
  \(results, goi, title) {
    results %>%
      mutate(
        label = case_when(
          padj <= 0.004 ~ gene,
          gene %in% goi ~ gene,
          .default = NA
        )
      ) %>%
      manhattan_plot(
        chr.colname = "SEQNAME",
        pval.colname = "padj",
        pos.colname = "GENESEQSTART",
        label.colname = "label",
        rescale = FALSE,
        signif = c(0.05),
        preserve.position = FALSE,
        plot.title = title,
        max.overlaps = 25,
        ylim = c(3, NA),
        force = 4,
        force_pull = 0.3,
      )
  }
)

grid_manhattan <- plots_manhattan %>%
  wrap_plots() +
  plot_layout(ncol = 1)

# Create volcano plots
plots_volcano2d <- map(
  results_mcsea,
  \(results) {
    EnhancedVolcano(
      results,
      lab = results$gene,
      x = "leadingEdge_mean",
      y = "padj",
      xlim = c(-0.6, 0.6),
      ylim = c(0, -log10(1e-13)),
      xlab = "Difference of beta values",
      ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
      axisLabSize = 8,
      title = NULL,
      titleLabSize = 8,
      subtitle = NULL,
      caption = NULL,
      pCutoff = 0.05,
      FCcutoff = 0.1,
      pointSize = 1.5,
      labSize = 2,
      boxedLabels = FALSE,
      legendPosition = "none",
      drawConnectors = TRUE,
      widthConnectors = 0.1,
      colConnectors = "grey60",
      arrowheads = FALSE,
      gridlines.major = FALSE,
      gridlines.minor = FALSE
    )
  }
)

grid_volcano2d <- plots_volcano2d %>%
  wrap_plots() +
  plot_layout(ncol = 1)

# Create combined grid
grid_timepoint <-
  wrap_plots(
    plots_manhattan[[1]],
    plots_volcano2d[[1]],
    plots_manhattan[[2]],
    plots_volcano2d[[2]]
  ) +
  plot_layout(
    ncol = 2,
    widths = c(3, 1)
  )

grid_treatment <- wrap_plots(
  plots_manhattan[[3]],
  plots_volcano2d[[3]],
  plots_manhattan[[4]],
  plots_volcano2d[[4]],
  plots_manhattan[[5]],
  plots_volcano2d[[5]],
  plots_manhattan[[6]],
  plots_volcano2d[[6]]
) +
  plot_layout(
    ncol = 2,
    widths = c(3, 1)
  )

# Export plots
## Manhattan plots
iwalk(
  plots_manhattan,
  \(x, idx) {
    ggsave(
      filename = str_c("manhattan_", idx, ".png"),
      plot = x,
      path = here::here(
        "output",
        "plots_explore"
      ),
      scale = 0.8,
      width = 8,
      height = 4,
      units = "in",
      dpi = 288
    )
  }
)

ggsave(
  filename = "manhattan_grid.png",
  plot = grid_manhattan,
  path = here::here(
    "output",
    "plots_explore"
  ),
  scale = 0.9,
  width = 8,
  height = 12,
  units = "in",
  dpi = 288
)

## Volcano plots
iwalk(
  plots_volcano2d,
  \(x, idx) {
    ggsave(
      filename = str_c("volcano2d_", idx, ".png"),
      plot = x,
      path = here::here(
        "output",
        "plots_explore"
      ),
      scale = 0.8,
      width = 8,
      height = 4,
      units = "in",
      dpi = 288
    )
  }
)

ggsave(
  filename = "volcano2d_grid.png",
  plot = grid_volcano2d,
  path = here::here(
    "output",
    "plots_explore"
  ),
  scale = 0.9,
  width = 8,
  height = 12,
  units = "in",
  dpi = 288
)

## Combined plots
ggsave(
  filename = "combined_timepoint.png",
  plot = grid_timepoint,
  path = here::here(
    "output",
    "plots_explore"
  ),
  scale = 1.1,
  width = 8,
  height = 5,
  units = "in",
  dpi = 288
)

ggsave(
  filename = "combined_treatment.png",
  plot = grid_treatment,
  path = here::here(
    "output",
    "plots_explore"
  ),
  scale = 1.1,
  width = 8,
  height = 10,
  units = "in",
  dpi = 288
)
