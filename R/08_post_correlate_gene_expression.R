here::i_am("R/08_post_correlate_gene_expression.R")

####################
# Build 2x2 way plots to correlate DNA methylation statistics with logFC
####################

# Import packages
library(cowplot)
library(DESeq2)
library(ggpmisc)
library(minfi)
library(patchwork)
library(purrr)
library(tidyverse)
library(vctrs)

# Load data
## Transcriptome
txome_quant_deseq2_batchcor <- readRDS(
  file = here::here(
    "input",
    "txome",
    "quant_deseq2_batchcor_nofilter.RDS"
  )
)

txome_results <- readRDS(
  file = here::here(
    "input",
    "txome",
    "results_deseq2_nofilter.RDS"
  )
)

## DNA methylation
results_mcsea <- readRDS(
  file = file.path(
    "output",
    "data_dmr",
    "results_mcsea_withmean.RDS"
  )
)

# Definitions
axis_titles <- list(
  timepoint_promoters = c(
    "P+5 vs. P+13 - Promoters",
    "P+5 vs. P+13"
  ),
  timepoint_genes = c(
    "P+5 vs. P+13 - Gene Body",
    "P+5 vs. P+13"
  ),
  treatment_early_promoters = c(
    "Hep vs. Ctrl @ P+5 - Promoters",
    "Hep vs. Ctrl @ P+5"
  ),
  treatment_early_genes = c(
    "Hep vs. Ctrl @ P+5 - Gene Body",
    "Hep vs. Ctrl @ P+5"
  ),
  treatment_late_promoters = c(
    "Hep vs. Ctrl @ P+13 - Promoters",
    "Hep vs. Ctrl @ P+13"
  ),
  treatment_late_genes = c(
    "Hep vs. Ctrl @ P+13 - Gene Body",
    "Hep vs. Ctrl @ P+13"
  )
)

axis_titles_poi <- list(
  axis_titles[[1]],
  axis_titles[[2]],
  axis_titles[[1]],
  axis_titles[[2]],
  axis_titles[[1]],
  axis_titles[[2]]
)

# Compile appropriate fold changes and deltabetas
# List of comparisons
# P+5 vs. P+13 UT D3 (P13vsP5_UT_D3)
# Treatment at P+5 D3 (Trt_P5_D3)
# Treatment at P+13 D3 (Trt_P13_D3)
## txome data
txome_results_filter <- list(
  timepoint_promoters = txome_results[["P13vsP5_UT_D3"]],
  timepoint_genes = txome_results[["P13vsP5_UT_D3"]],
  treatment_early_promoters = txome_results[["Trt_P5_D3"]],
  treatment_early_genes = txome_results[["Trt_P5_D3"]],
  treatment_late_promoters = txome_results[["Trt_P13_D3"]],
  treatment_late_genes = txome_results[["Trt_P13_D3"]]
) %>%
  map(
    \(x) {
      x %>%
        data.frame() %>%
        rownames_to_column(var = "gene_id") %>%
        left_join(
          rowRanges(txome_quant_deseq2_batchcor) %>%
            data.frame() %>%
            dplyr::select(gene_id, gene_name, entrezid),
          by = join_by(gene_id == gene_id)
        )
    }
  )

## DNAm data
mcsea_results_filter <- results_mcsea %>%
  map(\(x) rownames_to_column(x, var = "gene_name"))

# Merge lists
# Define colours based on signif level
results_merge <- map2(
  mcsea_results_filter,
  txome_results_filter,
  \(mcsea, deseqresults) {
    left_join(
      mcsea,
      deseqresults,
      by = join_by(gene_name == gene_name),
      suffix = c("_dnam", "_txome")
    ) %>%
      drop_na() %>%
      mutate(
        outcome_txome = case_when(
          padj_txome < 0.05 & log2FoldChange > 0 ~ 1,
          padj_txome < 0.05 & log2FoldChange < 0 ~ -1,
          .default = 0
        ),
        outcome_dnam = case_when(
          padj_dnam < 0.05 & NES > 0 ~ 1,
          padj_dnam < 0.05 & NES < 0 ~ -1,
          .default = 0
        )
      ) %>%
      mutate(
        outcome_combined = case_when(
          outcome_txome != 0 & outcome_dnam == 0 ~ "txome",
          outcome_txome == 0 & outcome_dnam != 0 ~ "dnam",
          outcome_txome != 0 & outcome_dnam != 0 ~ "both",
          outcome_txome == 0 & outcome_dnam == 0 ~ "none"
        ) %>%
          factor(
            levels = c("txome", "dnam", "both", "none"),
            labels = c("DEG", "DMR", "Both", "ns")
          )
      )
  }
)

# Plot quadrant plot
# All DMPs (All 6 comparisons, 3 cols x 2 rows grid)
list_plots_all_dmps <- map2(
  results_merge,
  axis_titles,
  \(results, titles) {
    results %>%
      dplyr::filter(!outcome_combined %in% c("ns")) %>%
      mutate(
        symbol_repel = case_when(
          outcome_combined %in% c("DEG", "ns") ~ NA,
          .default = gene_name
        )
      ) %>%
      ggplot(
        aes(
          x = NES,
          y = log2FoldChange,
          fill = outcome_combined,
          label = symbol_repel
        )
      ) +
      geom_point(
        colour = "black",
        shape = 21,
        alpha = 0.7,
        stroke = 0.5,
        size = 2
      ) +
      geom_text_repel(
        aes(colour = outcome_combined),
        min.segment.length = 2,
        max.overlaps = 30,
        force = 5,
        force_pull = 0.3,
        size = 3,
        bg.color = "gray90",
        bg.r = 0.15
      ) +
      geom_quadrant_lines(linetype = "dotted") +
      scale_fill_manual(values = palette_quadrant) +
      scale_colour_manual(values = palette_quadrant) +
      scale_x_continuous(limits = symmetric_limits(c(-3, 3))) +
      scale_y_continuous(limits = symmetric_limits(c(-8, 8))) +
      theme_classic() +
      guides(shape = "none") +
      labs(
        x = str_c("NES ", titles[[1]]),
        y = str_c("log2FC ", titles[[2]]),
        fill = "Outcome",
        colour = "Outcome"
      )
  }
)

grid_all_dmps <- list_plots_all_dmps %>%
  map(\(plot) plot + theme(legend.position = "none")) %>%
  wrap_plots(
    ncol = 2,
    byrow = TRUE
  )

grid_all_dmps <- (grid_all_dmps | get_legend(list_plots_all_dmps[[1]])) +
  plot_layout(widths = c(10, 1))

# Zoom into NFKB and TGFB, vs. timepoint comparison only (3x2 grid)
results_merge_poi <- list(
  results_merge[[1]],
  results_merge[[2]],
  results_merge[[1]] %>%
    dplyr::filter(entrezid %in% geneids_nfkb),
  results_merge[[2]] %>%
    dplyr::filter(entrezid %in% geneids_nfkb),
  results_merge[[1]] %>%
    dplyr::filter(entrezid %in% geneids_tgfb),
  results_merge[[2]] %>%
    dplyr::filter(entrezid %in% geneids_tgfb)
)

list_plots_poi <- map2(
  results_merge_poi,
  axis_titles_poi,
  \(results, titles) {
    results %>%
      dplyr::filter(!outcome_combined %in% c("ns")) %>%
      mutate(
        symbol_repel = case_when(
          outcome_combined %in% c("DEG", "ns") ~ NA,
          .default = gene_name
        )
      ) %>%
      ggplot(
        aes(
          x = NES,
          y = log2FoldChange,
          fill = outcome_combined,
          label = symbol_repel
        )
      ) +
      geom_point(
        colour = "black",
        shape = 21,
        alpha = 0.7,
        stroke = 0.5,
        size = 2
      ) +
      geom_text_repel(
        aes(colour = outcome_combined),
        min.segment.length = 2,
        max.overlaps = 30,
        force = 5,
        force_pull = 0.3,
        size = 3,
        bg.color = "gray90",
        bg.r = 0.15
      ) +
      geom_quadrant_lines(linetype = "dotted") +
      scale_fill_manual(values = palette_quadrant) +
      scale_colour_manual(values = palette_quadrant) +
      scale_x_continuous(limits = symmetric_limits(c(-3, 3))) +
      scale_y_continuous(limits = symmetric_limits(c(-8, 8))) +
      theme_classic() +
      guides(shape = "none") +
      labs(
        x = str_c("NES ", titles[[1]]),
        y = str_c("log2FC ", titles[[2]]),
        fill = "Outcome",
        colour = "Outcome"
      )
  }
)

grid_poi <- list_plots_poi %>%
  map(\(plot) plot + theme(legend.position = "none")) %>%
  wrap_plots(
    ncol = 2,
    byrow = TRUE
  )

grid_poi <- (grid_poi | get_legend(list_plots_poi[[1]])) +
  plot_layout(widths = c(10, 1))


# Save data
## Plot - all DMPs
ggsave(
  filename = here::here(
    "output",
    "plots_quadrant",
    "quadrant_all_dmps.png"
  ),
  plot = grid_all_dmps,
  height = 10,
  width = 9,
  scale = 1
)

## Plot - POIs
ggsave(
  filename = here::here(
    "output",
    "plots_quadrant",
    "quadrant_poi.png"
  ),
  plot = grid_poi,
  height = 10,
  width = 9,
  scale = 1
)
