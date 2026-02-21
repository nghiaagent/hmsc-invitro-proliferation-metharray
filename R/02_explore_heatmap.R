here::i_am("R/02_explore_heatmap.R")

####################
# Create exploratory heatmap of normalised M values
# For my methylationEPICv1 data
####################

# Import packages
library(biganalytics)
library(circlize)
library(ComplexHeatmap)
library(here)
library(palettetown)
library(tidyverse)
library(viridis)

# Load data
quant_ratioset_funnorm_filter <- readRDS(
  file = file.path(
    "output",
    "quant_ratioset_funnorm_filter.RDS"
  )
)

targets <- readRDS(
  file = file.path(
    "output",
    "targets.RDS"
  )
)

# Plots normalised M values as a heatmap.
# Get normalised M value matrix
m_heatmap <- getM(quant_ratioset_funnorm_filter) %>%
  t() %>%
  scale() %>%
  t()

# Set color scheme and breaks
## For methylation M vals (z-scores)
col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

## For annotation data
col_cell_line <- pokepal(6)[c(3, 4)]
col_treatment <- pokepal(283)[c(11, 5, 7)]
col_timepoint <- pokepal(191)[c(8, 3)]

# Build annotation; include only necessary metadata
## Dataframe of annotation data
anno <- tibble(
  `Cell population` = colData(quant_ratioset_funnorm_filter)[["cell_line"]],
  `Timepoint` = colData(quant_ratioset_funnorm_filter)[["timepoint"]],
  `Treatment` = colData(quant_ratioset_funnorm_filter)[["treatment"]]
)

## Define colour mapping
anno_cols <- list(
  `Cell population` = c(
    "hMSC" = col_cell_line[1],
    "hNSC" = col_cell_line[2]
  ),
  `Timepoint` = c("early" = col_timepoint[1], "late" = col_timepoint[2]),
  `Treatment` = c(
    "untreated" = col_treatment[1],
    "heparin" = col_treatment[2],
    "neurosphere" = col_treatment[3]
  )
)

## Build ComplexHeatmap metadata annotation object
anno_object <- HeatmapAnnotation(
  df = anno,
  which = "col",
  col = anno_cols,
  annotation_height = 0.6,
  annotation_legend_param = list(
    `Cell population` = list(
      nrow = 2,
      title = "Cell population",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),
    `Timepoint` = list(
      nrow = 2,
      title = "Timepoint",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),
    `Treatment` = list(
      nrow = 3,
      title = "Treatment",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12, fontface = "bold")
    )
  )
)

# Select probes to be labelled
anno_genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = seq(1, nrow(m_heatmap), 60),
    labels = rownames(m_heatmap)[seq(1, nrow(m_heatmap), 60)],
    labels_gp = gpar(fontsize = 8, fontface = "bold"),
    padding = 0.75
  ),
  width = unit(2.0, "cm") +
    max_text_width(
      rownames(m_heatmap)[seq(1, nrow(m_heatmap), 60)],
      gp = gpar(
        fontsize = 8,
        fontface = "bold"
      )
    )
)

# Perform supervised clustering on probes
clusters_probes <- bigkmeans(m_heatmap, centers = 10)

# Create heatmap object
heatmap <- Heatmap(
  m_heatmap,
  name = "M-value\nZ-\nscore",
  col = colorRamp2(breaks, col),
  border = FALSE,

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(8, "cm"),
    legend_height = unit(5.0, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12, fontface = "bold")
  ),

  # row (gene) parameters
  show_row_dend = TRUE,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  row_split = clusters_probes$cluster,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,

  # column (sample) parameters
  cluster_column_slices = TRUE,
  column_split = colData(quant_ratioset_funnorm_filter)[["cell_line"]],
  cluster_columns = FALSE,
  show_column_dend = TRUE,
  show_column_names = FALSE,

  # specify top and bottom annotations
  top_annotation = anno_object
)

# Export heatmap of M values
png(
  file = "./output/plots_heatmap/heatmap_m_vals.png",
  width = 12,
  height = 12,
  units = "in",
  res = 600
)

plot(heatmap)

dev.off()
