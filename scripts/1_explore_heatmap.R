# Load data

quant_ratioset_funnorm_filter <- readRDS(file = file.path(
    "output",
    "quant_ratioset_funnorm_filter.RDS"
))

targets <- readRDS(file = file.path(
    "output",
    "targets.RDS"
))

# Plots normalised M values as a heatmap.
# Get normalised M value matrix

m_heatmap <- t(scale(t(quant_ratioset_funnorm_filter@assays@data@listData[["M"]])))

# Set color scheme and breaks

## For methylation M vals (z-scores)

col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

## For annotation data

col_cell_line <- palettetown::pokepal(6)[c(3, 4)]

col_treatment <- palettetown::pokepal(283)[c(11, 5, 7)]

col_timepoint <- palettetown::pokepal(191)[c(8, 3)]

# Build annotation; include only necessary metadata

## Dataframe of annotation data

anno <- tibble(
    `Cell population` = quant_ratioset_funnorm_filter@colData@listData[["cell_line"]],
    `Timepoint` = quant_ratioset_funnorm_filter@colData@listData[["timepoint"]],
    `Treatment` = quant_ratioset_funnorm_filter@colData@listData[["treatment"]]
)

## Colour mapping

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

## ComplexHeatmap metadata annotation object

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

        max_text_width(rownames(m_heatmap)[seq(1, nrow(m_heatmap), 60)],
            gp = gpar(
                fontsize = 8, fontface = "bold"
            )
        )
)

# Cluster

clusters_probes <- bigkmeans(m_heatmap,
    centers = 10
)


# Create heatmap object

heatmap <- Heatmap(
    m_heatmap,
    name = "M-value\nZ-\nscore",
    col = colorRamp2(breaks, col),
    border = F,

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

    show_row_dend = T,
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_rot = 90,
    row_split = clusters_probes$cluster,
    show_row_names = F,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,

    # column (sample) parameters

    cluster_column_slices = T,
    column_split = quant_ratioset_funnorm_filter@colData@listData[["cell_line"]],
    cluster_columns = F,
    show_column_dend = T,
    show_column_names = F,

    # specify top and bottom annotations

    top_annotation = anno_object
)

# Export heatmap


# Export heatmap of raw data

png(
    file = "./output/plots_heatmap/heatmap_m_vals.png",
    width = 12,
    height = 12,
    units = "in",
    res = 600
)

export_heatmap <- plot(
    heatmap
)

dev.off()
