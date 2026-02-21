# Run preprocessing script if objects missing
if (
    !exists("quant_ratioset_funnorm") |
        !exists("targets") |
        !exists("quant_mset_none")
) {
    source(here::here("scripts", "0_preprocess.R"))
}

# Construct QC figure
# A: Per-sample detection p-value plots
# B: QC plot
# C: Density plot pre-norm
# D: Density plot post-norm

## Define features to keep
ntop <- 0.9

## Define hMSC samples
sel <- c(
    "200654430047_R01C01",
    "200654430047_R02C01",
    "200654430047_R03C01",
    "200654430047_R04C01",
    "200654430047_R05C01",
    "200654430047_R06C01"
)

# Construct EDA plots
# Top row: All samples Heatmap + PCA
# Middle row: MSC Heatmap + PCA
# Bottom row: Legend

## Construct Heatmap all samples

### Define params
#### Palette, number of steps

col <- viridis(n = 100)
breaks <- seq(0, 1, length.out = 100)

### Extract data, transpose so layout is probes by column and samples by row
### Get top variable genes
beta_heatmap <- getBeta(quant_ratioset_funnorm_filter)
beta_vars <- rowVars(beta_heatmap)
beta_heatmap_sel <-
    beta_heatmap[beta_vars >= quantile(beta_vars, ntop), ] %>%
    t()

### Cluster by fastcluster

# clusters_row <- fastcluster::hclust(dist(beta_heatmap_sel))
# clusters_col <- fastcluster::hclust(dist(t(beta_heatmap_sel)))

### Get shape mapping for sample

list_shape <-
    colData(quant_ratioset_funnorm_filter)[["condition_notreat"]] %>%
    case_match(
        "hMSCp5" ~ 16L,
        "hMSCp13" ~ 17L,
        "hNSCp6" ~ 15L,
        "hNSCp27" ~ 3L
    )

### Build annotation; include only necessary metadata

anno_object <-
    HeatmapAnnotation(
        `Cell type + Passage` = anno_simple(
            rep(1, times = 10),
            pch = list_shape,
            col = c("1" = pal2[2])
        ),
        `Treatment` = colData(quant_ratioset_funnorm_filter)[["treatment"]],
        col = list(
            "Treatment" = c(
                "untreated" = pal2[1],
                "heparin" = pal2[3],
                "neurosphere" = pal2[6]
            )
        ),
        which = "row",
        show_legend = c(FALSE, FALSE),
        show_annotation_name = FALSE
    )

### Plot

heatmap_all <- Heatmap(
    beta_heatmap_sel,
    col = colorRamp2(breaks, col),
    border = F,

    # parameters for the colour-bar that represents gradient of expression

    heatmap_legend_param = list(
        color_bar = "continuous",
        legend_direction = "horizontal",
        legend_width = unit(8, "cm"),
        legend_height = unit(5.0, "cm"),
        title = "Beta",
        title_position = "topleft",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),

    # row (gene) parameters

    show_row_dend = FALSE,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    # row_order = clusters_row$order,

    # column (sample) parameters

    show_column_dend = FALSE,
    show_column_names = FALSE,
    cluster_columns = TRUE,
    # column_order = clusters_col$order,

    # specify top and bottom annotations

    left_annotation = anno_object
)

### Export plot

png(
    file = "output/plots_QC/heatmap_all.png",
    width = 8,
    height = 4,
    units = "in",
    res = 300
)

draw(heatmap_all, heatmap_legend_side = "top")

dev.off()

## Construct Heatmap MSC only

### Extract data, transpose so layout is probes by column and samples by row
### Get top variable genes

beta_heatmap <- getBeta(quant_ratioset_funnorm_filter[, sel])
beta_vars <- rowVars(beta_heatmap)
beta_heatmap_sel <-
    beta_heatmap[beta_vars >= quantile(beta_vars, ntop), ] %>%
    t()

### Cluster by fastcluster

# clusters_row <- fastcluster::hclust(dist(beta_heatmap_sel))
# clusters_col <- fastcluster::hclust(dist(t(beta_heatmap_sel)))

### Get shape mapping for sample

list_shape <-
    colData(quant_ratioset_funnorm_filter[, sel])[["condition_notreat"]] %>%
    case_match(
        "hMSCp5" ~ 16L,
        "hMSCp13" ~ 17L,
        "hNSCp6" ~ 15L,
        "hNSCp27" ~ 3L
    )

### Build annotation; include only necessary metadata

anno_object <-
    HeatmapAnnotation(
        `Cell type + Passage` = anno_simple(
            rep(1, times = 6),
            pch = list_shape,
            col = c("1" = pal2[2])
        ),
        `Treatment` = colData(quant_ratioset_funnorm_filter[, sel])[[
            "treatment"
        ]],
        col = list(
            "Treatment" = c(
                "untreated" = pal2[1],
                "heparin" = pal2[3],
                "neurosphere" = pal2[6]
            )
        ),
        which = "row",
        show_legend = c(FALSE, FALSE),
        show_annotation_name = FALSE
    )

### Plot

heatmap_msc <- Heatmap(
    beta_heatmap_sel,
    col = colorRamp2(breaks, col),
    border = F,

    # parameters for the colour-bar that represents gradient of expression

    heatmap_legend_param = list(
        color_bar = "continuous",
        legend_direction = "horizontal",
        legend_width = unit(8, "cm"),
        legend_height = unit(5.0, "cm"),
        title = "Beta",
        title_position = "topleft",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),

    # row (gene) parameters

    show_row_dend = FALSE,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    # row_order = clusters_row$order,

    # column (sample) parameters

    show_column_dend = FALSE,
    show_column_names = FALSE,
    cluster_columns = TRUE,
    # column_order = clusters_col$order,

    # specify top and bottom annotations

    left_annotation = anno_object
)

### Export plot

png(
    file = "output/plots_QC/heatmap_hmsc.png",
    width = 8,
    height = 4,
    units = "in",
    res = 300
)

draw(heatmap_msc)

dev.off()

## Construct PCA plots

biplots <- map(
    .x = list(
        full = quant_ratioset_funnorm_filter,
        msc = quant_ratioset_funnorm_filter[, sel]
    ),
    .f = \(x) pca(getM(x), metadata = colData(x), removeVar = ntop)
) %>%
    map(
        .x = .,
        .f = \(x)
            biplot(
                x,
                x = "PC1",
                y = "PC2",
                lab = NULL,
                shape = "condition_notreat",
                shapeLegendTitle = "Cell type + Passage",
                colby = "treatment",
                colLegendTitle = "Treatment",
                colkey = c(
                    "untreated" = pal2[1],
                    "heparin" = pal2[3],
                    "neurosphere" = pal2[6]
                ),
                legendPosition = "bottom",
                labSize = 1,
                axisLabSize = 8
            )
    )

plots <-
    wrap_plots(
        fig(
            "output/plots_QC/heatmap_all.png",
            link_dim = TRUE,
            b_margin = ggplot2::margin(0, 0, 0, 0)
        ),
        biplots[[1]] +
            theme(
                legend.position = "none",
                plot.margin = margin(
                    t = 1,
                    r = 1,
                    b = 1,
                    l = 1
                )
            ),
        fig(
            "output/plots_QC/heatmap_hmsc.png",
            link_dim = TRUE,
            b_margin = ggplot2::margin(0, 0, 0, 0)
        ),
        biplots[[2]] +
            theme(
                legend.position = "none",
                plot.margin = margin(
                    t = 1,
                    r = 1,
                    b = 1,
                    l = 1
                )
            ),
        ggpubr::get_legend(biplots[[1]]),
        design = c(
            area(1, 1, 6, 4),
            area(1, 5, 6, 6),
            area(7, 1, 12, 4),
            area(7, 5, 12, 6),
            area(13, 1, 13, 6)
        )
    )

ggsave(
    filename = "output/plots_QC/EDA_combined.png",
    plot = plots,
    width = 12,
    height = 10,
    units = "in",
    scale = 0.8,
    dpi = 600
)
