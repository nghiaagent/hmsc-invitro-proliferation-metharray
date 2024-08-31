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

## Define plot layout and output

png(
  file = "output/plots_QC/QC_combined.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300
)

layout(matrix(
  c(1, 2,
    3, 4,
    5, 5),
  nrow = 3,
  ncol = 2,
  byrow = TRUE
), height = c(2, 2, 1))

par(mar = c(2, 4, 2, 4))

## Plot mean detection p-values across all samples to identify any failed samples

order <- targets$sample_name

barplot(
  colMeans(detP),
  col = pal2[order],
  las = 2,
  cex.names = 0.8,
  ylim = c(0, 0.0006),
  main = "Mean detection p-values",
  xaxt = "n"
)

## Plot mean channel intensity

quant_mset_none %>%
  getQC() %>%
  plotQC()

## Plot beta values density, before and after normalization

list("Raw beta" = quant_rg,
     "Funnorm beta" = getBeta(quant_ratioset_funnorm)) %$%
  walk2(
    .,
    names(.),
    \ (x, y) densityPlot(
      x,
      sampGroups = order,
      main = y,
      legend = FALSE,
      pal = pal2[order],
      lwd = 6
    )
  )

# Plot legend

plot(
  NULL,
  xaxt = 'n',
  yaxt = 'n',
  bty = 'n',
  ylab = '',
  xlab = '',
  xlim = 0:1,
  ylim = 0:1
)
legend("bottom",
       legend = levels(order),
       fill = pal2,
       ncol = 3)

dev.off()

# Construct EDA plots
# Top row: Heatmap
# Middle row: PC1-2; full dataset then hMSC only
# Bottom row: Legend

## Construct Heatmap

### Define params
#### Palette, number of steps

col <- viridis(n = 100)
breaks <- seq(0, 1, length.out = 100)

### Extract data, transpose so layout is probes by column and samples by row
### Get top variable genes

beta_heatmap <- getBeta(quant_ratioset_funnorm_filter)
beta_vars <- rowVars(beta_heatmap)
beta_heatmap_sel <-
  beta_heatmap[beta_vars >= quantile(beta_vars, ntop),] %>%
  t()

### Cluster by fastcluster

clusters_row <- fastcluster::hclust(dist(beta_heatmap_sel))
clusters_col <- fastcluster::hclust(dist(t(beta_heatmap_sel)))

### Get shape mapping for sample

list_shape <-
  colData(quant_ratioset_funnorm_filter)[["condition_notreat"]] %>%
  case_match("hMSCp5" ~ 16L,
             "hMSCp13" ~ 17L,
             "hNSCp6" ~ 15L,
             "hNSCp27" ~ 3L)

### Build annotation; include only necessary metadata

anno_object <-
  HeatmapAnnotation(
    `Cell type + Passage` = anno_simple(rep(1, times = 10),
                                        pch = list_shape,
                                        col = c("1" = pal2[2])),
    `Treatment` = colData(quant_ratioset_funnorm_filter)[["treatment"]],
    col = list(
      "Treatment" = c(
        'untreated'   = pal2[1],
        'heparin'     = pal2[3],
        'neurosphere' = pal2[6]
      )
    ),
    which = "row",
    show_legend = c(FALSE, FALSE),
    show_annotation_name = FALSE
  )

### Plot

heatmap <- Heatmap(
  beta_heatmap_sel,
  col = colorRamp2(breaks, col),
  border = F,
  
  # parameters for the colour-bar that represents gradient of expression
  
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'vertical',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title = "Beta",
    title_position = 'topleft',
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  
  # row (gene) parameters
  
  show_row_dend = FALSE,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  row_order = clusters_row$order,
  
  # column (sample) parameters
  
  show_column_dend = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_order = clusters_col$order,
  
  # specify top and bottom annotations
  
  left_annotation = anno_object
)

### Export plot

png(file = "output/plots_QC/heatmap.png", width = 10, height = 4, units = "in", res = 300)

draw(heatmap)

dev.off()

## Construct PCA plots

biplots <- map(
  .x = list(full = quant_ratioset_funnorm_filter,
            msc = quant_ratioset_funnorm_filter[, sel]),
  .f = \ (x) pca(getM(x), metadata = colData(x), removeVar = ntop)
) %>%
  map(
    .x = .,
    .f = \(x) biplot(
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
      legendPosition = 'bottom'
   )
  )

plots <-
  wrap_plots(
    fig("output/plots_QC/heatmap.png",
        link_dim = TRUE,
        b_margin = ggplot2::margin(0, 0, 0, 0)
),
    biplots[[1]] + theme(legend.position = "none"),
    biplots[[2]] + theme(legend.position = "none"),
    ggpubr::get_legend(biplots[[1]]),
    design = c(area(1, 1, 7, 4),
               area(8, 1, 12, 2),
               area(8, 3, 12, 4),
               area(13, 1, 13, 4))
  )

ggsave(filename = "output/plots_QC/EDA_combined.png",
       plot = plots,
       width = 12,
       height = 12,
       units = "in",
       dpi = 300)
