# Construct QC figure
# A: Per-sample detection p-value plots
# B: QC plot
# C: Density plot pre-norm
# D: Density plot post-norm

## Define features to keep

ntop <- 0.95

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

png(file = "output/plots_QC/combined.png",
    width = 8,
    height = 7,
    units = "in",
    res = 300)

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
walk2(., names(.),
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

plot(NULL,
     xaxt = 'n',
     yaxt = 'n',
     bty = 'n',
     ylab = '',
     xlab = '',
     xlim = 0:1,
     ylim = 0:1)
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

beta_heatmap <- getBeta(quant_ratioset_funnorm_filter)

### Get top 10% variable genes

beta_vars <- rowVars(beta_heatmap)

beta_heatmap_sel <- beta_heatmap[beta_vars >= quantile(beta_vars, ntop), ] %>%
  t()

### Cluster by seriation

clusters_row <- fastcluster::hclust(dist(beta_heatmap_sel))
clusters_col <- fastcluster::hclust(dist(t(beta_heatmap_sel)))

### Build annotation; include only necessary metadata

#### Dataframe of annotation data

anno <- colData(quant_ratioset_funnorm_filter) %$%
  tibble(
  "Cell typePassage" = .[["condition_notreat"]],
  "Treatment" = .[["treatment"]]
)

## Colour mapping

anno_cols <- list(
  "Cell typePassage" = c(
    'hMSCp5'  = pal2[1],
    'hMSCp13' = pal2[2],
    'hNSCp6'  = pal2[3],
    'hNSCp27' = pal2[4]
  ),
  "Treatment" = c(
    'untreated'   = pal2[1],
    'heparin'     = pal2[3],
    'neurosphere' = pal2[6]
  )
)

## ComplexHeatmap metadata annotation object

anno_object <-
  HeatmapAnnotation(
    `Cell type + Passage` = colData(quant_ratioset_funnorm_filter)[["condition_notreat"]],
    Treatment = colData(quant_ratioset_funnorm_filter)[["treatment"]],
    col = list(
      "Cell type + Passage" = c(
        'hMSCp5'  = pal2[2],
        'hMSCp13' = pal2[2],
        'hNSCp6'  = pal2[2],
        'hNSCp27' = pal2[2]
      ),
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

Heatmap(
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

  show_row_dend = T,
  show_row_names = F,
  cluster_rows = FALSE,
  row_order = clusters_row$order,

  # column (sample) parameters

  show_column_dend = T,
  show_column_names = F,
  cluster_columns = FALSE,
  column_order = clusters_col$order,
  
  # specify top and bottom annotations

  left_annotation = anno_object
  
) %>% plot()

## Construct PCA plots

quant_pca <- quant_ratioset_funnorm_filter %$%
  pca(getM(.), metadata = colData(.), removeVar = ntop)

quant_pca_reduced <- quant_ratioset_funnorm_filter[, sel] %$%
  pca(getM(.), metadata = colData(.), removeVar = ntop)

biplots_all <- map2(
  .x = c("PC1"),
  .y = c("PC2"),
  .f = \(x, y) biplot(quant_pca,
                      x = x,
                      y = y,
                      lab = NULL,
                      shape = "condition_notreat",
                      shapeLegendTitle = "Cell type + Passage",
                      colby = "treatment",
                      colLegendTitle = "Treatment",
                      colkey = c("untreated" = pal2[1],
                                 "heparin" = pal2[3],
                                 "neurosphere" = pal2[6]),
                      legendPosition = 'bottom')
)

biplots_msc <- map2(
  .x = c("PC1"),
  .y = c("PC2"),
  .f = \(x, y) biplot(quant_pca_reduced,
                      x = x,
                      y = y,
                      lab = NULL,
                      shape = "condition_notreat",
                      shapeLegendTitle = "Cell type + Passage",
                      colby = "treatment",
                      colLegendTitle = "Treatment",
                      colkey = c("untreated" = pal2[1],
                                 "heparin" = pal2[3],
                                 "neurosphere" = pal2[6]),
                      legendPosition = 'bottom')
)

plots <- c(
#  [HEATMAP HERE],
  map(c(biplots_all, biplots_msc),
      \ (x) x %<>% + theme(legend.position="none")),
  list(ggpubr::get_legend(biplots_all[[1]]))) %>%
  wrap_plots(design = c(area(1, 1, 2, 4),
                        area(3, 1, 4, 2),
                        area(3, 3, 4, 4),
                        area(5, 1, 5, 4)))