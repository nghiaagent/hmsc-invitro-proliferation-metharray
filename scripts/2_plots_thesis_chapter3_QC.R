# Construct QC figure
# A: Per-sample detection p-value plots
# B: QC plot
# C: Density plot pre-norm
# D: Density plot post-norm


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