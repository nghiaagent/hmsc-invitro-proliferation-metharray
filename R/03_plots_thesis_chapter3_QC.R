here::i_am("R/03_plots_thesis_chapter3_QC.R")

####################
# Build EDA plots for thesis chapter 3
####################

# Import packages
library(here)
library(minfi)
library(tidyverse)

# Load data
## Sample sheet and detection p-values from preprocessing script
targets <- readRDS(
  file = here::here(
    "output",
    "targets.RDS"
  )
)

det_p <- readRDS(
  file = here::here(
    "output",
    "det_p.RDS"
  )
)

## Ratiosets
quant_rg <- readRDS(
  file = here::here(
    "output",
    "quant_rg.RDS"
  )
)

quant_mset_none <- readRDS(
  file = here::here(
    "output",
    "quant_mset_none.RDS"
  )
)

quant_ratioset_funnorm <- readRDS(
  file = here::here(
    "output",
    "quant_ratioset_funnorm.RDS"
  )
)

quant_ratioset_funnorm_filter <- readRDS(
  file = here::here(
    "output",
    "quant_ratioset_funnorm_filter.RDS"
  )
)

# Construct QC figure
# A: Per-sample detection p-value plots
# B: QC plot
# C: Density plot pre-norm
# D: Density plot post-norm
## Define plot layout and output
png(
  file = "output/plots_QC/QC_combined.png",
  width = 10,
  height = 7,
  units = "in",
  res = 300
)

layout(
  matrix(
    c(1, 2, 3, 4, 5, 6, 7, 7),
    nrow = 4,
    ncol = 2,
    byrow = TRUE
  ),
  height = c(2, 2, 2, 1)
)

par(mar = c(2, 4, 2, 4))

## Plot mean detection p-values across all samples to identify any failed samples
order <- targets$sample_name
order2 <- targets$condition_notreat

barplot(
  colMeans(det_p),
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
list(
  "Raw beta" = quant_rg,
  "Funnorm beta" = getBeta(quant_ratioset_funnorm)
) %>%
  iwalk(\(x, y) {
    densityPlot(
      x,
      sampGroups = order,
      main = y,
      legend = FALSE,
      pal = pal2[order],
      lwd = 6
    )
  })

list(
  "Raw beta (color by cell type)" = quant_rg,
  "Funnorm beta (color by cell type)" = getBeta(quant_ratioset_funnorm)
) %>%
  iwalk(\(x, y) {
    densityPlot(
      x,
      sampGroups = order2,
      main = y,
      legend = FALSE,
      pal = pal2[11:14],
      lwd = 6
    )
    legend("top", legend = levels(order2), text.col = pal2[11:14])
  })

# Plot legend
plot(
  NULL,
  xaxt = "n",
  yaxt = "n",
  bty = "n",
  ylab = "",
  xlab = "",
  xlim = 0:1,
  ylim = 0:1
)
legend("bottom", legend = levels(order), fill = pal2, ncol = 5)

dev.off()
