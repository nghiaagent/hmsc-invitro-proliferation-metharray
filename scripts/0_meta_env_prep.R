# Install pacman, a package manager to install and load multiple packages at once

if (!requireNamespace("pacman")) {
  install.packages("pacman")
}
library(pacman)

# Install Bioconductor package manager

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

# List Bioconductor packages to load/install

list_Bioc_Pkg <- c(
  "BiocHubsShiny",
  "AnnotationHub",
  "Homo.sapiens",
  "limma",
  "Glimma",
  "ComplexHeatmap",
  "SummarizedExperiment",
  "ensembldb",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "BiocStyle",
  "mCSEA",
  "missMethyl",
  "minfi",
  "Gviz",
  "DMRcate",
  "IlluminaHumanMethylationEPICmanifest",
  "BiocParallel",
  "shinyMethyl",
  "edgeR",
  "clusterProfiler"
)

# Install Bioconductor packages, if they are not yet installed
# Swap update = TRUE / FALSE depending on need to update pkg

BiocManager::install(
  pkgs = list_Bioc_Pkg,
  update = FALSE
)

# Load Bioconductor packages

invisible(lapply(list_Bioc_Pkg, function(x) {
  library(x, character.only = TRUE)
}))

# Load all CRAN packages

p_load(
  RColorBrewer,
  tidyverse,
  ggrepel,
  ggpubr,
  ggalt,
  ggplotify,
  cowplot,
  gridExtra,
  statmod,
  writexl,
  cowplot,
  viridis,
  circlize,
  magick,
  cluster,
  rmdformats,
  future,
  shinybusy,
  data.table,
  DT,
  matrixStats,
  biganalytics,
  openxlsx,
  magrittr,
  ggVennDiagram,
  vctrs,
  furrr,
  here
)

# Clean up package list

rm(list_Bioc_Pkg)

# Limit number of cores used due to memory issues on laptops.

BiocParallel::register(SnowParam(workers = 8),
                       default = TRUE)
bpparam()

# Set up multicore processing with furrr::map

future::plan(strategy = multisession, workers = 8)

# TO ADD: Create folders needed for outputting files

# Add color palettes

pal2 <- palette.colors(n = 10, palette = "Polychrome 36")
