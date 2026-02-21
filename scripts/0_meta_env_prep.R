# Install pacman, a package manager to install and load multiple packages at once
if (!requireNamespace("pacman")) {
  install.packages("pacman")
}
library(pacman)

# Install Bioconductor package manager
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

# Load all CRAN packages
p_load(
  renv,
  colorspace,
  RColorBrewer,
  tidyverse,
  ggrepel,
  ggpubr,
  ggalt,
  ggplotify,
  GSEABase,
  cowplot,
  gridExtra,
  statmod,
  writexl,
  ggpmisc,
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
  here,
  matrixStats,
  biganalytics,
  openxlsx,
  magrittr,
  ggVennDiagram,
  vctrs,
  furrr,
  here,
  patchwork,
  fastcluster,
  gprofiler2,
  CMplot,
  PupillometryR
)

p_load_gh(
  "BradyAJohnston/figpatch"
)

# List Bioconductor packages to load/install
list_Bioc_Pkg <- c(
  "EnhancedVolcano",
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
  "mCSEA",
  "minfi",
  "Gviz",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "BiocParallel",
  "shinyMethyl",
  "edgeR",
  "clusterProfiler",
  "PCAtools",
  "trackViewer",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "biomaRt",
  "EnsDb.Hsapiens.v86",
  "ggmanh"
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


# Clean up package list

rm(list_Bioc_Pkg)

# Limit number of cores used due to memory issues on laptops.

BiocParallel::register(SnowParam(workers = 8), default = TRUE)
bpparam()

# Set up multicore processing with furrr::map
future::plan(strategy = multisession, workers = 8)

# TO ADD: Create folders needed for outputting files
# Add color palettes
pal2 <- palette.colors(n = 36, palette = "Polychrome 36")

# Source other scripts
source(
  here::here(
    "scripts",
    "0_define_GOIs.R"
  )
)

source(
  here::here(
    "scripts",
    "0_define_colours.R"
  )
)

# Set up biomaRt
mart <- useMart("ensembl", "hsapiens_gene_ensembl")

t2g <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "ensembl_gene_id_version",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  mart = mart
)
