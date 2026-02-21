# Run ORA

# Load data

mcsea_heparin <- readRDS(
  file = file.path(
    "output",
    "data_dmr",
    "mcsea_deltabeta",
    "dmr_hMSC_heparin.RDS"
  )
)

mcsea_timepoint <- readRDS(
  file = file.path(
    "output",
    "data_dmr",
    "mcsea_deltabeta",
    "dmr_hMSC_timepoint.RDS"
  )
)

# Get list of genes
# ENTREZID required for the ORA function

# Define lists of DMRs to be used

list_dmr <- list(
  heparin_promoters = mcsea_heparin$promoters,
  heparin_genes = mcsea_heparin$genes,
  timepoint_promoters = mcsea_timepoint$promoters,
  timepoint_genes = mcsea_timepoint$genes
)

# Extract ENTREZIDs

list_entrez <- purrr::map(list_dmr, .f = \(list) {
  filter(list, padj < 0.05) %>%
    rownames() %>%
    AnnotationDbi::select(org.Hs.eg.db, keys = ., "ENTREZID", "SYMBOL") %>%
    .[["ENTREZID"]]
}) %>%
  list_drop_empty()

# ORA between passages

ora_results <- map2(.x = list_entrez, .y = names(list_entrez), .f = \(x, y) {
  run_ORA(x, y)
})
