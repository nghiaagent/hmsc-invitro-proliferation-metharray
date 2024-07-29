# Run relevant scripts beforehand
# 3_mcsea_deltabeta_ranked
# dge_extract_venn_genes

# List of gene sets to be used for GSEA
# GO (CC, BP, MF)
# KEGG
# ReactomePA
# MSigDB h, c2, c3

# Obtain gene sets

## GO

### GOBP
### Represented by MSigDB c5/GO/BP

msigdb_GOBP <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt")

### GOMF
### Represented by MSigDB c5/GO/MF

msigdb_GOMF <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt")

### GOCC
### Represented by MSigDB c5/GO/CC

msigdb_GOCC <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.cc.v2023.2.Hs.entrez.gmt")

## MSigDB

msigdb_h  <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt")

msigdb_c2 <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cgp.v2023.2.Hs.entrez.gmt")

msigdb_c3 <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/c3.all.v2023.2.Hs.entrez.gmt")

## ReactomePA
## Represented by MSigDB c2/CP/Reactome

msigdb_reactome <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.entrez.gmt")

## KEGG
## Represented by MSigDB c2/CP/KEGG_LEGACY

msigdb_KEGG <-
  read.gmt(
    "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt"
  )

# Define function to convert EnrichResult to something that EnrichmentMap accepts

convert_EnrichResult_to_EnrichmentMap_table <-
  function(EnrichResult) {
    EnrichResult@result %>%
      # # filter for term size to keep only gene count >= 5
      # filter(Count >= 5) %>%
      # format gene list column
      mutate(geneID = gsub("/", ",", .$geneID)) %>%
      # add column for phenotype
      mutate(phenotype = 1) %>%
      # Select needed columns for EnrichmentMap in Cytoscape, rename them.
      select(c(
        "ID",
        "Description",
        "pvalue",
        "qvalue",
        "phenotype",
        "geneID"
      )) %>%
      rename('Name' = 'ID', 'genes' = 'geneID')
  }

# Define function to run ORA on list of genes extracted from Venn diagram
# Takes gene list as input
# Outputs folder in ./output/data_enrichment with all enrichment RDS objects
# And gene lists ready for Cytoscape

run_ORA <- function(genes_list, name_output) {
  # Make prerequisite folders under ./output
  
  name_output <- as.character(name_output)
  
  path_output <- file.path(getwd(), 'output', 'data_enrichment', 'ORA', name_output)
  
  message(str_c("Output ORA results to", path_output, sep = " "))
  
  if (!dir.exists(path_output)) {
    dir.create(path_output, recursive = TRUE)
  }
  
  # Run ORA on GO gene sets
  
  message(str_c("Running GOMF enrichment for", name_output, sep = " "))
  
  ORA_GOMF <- enricher(
    genes_list,
    TERM2GENE = msigdb_GOMF,
    minGSSize = 25,
    maxGSSize = 500
  ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  message(str_c("Running GOBP enrichment for", name_output, sep = " "))
  
  ORA_GOBP <- enricher(
    genes_list,
    TERM2GENE = msigdb_GOBP,
    minGSSize = 25,
    maxGSSize = 500
  ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  message(str_c("Running GOCC enrichment for", name_output, sep = " "))
  
  ORA_GOCC <- enricher(
    genes_list,
    TERM2GENE = msigdb_GOCC,
    minGSSize = 25,
    maxGSSize = 500
  ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  # Run ORA on KEGG
  
  message(str_c("Running KEGG enrichment for", name_output, sep = " "))
  
  ORA_KEGG <- enricher(
    genes_list,
    TERM2GENE = msigdb_KEGG,
    minGSSize = 25,
    maxGSSize = 500
  ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  # Run ORA on ReactomePA
  
  message(str_c("Running Reactome enrichment for", name_output, sep = " "))
  
  ORA_Reactome <-
    enricher(
      genes_list,
      TERM2GENE = msigdb_reactome,
      minGSSize = 25,
      maxGSSize = 500
    ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  # Run ORA on MSigDB
  
  message(str_c("Running MSigDB h enrichment for", name_output, sep = " "))
  
  ORA_MSigDB_h <-
    enricher(
      genes_list,
      TERM2GENE = msigdb_h,
      minGSSize = 25,
      maxGSSize = 500
    ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  message(str_c("Running MSigDB c2 enrichment for", name_output, sep = " "))
  
  ORA_MSigDB_c2 <-
    enricher(
      genes_list,
      TERM2GENE = msigdb_c2,
      minGSSize = 25,
      maxGSSize = 500
    ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  message(str_c("Running MSigDB c3 enrichment for", name_output, sep = " "))
  
  ORA_MSigDB_c3 <-
    enricher(
      genes_list,
      TERM2GENE = msigdb_c3,
      minGSSize = 25,
      maxGSSize = 500
    ) %>%
    setReadable('org.Hs.eg.db', 'ENTREZID')
  
  # Export RDS
  
  message(str_c("Saving RDS to", name_output, sep = " "))
  
  output <- list(
    "GOBP" = ORA_GOBP,
    "GOCC" = ORA_GOCC,
    "GOMF" = ORA_GOMF,
    "KEGG" = ORA_KEGG,
    "MSigDB_h" = ORA_MSigDB_h,
    "MSigDB_c2" = ORA_MSigDB_c2,
    "MSigDB_c3" = ORA_MSigDB_c3,
    "Reactome" = ORA_Reactome,
    "path" = as.character(path_output)
  )
  
  saveRDS(output, file = file.path(path_output, "ORA_results.RDS"))
  
  # Export table for Cytoscape
  
  message(str_c(
    "Saving table for Cytoscape/EnrichmentMap to",
    name_output,
    sep = " "
  ))
  
  ## GO gene sets
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_GOMF),
    file.path(path_output, "ORA_GOMF_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_GOBP),
    file.path(path_output, "ORA_GOBP_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_GOCC),
    file.path(path_output, "ORA_GOCC_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  ## KEGG
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_KEGG),
    file.path(path_output, "ORA_KEGG_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  ## ReactomePA
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_Reactome),
    file.path(path_output, "ORA_Reactome_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  ## MSigDB
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_MSigDB_h),
    file.path(path_output, "ORA_MSigDB_h_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_MSigDB_c2),
    file.path(path_output, "ORA_MSigDB_c2_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  write.table(
    convert_EnrichResult_to_EnrichmentMap_table(ORA_MSigDB_c3),
    file.path(path_output, "ORA_MSigDB_c3_table.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  return(output)
  
}