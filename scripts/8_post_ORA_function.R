# Obtain gene sets

list_msigdb <- list(
    GOBP = "c5.go.bp.v2023.2.Hs.entrez.gmt",
    GOMF = "c5.go.mf.v2023.2.Hs.entrez.gmt",
    GOCC = "c5.go.cc.v2023.2.Hs.entrez.gmt",
    h = "h.all.v2023.2.Hs.entrez.gmt",
    c2 = "c2.cgp.v2023.2.Hs.entrez.gmt",
    c3 = "c3.all.v2023.2.Hs.entrez.gmt",
    reactome = "c2.cp.reactome.v2023.2.Hs.entrez.gmt",
    KEGG = "c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt"
)

msigdb <- map(
    list_msigdb,
    \(x) read.gmt(here(
        "input",
        "genesets",
        "msigdb_v2023.2.Hs_GMTs",
        x
    ))
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
            rename("Name" = "ID", "genes" = "geneID")
    }

# Define function to run ORA on list of genes extracted from Venn diagram
# Takes gene list as input
# Outputs folder in ./output/data_enrichment with all enrichment RDS objects
# And gene lists ready for Cytoscape

run_ORA <- function(list_genes, name_output) {
    # Make prerequisite folders under ./output

    name_output <- as.character(name_output)

    path_output <- here(
        "output",
        "data_enrichment",
        "ORA",
        name_output
    )

    message(str_c("Output ORA results to", path_output, sep = " "))

    if (!dir.exists(path_output)) {
        dir.create(path_output, recursive = TRUE)
    }

    # Run ORA on GO gene sets

    output <- map2(
        msigdb,
        names(msigdb),
        \(x, y) {
            message(str_c(
                "Running ",
                y,
                " enrichment for ",
                name_output
            ))
            enricher(
                list_genes,
                TERM2GENE = x,
                minGSSize = 15,
                maxGSSize = 1000
            ) %>%
                setReadable("org.Hs.eg.db", "ENTREZID")
        }
    )

    # Export RDS

    message(str_c("Saving RDS to", name_output, sep = " "))

    saveRDS(output, file = file.path(path_output, "ORA_results.RDS"))

    # Export table for Cytoscape

    message(str_c(
        "Saving table for Cytoscape/EnrichmentMap to",
        name_output,
        sep = " "
    ))

    return(output)
}
