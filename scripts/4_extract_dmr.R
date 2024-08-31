# Define function to extract genes (rownames) with padj < 0.05

extract_siggenes <- function(fgsea_results, padj = 0.05) {
  filter(fgsea_results, padj < 0.05)
}

# Load data
results_mcsea <- readRDS(here::here("output",
                                    "data_dmr",
                                    "results_mcsea.RDS"))

# Get list of DMRs

list_dmr <- map(results_mcsea,
                \ (results) {
                  map(c(promoters = "promoters",
                        genes = "genes"),
                      \ (region) extract_siggenes(results[[region]]))
                }) %>%
  unlist(recursive = FALSE)

names_dmr <- map(list_dmr,
                 \ (x) {
                   if (nrow(x) == 0) {
                     return(character(0))
                   } else {
                     rownames(x)
                   }
                 })

entrezid_dmr <- map(list_dmr,
                    \ (x) {
                      if (nrow(x) == 0) {
                        return(character(0))
                      } else {
                      rownames(x) %>%
                        AnnotationDbi::select(org.Hs.eg.db,
                                              .,
                                              columns = "ENTREZID",
                                              keytype = "ALIAS") %>%
                        .[, 2]}
                    })

# Save data

walk2(list(list_dmr, names_dmr, entrezid_dmr),
     c("list_dmr", "names_dmr", "entrezid_dmr"),
     \ (x, y) saveRDS(x, file = here("output", "data_dmr", str_c(y, ".RDS"))))

# Save xlsx

openxlsx::write.xlsx(list_dmr,
                     file = here("output",
                                 "data_dmr",
                                 "list_dmr.xlsx"))


map2(entrezid_dmr, names(entrezid_dmr),
     \ (x, y) {
       x %>%
         list %>%
         fwrite(file = here("output",
                            "data_dmr",
                            str_c("list_dmr_", y, ".txt")))
     })
