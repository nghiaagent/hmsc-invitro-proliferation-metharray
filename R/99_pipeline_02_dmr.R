here::i_am("R/99_pipeline_02_dmr.R")

########################
# Run pre-processing, QC, normalisation
# And plotting prior to DMR analysis
########################

# Run scripts of pipeline
c(
  "04_mcsea_deltabeta_ranked.R",
  "05_extract_dmr.R",
  "06_get_leading_edge_mean.R"
) %>%
  walk(
    \(x) {
      message(paste0("Sourcing ", here::here("scripts", x)))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
