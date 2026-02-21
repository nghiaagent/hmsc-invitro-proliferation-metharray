here::i_am("R/99_pipeline_03_postdmr_plots.R")

########################
# Run pre-processing, QC, normalisation
# And plotting prior to DMR analysis
########################

# Run scripts of pipeline
c(
  "07_post_manhattan.R",
  "08_post_correlate_gene_expression.R",
  "09_post_gprofiler.R"
) %>%
  walk(
    \(x) {
      message(paste0("Sourcing ", here::here("scripts", x)))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
