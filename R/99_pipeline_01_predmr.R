here::i_am("R/99_pipeline_01_predmr.R")

########################
# Run pre-processing, QC, normalisation
# And plotting prior to DMR analysis
########################

# Run scripts of pipeline
c(
  "01_preprocess.R",
  "02_explore_heatmap.R",
  "03_plots_chapter6_explore.R",
  "03_plots_thesis_chapter3_EDA.R",
  "03_plots_thesis_chapter3_QC.R"
) %>%
  walk(
    \(x) {
      message(paste0("Sourcing ", here::here("scripts", x)))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
