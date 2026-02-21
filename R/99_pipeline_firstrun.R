here::i_am("R/99_pipeline_firstrun.R")

########################
# Run entire pipeline
########################

# Run scripts of pipeline
c(
  "99_pipeline_01_predmr.R",
  "99_pipeline_02_dmr.R",
  "99_pipeline_03_postdmr_plots.R"
) %>%
  walk(
    \(x) {
      message(paste0("Sourcing ", here::here("scripts", x)))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
