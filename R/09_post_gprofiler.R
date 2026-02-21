here::i_am("R/09_post_gprofiler.R")

####################
# Run GO enrichment using g:gOST
####################

# Import packages
library(gprofiler2)
library(purrr)
library(tidyverse)

# Load data
names_dmr <- readRDS(
  file = here::here(
    "output",
    "data_dmr",
    "names_dmr.RDS"
  )
)

# Run g:OST GO enrichment
results_gost <- names_dmr %>%
  map(\(x) gost(x, organism = "hsapiens")) %>%
  compact()

plots_gost <- results_gost %>%
  map(\(x) gostplot(x, interactive = FALSE))

# Save data
walk2(
  list(results_gost, plots_gost),
  c("results_gost", "plots_gost"),
  \(x, y) saveRDS(x, file = here("output", "data_gprofiler", str_c(y, ".RDS")))
)
