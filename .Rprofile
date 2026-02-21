source("renv/activate.R")

# Source local function files
library(tidyverse)

c(
  "00_define_colours.R",
  "00_define_GOIs.R"
) %>%
  purrr::walk(\(x) {
    message(paste0("Sourcing ", here::here("scripts", x)))
    source(here::here("R", x), echo = FALSE, verbose = FALSE)
  })
