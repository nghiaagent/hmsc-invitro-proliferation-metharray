source("renv/activate.R")

# Source local function files
c(
  "00_define_colours.R",
  "00_define_GOIs.R",
  "00_post_ORA_function.R"
) %>%
  purrr::walk(\(x) {
    message(paste0("Sourcing ", here::here("scripts", x)))
    source(here::here("R", x), echo = FALSE, verbose = FALSE)
  })
