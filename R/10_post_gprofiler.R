# Load data

names_dmr <- readRDS(here(
  "output",
  "data_dmr",
  "names_dmr.RDS"
))

# Run g:OST GO enrichment

results_gost <- map(
  names_dmr,
  \(x) gost(x, organism = "hsapiens")
) %>%
  compact()

plots_gost <- map(
  results_gost,
  \(x) gostplot(x, interactive = FALSE)
)

# Save data

walk2(
  list(results_gost, plots_gost),
  c("results_gost", "plots_gost"),
  \(x, y) saveRDS(x, file = here("output", "data_gprofiler", str_c(y, ".RDS")))
)
