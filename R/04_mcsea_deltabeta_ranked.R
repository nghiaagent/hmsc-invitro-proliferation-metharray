here::i_am("R/04_mcsea_deltabeta_ranked.R")

####################
# Calculate scores for mCSEA using rank of beta differences
####################

# Import packages
library(here)
library(mCSEA)
library(minfi)
library(tidyverse)

# Load data
quant_ratioset_funnorm_filter <- readRDS(
  file = file.path(
    "output",
    "quant_ratioset_funnorm_filter.RDS"
  )
)
## Filter for hMSC, untreated vs treated samples only
quant_ratioset_funnorm_filter <- quant_ratioset_funnorm_filter[,
  colData(quant_ratioset_funnorm_filter)$cell_line == "hMSC"
]
quant_ratioset_funnorm_filter <- quant_ratioset_funnorm_filter[,
  colData(quant_ratioset_funnorm_filter)$treatment != "neurosphere"
]

# Make design matrix
table_design <- colData(quant_ratioset_funnorm_filter) %>%
  as.data.frame() %>%
  mutate(
    sample_name = factor(
      sample_name,
      levels = c(
        "hMSC_p5_untreated",
        "hMSC_p5_heparin",
        "hMSC_p13_untreated",
        "hMSC_p13_heparin"
      )
    )
  ) %>%
  mutate(
    slide = slide %>%
      factor(),
    timepoint = timepoint %>%
      factor(levels = c("early", "late")),
    treatment = treatment %>%
      factor(levels = c("untreated", "heparin"))
  )

# Create limma EList object containing beta values
quant_beta_vals <- new("EList")
quant_beta_vals$E <- getBeta(quant_ratioset_funnorm_filter)
quant_beta_vals$targets <- table_design
quant_beta_vals$genes <- rowRanges(quant_ratioset_funnorm_filter)

# Extract beta values and relevant ranks
beta <- list(
  early_untreated = "200654430047_R01C01",
  early_treated = "200654430047_R02C01",
  late_untreated = "200654430047_R03C01",
  late_treated = "200654430047_R04C01"
) %>%
  map(\(x) quant_beta_vals$E[, x])

# Run mCSEATest
results_mcsea <- list(
  timepoint = beta[["late_untreated"]] - beta[["early_untreated"]],
  treatment_early = beta[["early_treated"]] - beta[["early_untreated"]],
  treatment_late = beta[["late_treated"]] - beta[["late_untreated"]]
) %>%
  map(\(x) {
    mCSEATest(
      rank = x,
      methData = quant_beta_vals$E,
      pheno = quant_beta_vals$targets,
      minCpGs = 5,
      nproc = 1,
      platform = "EPIC"
    )
  })

# Save data
## RDS
saveRDS(
  results_mcsea,
  file = here::here(
    "output",
    "data_dmr",
    "results_mcsea.RDS"
  )
)

## Save xlsx
map2(
  .x = results_mcsea,
  .y = names(results_mcsea),
  .f = \(x, y) {
    openxlsx::write.xlsx(
      list(
        promoters = x[["promoters"]],
        genes = x[["genes"]],
        CGI = x[["CGI"]]
      ),
      file = file.path(
        "output",
        "data_dmr",
        str_c("results_mcsea_", y, ".xlsx")
      ),
      asTable = TRUE,
      rowNames = TRUE
    )
  }
)
