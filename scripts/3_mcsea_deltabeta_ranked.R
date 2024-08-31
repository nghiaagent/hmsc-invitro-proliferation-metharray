# Load data

quant_ratioset_funnorm_filter <- readRDS(file = file.path("output", "quant_ratioset_funnorm_filter.RDS")) %>%
  .[, .@colData$cell_line == "hMSC"] %>%
  .[, .@colData$treatment != "neurosphere"]

# Make design matrix

table_design <- quant_ratioset_funnorm_filter@colData %>%
  as.data.frame() %>%
  mutate(sample_name = factor(
    sample_name,
    levels = c(
      "hMSC_p5_untreated",
      "hMSC_p5_heparin",
      "hMSC_p13_untreated",
      "hMSC_p13_heparin"
    )
  )) %>%
  mutate(slide = factor(slide)) %>%
  mutate(timepoint = factor(timepoint, levels = c("early", "late"))) %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "heparin")))

# Create limma EList object containing beta values

quant_beta_vals <- new("EList")
quant_beta_vals$E <- getBeta(quant_ratioset_funnorm_filter)
quant_beta_vals$targets <- table_design
quant_beta_vals$genes <- quant_ratioset_funnorm_filter@rowRanges

# Extract beta values and relevant ranks

beta <- map(list(early_untreated = "200654430047_R01C01",
                 early_treated = "200654430047_R02C01",
                 late_untreated = "200654430047_R03C01"),
            \ (x) quant_beta_vals$E[, x])

deltabeta <-
  list(timepoint = beta[["late_untreated"]] - beta[["early_untreated"]],
       treatment = beta[["early_treated"]]  - beta[["early_untreated"]])

# Run mCSEATest

results_mcsea <- map(
  deltabeta,
  \ (x) mCSEATest(
    rank = x,
    methData = quant_beta_vals$E,
    pheno = quant_beta_vals$targets,
    minCpGs = 5,
    nproc = 16,
    platform = "EPIC"
  )
)

# Save data

## RDS

saveRDS(results_mcsea,
        file = file.path("output",
                         "data_dmr",
                         "results_mcsea.RDS"))

## Save xlsx

map2(.x = results_mcsea,
     .y = names(results_mcsea),
     .f = \ (x, y) {
       openxlsx::write.xlsx(
         list(
           promoters = x[["promoters"]],
           genes = x[["genes"]],
           CGI = x[["CGI"]]
         ),
         file = file.path("output",
                          "data_dmr",
                          str_c("restults_mcsea_", y, ".xlsx")),
         asTable = TRUE,
         rowNames = TRUE
       )
     })
