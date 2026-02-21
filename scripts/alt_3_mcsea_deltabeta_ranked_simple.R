# Load data (replace with relevant path for your minfi RatioSet)
quant_ratioset <- readRDS(
    file = here::here(
        "quant_ratioset_funnorm_filter.RDS"
    )
)

# Create limma EList object containing beta values for mCSEA
quant_beta_vals <- new("EList")
quant_beta_vals$E <- getBeta(quant_ratioset)
quant_beta_vals$targets <- as.data.frame(quant_ratioset@colData)
quant_beta_vals$genes <- quant_ratioset@rowRanges

# Extract beta values and calculate change in beta (select your columns based on project)
beta_late_untreated <- quant_beta_vals$E[, "200654430047_R03C01"]
beta_early_untreated <- quant_beta_vals$E[, "200654430047_R01C01"]

deltabeta <- beta_late_untreated - beta_early_untreated

# Run mCSEATest (change nproc depending on computer, change platform depending on chip in use)
results_mcsea <- mCSEATest(
    rank = deltabeta,
    methData = quant_beta_vals$E,
    pheno = quant_beta_vals$targets,
    minCpGs = 5,
    nproc = 16,
    platform = "EPIC"
)

# Save data
# ## RDS
# saveRDS(
#     results_mcsea,
#     file = here::here(
#         "results_mcsea.RDS"
#     )
# )

# ## Save xlsx
# openxlsx::write.xlsx(
#     list(
#         promoters = results_mcsea[["promoters"]],
#         genes = results_mcsea[["genes"]],
#         CGI = results_mcsea[["CGI"]]
#     ),
#     file = here::here(
#         str_c("results_mcsea", ".xlsx")
#     ),
#     asTable = TRUE,
#     rowNames = TRUE
# )
