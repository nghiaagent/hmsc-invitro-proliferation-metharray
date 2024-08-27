# Load data

quant_ratioset_funnorm_filter <- readRDS(file = file.path("output", "quant_ratioset_funnorm_filter.RDS"))
quant_m_vals <- readRDS(file = file.path("output", "data_meth", "limma", "quant_m_vals_hmsc.RDS"))
beta_sel <- quant_ratioset_funnorm_filter[, c("200654430047_R01C01",
                                               "200654430047_R02C01",
                                               "200654430047_R03C01",
                                               "200654430047_R04C01")] %>%
  getBeta()

# Extract beta values
# Calculate differences in beta values between P5 UT and P13 UT
# hMSC early UT: 200654430047_R01C01
# hMSC late UT: 200654430047_R03C01

beta_early_untreated <- quant_ratioset_funnorm_filter[, "200654430047_R01C01"] %>%
  getBeta() %>%
  .[, 1]

beta_early_treated   <- quant_ratioset_funnorm_filter[, "200654430047_R02C01"] %>%
  getBeta() %>%
  .[, 1]

beta_late_untreated  <- quant_ratioset_funnorm_filter[, "200654430047_R03C01"] %>%
  getBeta() %>%
  .[, 1]

deltabeta_timepoint <- beta_late_untreated - beta_early_untreated
deltabeta_treatment <- beta_early_treated - beta_early_untreated

# Run mCSEATest

mcsea_timepoint <- mCSEATest(
  rank = deltabeta_timepoint,
  methData = beta_sel,
  pheno = quant_m_vals$targets,
  minCpGs = 5,
  nproc = 16,
  platform = "EPIC"
)

mcsea_heparin <- mCSEATest(
  rank = deltabeta_treatment,
  methData = beta_sel,
  pheno = quant_m_vals$targets,
  minCpGs = 5,
  nproc = 16,
  platform = "EPIC"
)

# Save data

saveRDS(mcsea_heparin,
        file = file.path("output",
                         "data_dmr",
                         "mcsea_deltabeta",
                         "dmr_hMSC_heparin.RDS"))

saveRDS(
  mcsea_timepoint,
  file = file.path(
    "output",
    "data_dmr",
    "mcsea_deltabeta",
    "dmr_hMSC_timepoint.RDS"
  )
)

openxlsx::write.xlsx(
  list(
    promoters = mcsea_heparin[["promoters"]],
    genes = mcsea_heparin[["genes"]],
    CGI = mcsea_heparin[["CGI"]]
  ),
  file = file.path("output",
                   "data_dmr",
                   "mcsea_deltabeta",
                   "dmr_hMSC_heparin.xlsx"),
  asTable = TRUE,
  rowNames = TRUE
)

openxlsx::write.xlsx(
  list(
    promoters = mcsea_timepoint[["promoters"]],
    genes = mcsea_timepoint[["genes"]],
    CGI = mcsea_timepoint[["CGI"]]
  ),
  file = file.path(
    "output",
    "data_dmr",
    "mcsea_deltabeta",
    "dmr_hMSC_timepoint.xlsx"
  ),
  asTable = TRUE,
  rowNames = TRUE,
)