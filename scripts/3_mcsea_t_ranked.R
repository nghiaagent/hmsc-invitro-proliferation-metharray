# Run mCSEA

# Load data

quant_beta_vals <- readRDS(file.path("output", "data_meth", "quant_beta_vals_hmsc.RDS"))

fit_contrasts <- readRDS(file.path("output", "data_meth", "fit_contrasts.RDS"))

t_timepoint <- readRDS(file.path("output", "data_meth", "t_timepoint.RDS"))

t_heparin <- readRDS(file.path("output", "data_meth", "t_heparin.RDS"))


# Run mCSEATest

mcsea_timepoint <- mCSEATest(
  rank = t_timepoint,
  methData = quant_beta_vals$E,
  pheno = quant_beta_vals$targets,
  minCpGs = 5,
  nproc = 16,
  platform = "EPIC"
)

mcsea_heparin <- mCSEATest(
  rank = t_heparin,
  methData = quant_beta_vals$E,
  pheno = quant_beta_vals$targets,
  minCpGs = 5,
  nproc = 16,
  platform = "EPIC"
)


# Get associations table

associations <- mcsea_timepoint[["promoters_association"]]

associations_excl <- associations[as_vector(map(associations, length)) < 5]

mCSEAPlot(
  mcsea_timepoint,
  regionType = "promoters",
  dmrName = "DCN",
  makePDF = FALSE
)

# Save data

saveRDS(mcsea_heparin,
        file = file.path("output", "data_dmr", "dmr_hMSC_heparin.RDS"))

saveRDS(mcsea_timepoint,
        file = file.path("output", "data_dmr", "dmr_hMSC_timepoint.RDS"))

openxlsx::write.xlsx(
  list(
    promoters = mcsea_heparin[["promoters"]],
    genes = mcsea_heparin[["promoters"]],
    CGI = mcsea_heparin[["CGI"]]
  ),
  file = file.path("output", "data_dmr", "dmr_hMSC_heparin.xlsx"),
  asTable = TRUE,
  rowNames = TRUE
)

openxlsx::write.xlsx(
  list(
    promoters = mcsea_timepoint[["promoters"]],
    genes = mcsea_timepoint[["promoters"]],
    CGI = mcsea_timepoint[["CGI"]]
  ),
  file = file.path("output", "data_dmr", "dmr_hMSC_timepoint.xlsx"),
  asTable = TRUE,
  rowNames = TRUE,
)