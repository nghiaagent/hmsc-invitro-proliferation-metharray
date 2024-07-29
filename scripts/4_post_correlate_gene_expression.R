# Run ORA

# Load data

mcsea_heparin <- readRDS(file = file.path(
  "output",
  "data_dmr",
  "mcsea_deltabeta",
  "dmr_hMSC_heparin.RDS"
))

mcsea_timepoint <- readRDS(file = file.path(
  "output",
  "data_dmr",
  "mcsea_deltabeta",
  "dmr_hMSC_timepoint.RDS"
))
