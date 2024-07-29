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
  # Uncomment below line if analysing hNSC too
  # mutate(cell_line = factor(cell_line, levels = c("hMSC", "hNSC"))) %>%
  mutate(timepoint = factor(timepoint, levels = c("early", "late"))) %>%
  # Uncomment below line if analysing hNSC too
  # mutate(passage = factor(passage, levels = c("p5", "p13"))) %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "heparin")))

# Create limma EList object containing M values

quant_m_vals <- new("EList")
quant_m_vals$E <- getM(quant_ratioset_funnorm_filter)
quant_m_vals$targets <- table_design
quant_m_vals$genes <- quant_ratioset_funnorm_filter@rowRanges

# Define design matrix for limma
## Treat time points and treatments (UT vs. Hep) as fixed effects

design <- model.matrix( ~ timepoint + treatment, data = table_design)

# Define contrasts

matrix_contrasts <- makeContrasts(P13vsP5 = timepointlate - 0,
                                  TvsUT = treatmentheparin - 0,
                                  levels = design)

rownames(matrix_contrasts)[1] <- colnames(design)[1]

# Run DMRcate to calculate per gene t-statistic

t_stat_timepoint <- cpg.annotate(
  datatype = "array",
  quant_m_vals$E,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  design = design,
  constrasts = TRUE,
  cont.matrix = matrix_contrasts,
  coef = 2
)

t_stat_heparin <- cpg.annotate(
  datatype = "array",
  quant_m_vals$E,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  design = design,
  constrasts = TRUE,
  cont.matrix = matrix_contrasts,
  coef = 3
)

dmr_timepoint <- dmrcate(t_stat_timepoint)

ranges_dmr_timepoint <- extractRanges(dmr_timepoint, genome = "hg19")

# Save data

saveRDS(quant_m_vals,
        file.path("output", "data_meth", "dmrcate", "quant_m_vals_hmsc.RDS"))

saveRDS(t_stat_timepoint,
        file.path("output", "data_meth", "dmrcate", "t_stat_timepoint.RDS"))

saveRDS(t_stat_heparin, 
        file.path("output", "data_meth", "dmrcate", "t_stat_heparin.RDS"))

saveRDS(dmr_timepoint,
        file.path("output", "data_meth", "dmrcate", "dmr_timepoint.RDS"))

saveRDS(ranges_dmr_timepoint,
        file.path("output", "data_meth", "dmrcate", "ranges_dmr_timepoint.RDS"))
