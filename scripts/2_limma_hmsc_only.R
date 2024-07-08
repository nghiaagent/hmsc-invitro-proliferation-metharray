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

# Create limma EList object containing beta values

quant_beta_vals <- new("EList")
quant_beta_vals$E <- getM(quant_ratioset_funnorm_filter)
quant_beta_vals$targets <- table_design
quant_beta_vals$genes <- quant_ratioset_funnorm_filter@rowRanges

# Define design matrix for limma
## Treat time points and treatments (UT vs. Hep) as fixed effects

design <- model.matrix( ~ timepoint + treatment, data = table_design)

colnames(design) %<>% make.names(.)

# Define contrasts

matrix_contrasts <- makeContrasts(P13vsP5 = timepointlate - 0,
                                  TvsUT = treatmentheparin - 0,
                                  levels = design)

# Apply limma model fit

fit_contrasts <- lmFit(quant_beta_vals, design) %>%
  eBayes() %>%
  contrasts.fit(contrasts = matrix_contrasts) %>%
  eBayes()

# Get t-statistics for mCSEA downstream

t_timepoint <- topTable(fit_contrasts,
                        coef = 1,
                        number = Inf,
                        sort = "none")$t
names(t_timepoint) <- rownames(topTable(
  fit_contrasts,
  coef = 1,
  number = Inf,
  sort = "none"
))

t_heparin <- topTable(fit_contrasts,
                      coef = 2,
                      number = Inf,
                      sort = "none")$t
names(t_heparin) <- rownames(topTable(
  fit_contrasts,
  coef = 2,
  number = Inf,
  sort = "none"
))

# Save data

saveRDS(quant_beta_vals,
        file.path("output", "data_meth", "quant_beta_vals_hmsc.RDS"))

saveRDS(fit_contrasts,
        file.path("output", "data_meth", "fit_contrasts.RDS"))

saveRDS(t_timepoint,
        file.path("output", "data_meth", "t_timepoint.RDS"))

saveRDS(t_heparin, 
        file.path("output", "data_meth", "t_heparin.RDS"))
