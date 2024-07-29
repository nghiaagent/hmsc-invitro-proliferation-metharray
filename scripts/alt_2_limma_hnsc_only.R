# Load data

quant_ratioset_funnorm_filter <- readRDS(file = file.path("output",
                                                        "quant_ratioset_funnorm_filter.RDS"))

quant_ratioset_funnorm_filter <- quant_ratioset_funnorm_filter[, quant_ratioset_funnorm_filter@colData$cell_line == "hNSC"]

# Make design matrix

table_design <- quant_ratioset_funnorm_filter@colData %>%
  as.data.frame() %>%
  mutate(sample_name = factor(
    sample_name,
    levels = c(
      "hNSC_p6_untreated",
      "hNSC_p6_heparin",
      "hNSC_p27_untreated",
      "hNSC_p27_heparin"
    )
  )) %>%
  mutate(slide = factor(slide)) %>%
  # Uncomment below line if analysing hNSC too
  # mutate(cell_line = factor(cell_line, levels = c("hMSC", "hNSC"))) %>%
  mutate(timepoint = factor(timepoint, levels = c("early", "late"))) %>%
  # Uncomment below line if analysing hNSC too
  # mutate(passage = factor(passage, levels = c("p6", "p27"))) %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "heparin")))


# Create limma EList object containing beta values

quant_beta_vals <- new("EList")
quant_beta_vals$E <- getM(quant_ratioset_funnorm_filter)
quant_beta_vals$targets <- table_design
quant_beta_vals$genes <- quant_ratioset_funnorm_filter@rowRanges

# Define design matrix for limma
## Simple design allowing only pairwise comparisons
## Treat time points and treatments (condition) as fixed effects
## Treat cell line as an additive factor
## Include batch as an additive factor

design <- model.matrix(~ timepoint*treatment,
                       data = table_design)

colnames(design) <- make.names(colnames(design))

# Define contrasts

matrix_contrasts <- makeContrasts(
  P13vsP5 = timepointlate - 0,
  TvsUT = treatmentheparin - 0,
  levels = design
)

# Apply limma model fit

fit_contrasts <- lmFit(quant_beta_vals,
             design) %>%
  eBayes() %>%
  contrasts.fit(contrasts = matrix_contrasts) %>%
  eBayes()

