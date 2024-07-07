# Load data

quant_ratioset_funnorm_filter <- readRDS(file = file.path("output",
                                                        "quant_ratioset_funnorm_filter.RDS"))

quant_ratioset_funnorm_filter[,cell_line = "hMSC"]

# Make design matrix

table_design <- quant_ratioset_funnorm_filter@colData %>%
  as_data_frame() %>%
  mutate(sample_name = factor(
    sample_name,
    levels = c(
      "hMSC_p5_untreated",
      "hMSC_p5_heparin",
      "hMSC_p5_neurosphere",
      "hMSC_p13_untreated",
      "hMSC_p13_heparin",
      "hMSC_p13_neurosphere",
      "hNSC_p6_untreated",
      "hNSC_p6_heparin",
      "hNSC_p27_untreated",
      "hNSC_p27_heparin"
    )
  )) %>%
  mutate(slide = factor(slide)) %>%
  mutate(cell_line = factor(cell_line, levels = c("hMSC", "hNSC"))) %>%
  mutate(timepoint = factor(timepoint, levels = c("early", "late"))) %>%
  mutate(passage = factor(passage, levels = c("p5", "p13", "p27"))) %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "heparin", "neurosphere")))

# Define design matrix for limma

## Simple design allowing only pairwise comparisons
## Treat time points and treatments (condition) as fixed effects
## Treat cell line as an additive factor
## Include batch as an additive factor

design <- model.matrix(~ sample_name,
                       data = table_design)

colnames(design) <- make.names(colnames(design))

# Create limma EList object containing beta values

quant_beta_vals <- new("EList")
quant_beta_vals$E <- getBeta(quant_ratioset_funnorm_filter)
quant_beta_vals$targets <- table_design
quant_beta_vals$genes <- quant_ratioset_funnorm_filter@rowRanges

# Apply limma model fit

fit <- lmFit(quant_beta_vals,
             design) %>%
  eBayes()