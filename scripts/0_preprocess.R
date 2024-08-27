# Load samples

targets <- read.csv(file.path("input", "annotation", "sample_sheet_hmsc.csv"),
                    stringsAsFactors = FALSE) %>%
  mutate(Basename = file.path(
    getwd(),
    "input",
    "dnam_epic_v1",
    slide,
    str_c(slide, array, sep = "_")
  )) %>%
  # Order factors
  mutate(cell_line = factor(cell_line, levels = c("hMSC", "hNSC"))) %>%
  mutate(timepoint = factor(timepoint, levels = c("early", "late"))) %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "heparin", "neurosphere"))) %>%
  arrange(cell_line, timepoint, treatment) %>%
  
  # Order combined factors
  mutate(condition_notreat = str_c(cell_line, passage) %>% factor() %>% fct_inorder) %>%
  mutate(condition_nocell = str_c(timepoint, treatment) %>% factor() %>% fct_inorder) %>%
  mutate(sample_name = factor(sample_name) %>% fct_inorder())

quant_rg <- read.metharray.exp(targets = targets,
                               extended = TRUE,
                               verbose = TRUE)

# Load cross reactive and variant probes

epic.cross1 <- read.csv(
  file.path(
    "input", "annotation",
    pattern = "13059_2016_1066_MOESM1_ESM.csv"
  )
)
epic.cross2 <- read.csv(
  file.path(
    "input", "annotation",
    pattern = "13059_2016_1066_MOESM2_ESM.csv"
  )
)
epic.cross3 <- read.csv(
  file.path(
    "input", "annotation",
    pattern = "13059_2016_1066_MOESM3_ESM.csv"
  )
)
epic.variants1 <- read.csv(
  file.path(
    "input", "annotation",
    pattern = "13059_2016_1066_MOESM4_ESM.csv"
  )
)
epic.variants2 <- read.csv(
  file.path(
    "input", "annotation",
    pattern = "13059_2016_1066_MOESM5_ESM.csv"
  )
)
epic.variants3 <- read.csv(
  file.path(
    "input", "annotation",
    pattern = "13059_2016_1066_MOESM6_ESM.csv"
  )
)

excl_probes <- c(
  as.character(epic.cross1$X),
  as.character(epic.variants1$PROBE),
  as.character(epic.variants2$PROBE),
  as.character(epic.variants3$PROBE)
) %>% unique()

# EXCLUDE BAD SAMPLES by detection p-val
## Calculate mean detection p-val in each sample
## Remove samples with mean detection p-val > 0.01

keep <- colMeans(detectionP(quant_rg)) < 0.01
quant_rg <- quant_rg[, keep]

# NORMALISATION
#### Quantile if not so different tissue types, Funnorm if very different tissue types.
#### Sometimes quantile doesnt work for plotting (NA probes in some cases) so redo with Funnorm

quant_ratioset_funnorm <- preprocessFunnorm(quant_rg, ratioConvert = FALSE) %>%
  ratioConvert(what = "both", 
               keepCN = TRUE)
quant_mset_none <- preprocessRaw(quant_rg)

# Quality control

## Draw QC plots

plotQC(getQC(quant_mset_none))

densityPlot(
  quant_mset_none,
  sampGroups = str_c(
    quant_mset_none@colData$cell_line,
    quant_mset_none@colData$timepoint,
    sep = "_"
  )
)

# Keep only probes detected in all samples
# Keep only probes not associated with snps

detP <- detectionP(quant_rg)[match(featureNames(quant_ratioset_funnorm), rownames(detectionP(quant_rg))), ]
keep <- rowSums(detP < 0.01) == ncol(quant_ratioset_funnorm)

table(keep)

quant_ratioset_funnorm_filter <- quant_ratioset_funnorm[keep, ] %>%
  dropLociWithSnps()

# Keep only probes that are not cross-reactive

keep <- !(featureNames(quant_ratioset_funnorm_filter) %in% excl_probes)

table(keep)

quant_ratioset_funnorm_filter <- quant_ratioset_funnorm_filter[keep, ]

# Save data

saveRDS(quant_ratioset_funnorm_filter, file = file.path("output",
                                                        "quant_ratioset_funnorm_filter.RDS"))

saveRDS(targets, file = file.path("output",
                                  "targets.RDS"))