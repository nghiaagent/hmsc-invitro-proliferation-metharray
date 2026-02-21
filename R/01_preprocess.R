here::i_am("R/01_preprocess.R")

####################
# Load, process, and filter MethylationEPICv1 data
####################

# Import packages
library(dplyr)
library(here)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(tidyverse)

# Load samples
targets <- read.csv(
  here::here(
    "input",
    "annotation",
    "sample_sheet_hmsc.csv"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Basename = here::here(
      getwd(),
      "input",
      "dnam_epic_v1",
      slide,
      str_c(slide, array, sep = "_")
    ),
    # Order factors
    cell_line = cell_line %>%
      factor(levels = c("hMSC", "hNSC")),
    timepoint = timepoint %>%
      factor(levels = c("early", "late")),
    treatment = treatment %>%
      factor(levels = c("untreated", "heparin", "neurosphere"))
  ) %>%
  arrange(cell_line, timepoint, treatment) %>%
  # Order combined factors
  mutate(
    condition_notreat = str_c(cell_line, passage) %>%
      factor() %>%
      fct_inorder(),
    condition_nocell = str_c(timepoint, treatment) %>%
      factor() %>%
      fct_inorder(),
    sample_name = factor(sample_name) %>%
      fct_inorder()
  )

quant_rg <- read.metharray.exp(
  targets = targets,
  extended = TRUE,
  verbose = TRUE
)

# Load cross reactive and variant probes
list_pidsley_probes <- map(
  list(
    epic.cross1 = "13059_2016_1066_MOESM1_ESM.csv",
    epic.variants1 = "13059_2016_1066_MOESM4_ESM.csv",
    epic.variants2 = "13059_2016_1066_MOESM5_ESM.csv",
    epic.variants3 = "13059_2016_1066_MOESM6_ESM.csv"
  ),
  \(x) read_csv(here("input", "annotation", x))
)

excl_probes <- list_pidsley_probes %$%
  c(
    as.character(.$epic.cross1[[1]]),
    as.character(.$epic.variants1$PROBE),
    as.character(.$epic.variants2$PROBE),
    as.character(.$epic.variants3$PROBE)
  ) %>%
  unique()

# EXCLUDE BAD SAMPLES by detection p-val
## Calculate mean detection p-val in each sample
## Remove samples with mean detection p-val > 0.01
keep <- colMeans(detectionP(quant_rg)) < 0.01
quant_rg <- quant_rg[, keep]

# NORMALISATION
#### Quantile if not so different tissue types
#### Funnorm if very different tissue types.
#### Sometimes quantile doesnt work for plotting
#### (NA probes in some cases) so redo with Funnorm
quant_ratioset_funnorm <- quant_rg %>%
  preprocessFunnorm(ratioConvert = FALSE) %>%
  addQC(
    .,
    getQC(.)
  ) %>%
  ratioConvert(
    what = "both",
    keepCN = TRUE
  )

quant_mset_none <- quant_rg %>%
  preprocessRaw() %>%
  addQC(
    .,
    getQC(.)
  )

# Keep only probes detected in all samples
# Keep only probes not associated with snps
det_p <- detectionP(quant_rg)[
  match(
    featureNames(quant_ratioset_funnorm),
    rownames(detectionP(quant_rg))
  ),
]

keep <- rowSums(det_p < 0.01) == ncol(quant_ratioset_funnorm)
table(keep)

quant_ratioset_funnorm_filter <- quant_ratioset_funnorm[keep, ] %>%
  dropLociWithSnps()

# Keep only probes that are not cross-reactive
keep <- !(featureNames(quant_ratioset_funnorm_filter) %in% excl_probes)
table(keep)

quant_ratioset_funnorm_filter <- quant_ratioset_funnorm_filter[keep, ]

# Save data
# RGset, Mset raw, Mset funnorm, Mset funnorm filtered, detection p-values, targets
## RGset
saveRDS(
  quant_rg,
  file = here::here(
    "output",
    "quant_rg.RDS"
  )
)

## Mset raw
saveRDS(
  quant_mset_none,
  file = here::here(
    "output",
    "quant_mset_none.RDS"
  )
)

## Ratioset funnorm
saveRDS(
  quant_ratioset_funnorm,
  file = here::here(
    "output",
    "quant_ratioset_funnorm.RDS"
  )
)

## Ratioset funnorm filtered
saveRDS(
  quant_ratioset_funnorm_filter,
  file = here::here(
    "output",
    "quant_ratioset_funnorm_filter.RDS"
  )
)

## Detection p vals
saveRDS(
  det_p,
  file = here::here(
    "output",
    "det_p.RDS"
  )
)

## Sample sheet
saveRDS(
  targets,
  file = here::here(
    "output",
    "targets.RDS"
  )
)
