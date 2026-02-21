here::i_am("R/05_post_correlate_gene_expression.R")

####################
# Correlate DNA methylation statistics with logFC
# Using mCSEA builtin function
####################

# Import packages
library(here)
library(mCSEA)
library(tidyverse)

# Correlate DNA methylation statistics with logFC
# DNA methylation statistic is defined (following mCSEA pub) as the mean beta-value of CpG probes forming the leading edge in each DMR.

# Define function that calculates mean delta-beta of leading edges (provided as a column in an mCSEAResults object)

calc_deltabeta <- function(leadingEdge, deltabetas) {
  map(
    .x = leadingEdge,
    .f = \(x) {
      strsplit(x, ", ") %>%
        unlist()
    },
    .progress = TRUE
  ) %>%
    map(
      .x = .,
      .f = \(x) {
        deltabetas %>%
          .[names(.) %in% x] %>%
          mean()
      },
      .progress = TRUE
    ) %>%
    unlist()
}

# Define function that tabulates methylation and expression of genes with DMRs of selected type
## Includes gg correlation plot too!

integrate_meth_expr <- function(
  coef,
  mCSEAResults,
  deltabetas,
  region = c("promoters", "genes")
) {
  table_expr <- topTable(fit_contrasts, coef = coef, number = Inf) %>%
    dplyr::select(ENTREZID, GENENAME, logFC)

  if (length(mCSEAResults[[region]] %>% filter(padj < 0.05)) < 1) {
    stop("There are no DMRs in this list")
  }

  data <- mCSEAResults[[region]] %>%
    filter(padj < 0.05) %>%
    mutate(GENENAME = rownames(.)) %>%
    select(!c("log2err", "pval")) %>%
    left_join(., table_expr, by = join_by(GENENAME == GENENAME)) %>%
    drop_na(logFC) %>%
    mutate(leadingEdge_mean = calc_deltabeta(leadingEdge, deltabetas))

  corrplot <- ggplot(data, aes(x = logFC, y = NES)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(
      "text",
      x = 0,
      y = 0.3,
      label = paste0(
        "r = ",
        round(
          cor(
            data$logFC,
            data$NES
          ),
          5
        )
      ),
      hjust = 0
    ) +
    annotate(
      "text",
      x = 0,
      y = -0.3,
      label = paste0(
        "p = ",
        round(
          cor.test(data$logFC, data$NES)$p.value,
          5
        )
      ),
      hjust = 0
    )
  return(list(data = data, corrplot = corrplot))
}

# Load data

mcsea_heparin <- readRDS(
  file = here("output", "data_dmr", "mcsea_deltabeta", "dmr_hMSC_heparin.RDS")
)
mcsea_timepoint <- readRDS(
  file = here("output", "data_dmr", "mcsea_deltabeta", "dmr_hMSC_timepoint.RDS")
)
fit_contrasts <- readRDS(
  file = here("output", "data_expression", "post_DGE", "fit_contrasts.RDS")
)
quant_ratioset_funnorm_filter <- readRDS(
  file = here("output", "quant_ratioset_funnorm_filter.RDS")
)
quant_m_vals <- readRDS(
  file = here("output", "data_meth", "limma", "quant_m_vals_hmsc.RDS")
)
beta_sel <- quant_ratioset_funnorm_filter[, c(
  "200654430047_R01C01",
  "200654430047_R02C01",
  "200654430047_R03C01",
  "200654430047_R04C01"
)] %>%
  getBeta()

# Extract beta values
# Calculate differences in beta values between P5 UT and P13 UT
# hMSC early UT: 200654430047_R01C01
# hMSC early T : 200654430047_R02C01
# hMSC late  UT: 200654430047_R03C01

beta <- list(
  # Define samples to extract
  early_untreated = "200654430047_R01C01",
  early_treated = "200654430047_R02C01",
  late_untreated = "200654430047_R03C01"
) %>%
  purrr::map(.x = ., .f = \(x) {
    quant_ratioset_funnorm_filter[, x] %>%
      getBeta() %>%
      .[, 1]
  })

deltabeta_timepoint <- beta$late_untreated - beta$early_untreated
deltabeta_treatment <- beta$early_treated - beta$early_untreated

# Make logFC vs NES table for between passages data

integrate_timepoint_promoters <- integrate_meth_expr(
  coef = 15,
  mCSEAResults = mcsea_timepoint,
  deltabetas = deltabeta_timepoint,
  region = "promoters"
)

integrate_timepoint_genebodies <- integrate_meth_expr(
  coef = 15,
  mCSEAResults = mcsea_timepoint,
  deltabetas = deltabeta_timepoint,
  region = "genes"
)

# Make logFC vs NES table for between treatments data

integrate_heparin_genebodies <- integrate_meth_expr(
  coef = 1,
  mCSEAResults = mcsea_heparin,
  deltabetas = deltabeta_treatment,
  region = "genes"
)
