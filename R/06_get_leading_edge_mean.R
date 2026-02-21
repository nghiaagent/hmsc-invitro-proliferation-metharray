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
    region = c("promoters", "genes")) {
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
            label = paste0("r = ", round(cor(
                data$logFC, data$NES
            ), 5)),
            hjust = 0
        ) +
        annotate(
            "text",
            x = 0,
            y = -0.3,
            label = paste0("p = ", round(
                cor.test(data$logFC, data$NES)$p.value, 5
            )),
            hjust = 0
        )
    return(list(data = data, corrplot = corrplot))
}

# Load data
results_mcsea <- readRDS(
    file = here(
        "output",
        "data_dmr",
        "results_mcsea.RDS"
    )
)

quant_ratioset_funnorm_filter <- readRDS(
    file = here(
        "output",
        "quant_ratioset_funnorm_filter.RDS"
    )
)

## Rearrange mCSEA results
results_mcsea <- list(
    timepoint_promoters = results_mcsea$timepoint$promoters,
    timepoint_genes = results_mcsea$timepoint$genes,
    treatment_early_promoters = results_mcsea$treatment_early$promoters,
    treatment_early_genes = results_mcsea$treatment_early$genes,
    treatment_late_promoters = results_mcsea$treatment_late$promoters,
    treatment_late_genes = results_mcsea$treatment_late$genes
)

# Extract beta values
# Calculate differences in beta values between P5 UT and P13 UT
# hMSC early UT: 200654430047_R01C01
# hMSC early T : 200654430047_R02C01
# hMSC late  UT: 200654430047_R03C01

beta <- map(
    list(
        early_untreated = "200654430047_R01C01",
        early_treated = "200654430047_R02C01",
        late_untreated = "200654430047_R03C01",
        late_treated = "200654430047_R04C01"
    ),
    \(x) getBeta(quant_ratioset_funnorm_filter)[, x]
)

deltabeta <- list(
    timepoint_promoters = beta[["late_untreated"]] - beta[["early_untreated"]],
    timepoint_genes = beta[["late_untreated"]] - beta[["early_untreated"]],
    treatment_early_promoters = beta[["early_treated"]] - beta[["early_untreated"]],
    treatment_early_genes = beta[["early_treated"]] - beta[["early_untreated"]],
    treatment_late_promoters = beta[["late_treated"]] - beta[["late_untreated"]],
    treatment_late_genes = beta[["late_treated"]] - beta[["late_untreated"]]
)

# Get mean deltabeta of leading edge
results_mcsea <- map2(
    results_mcsea, deltabeta,
    \(results, deltabeta) {
        results %>%
            mutate(leadingEdge_mean = calc_deltabeta(leadingEdge, deltabeta))
    }
)

# Save data
## RDS
saveRDS(results_mcsea,
    file = file.path(
        "output",
        "data_dmr",
        "results_mcsea_withmean.RDS"
    )
)

## Save xlsx
openxlsx::write.xlsx(
    results_mcsea,
    file = file.path(
        "output",
        "data_dmr",
        "results_mcsea_withmean.xlsx"
    ),
    asTable = TRUE,
    rowNames = TRUE
)