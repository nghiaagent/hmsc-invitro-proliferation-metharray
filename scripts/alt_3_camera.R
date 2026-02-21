# Apply camera to a function similar to fgsea in mCSEA

# Import data

quant_beta_vals <- readRDS(file.path(
  "output",
  "data_meth",
  "quant_beta_vals_hmsc.RDS"
))

fit_contrasts <- readRDS(file.path("output", "data_meth", "fit_contrasts.RDS"))

t_timepoint <- readRDS(file.path("output", "data_meth", "t_timepoint.RDS"))

mcsea_timepoint <- readRDS(
  file = file.path("output", "data_dmr", "mcsea", "dmr_hMSC_timepoint.RDS")
)

cpg2region_promoters <- mcsea_timepoint$promoters_association %>%
  ids2indices(identifiers = fit_contrasts$genes@ranges@NAMES)

cpg2region_genes <- mcsea_timepoint$genes_association %>%
  ids2indices(identifiers = fit_contrasts$genes@ranges@NAMES)

cpg2region_CGI <- mcsea_timepoint$CGI_association %>%
  ids2indices(identifiers = fit_contrasts$genes@ranges@NAMES)

run_camera <- function(fit, coefficient) {
  name_output <- as.character(colnames(matrix_contrasts)[coefficient])

  path_output <- file.path(getwd(), 'output', 'data_dmr', 'camera', name_output)

  message(str_c("Output camera results to", path_output, sep = " "))

  if (!dir.exists(path_output)) {
    dir.create(path_output, recursive = TRUE)
  }

  # Run camera on promoters

  message(str_c("Running camera for", name_output, sep = " "))

  camera_promoters <- camera(
    quant_beta_vals$E,
    cpg2region_promoters,
    fit_contrasts$design,
    fit_contrasts$contrasts[, coefficient],
    sort = FALSE,
    inter.gene.cor = NULL
  ) %>%
    mutate(name = rownames(.))

  camera_genes <- camera(
    quant_beta_vals$E,
    cpg2region_genes,
    fit_contrasts$design,
    fit_contrasts$contrasts[, coefficient],
    sort = FALSE,
    inter.gene.cor = NULL
  ) %>%
    mutate(name = rownames(.))

  camera_CGI <- camera(
    quant_beta_vals$E,
    cpg2region_CGI,
    fit_contrasts$design,
    fit_contrasts$contrasts[, coefficient],
    sort = FALSE,
    inter.gene.cor = NULL
  ) %>%
    mutate(name = rownames(.))

  camera_list <- list(
    "promoters" = camera_promoters,
    "genes" = camera_genes,
    "CGI" = camera_CGI
  )

  saveRDS(camera_list, file = file.path(path_output, "camera_results.RDS"))

  return(camera_list)
}

dmr_camera <- purrr::map(1:ncol(fit_contrasts$contrasts), \(x) {
  run_camera(fit_contrasts, x)
})

names(dmr_camera) <- colnames(dmr_camera)

saveRDS(dmr_camera, file.path('output', 'data_dmr', 'camera', "camera_all.RDS"))
