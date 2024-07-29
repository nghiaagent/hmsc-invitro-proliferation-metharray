# Compare results of calculating DMR via t-statistics and deltabeta

mcsea_timepoint_tstat <- readRDS(file.path("output", "data_dmr", "mcsea", "dmr_hMSC_timepoint.RDS"))[["promoters"]] %>%
  mutate(gene = rownames(.)) %>%
  arrange(gene) %>%
  relocate(gene)

mcsea_timepoint_deltabeta <- readRDS(file.path(
  "output",
  "data_dmr",
  "mcsea_deltabeta",
  "dmr_hMSC_timepoint.RDS"
))[["promoters"]] %>%
  mutate(gene = rownames(.)) %>%
  arrange(gene) %>%
  relocate(gene)

df_cor_nes <- left_join(
  mcsea_timepoint_tstat,
  mcsea_timepoint_deltabeta,
  by = join_by(gene == gene),
  suffix = c("tstat", "deltabeta")
) %>%
  drop_na()

# Plot correlation of ES

ggplot(df_cor_nes, aes(x = EStstat, y = ESdeltabeta)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate(
    "text",
    x = 0,
    y = 1.25,
    label = paste0("r = ", round(
      cor(df_cor_nes$EStstat, df_cor_nes$ESdeltabeta), 2
    )),
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0,
    y = 1,
    label = paste0("p = ", round(
      cor.test(df_cor_nes$EStstat, df_cor_nes$ESdeltabeta)$p.value,
      3
    )),
    hjust = 0
  ) +
  theme_classic()

# Plot correlation of NES

ggplot(df_cor_nes, aes(x = NEStstat, y = NESdeltabeta)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate(
    "text",
    x = 2,
    y = 4,
    label = paste0("r = ", round(
      cor(df_cor_nes$NEStstat, df_cor_nes$NESdeltabeta), 2
    )),
    hjust = 0
  ) +
  annotate(
    "text",
    x = 2,
    y = 3.75,
    label = paste0("p = ", round(
      cor.test(df_cor_nes$NEStstat, df_cor_nes$NESdeltabeta)$p.value,
      3
    )),
    hjust = 0
  ) +
  theme_classic()

# Plot correlation of p-vals

ggplot(df_cor_nes, aes(x = pvaltstat, y = pvaldeltabeta)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate(
    "text",
    x = 0,
    y = 1.1,
    label = paste0("r = ", round(
      cor(df_cor_nes$pvaltstat, df_cor_nes$pvaldeltabeta), 2
    )),
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0,
    y = 1,
    label = paste0("p = ", round(
      cor.test(df_cor_nes$pvaltstat, df_cor_nes$pvaldeltabeta)$p.value,
      3
    )),
    hjust = 0
  ) +
  theme_classic()

# Venn diagram

ggVennDiagram(list(
  tstat = df_cor_nes[df_cor_nes$padjtstat < 0.05, ]$gene,
  beta =  df_cor_nes[df_cor_nes$padjdeltabeta < 0.05, ]$gene
))

# Save data

saveRDS(df_cor_nes, file = file.path("output",
                                     "data_dmr",
                                     "df_cor_nes.RDS"))
