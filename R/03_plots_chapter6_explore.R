# Load data
quant_ratioset_funnorm_filter <- readRDS(
  file = here::here(
    "output",
    "quant_ratioset_funnorm_filter.RDS"
  )
)

samples_msc <- c(
  early_untreated = "200654430047_R01C01",
  early_treated = "200654430047_R02C01",
  late_untreated = "200654430047_R03C01",
  late_treated = "200654430047_R04C01"
)

# Get annotation for probes
annotation <- quant_ratioset_funnorm_filter %>%
  getAnnotation() %>%
  as.data.frame()

# Create dataframe suitable for plotting probes by genomic region
# and relation to gene
df_beta <- getBeta(quant_ratioset_funnorm_filter) %>%
  as.data.frame() %>%
  # Select relevane samples
  dplyr::select(all_of(samples_msc)) %>%
  # Add annotation data
  rownames_to_column("Name") %>%
  left_join(
    annotation,
    by = join_by(Name == Name)
  ) %>%
  # Remove probes with no annotation
  drop_na(
    UCSC_RefGene_Group,
    Relation_to_Island
  ) %>%
  dplyr::filter(
    UCSC_RefGene_Group != "",
    Relation_to_Island != ""
  ) %>%
  # Calculate sd, filter for probes sd > 1
  rowwise() %>%
  mutate(
    beta_mean = mean(c(
      early_untreated,
      early_treated,
      late_untreated,
      late_treated
    )),
    beta_sd = sd(c(
      early_untreated,
      early_treated,
      late_untreated,
      late_treated
    ))
  ) %>%
  ungroup() %>%
  dplyr::filter(beta_sd > 0.1) %>%
  # Clean up probe sd data
  mutate(
    UCSC_RefGene_Group = UCSC_RefGene_Group %>%
      str_replace("\\;.*$", "")
  )

df_beta_long <- df_beta %>%
  # Pivot to long table
  pivot_longer(
    cols = c(
      early_untreated,
      early_treated,
      late_untreated,
      late_treated
    ),
    names_to = "condition_id",
    values_to = "beta"
  ) %>%
  # Create factors for easy plotting
  mutate(
    condition_id = condition_id %>%
      factor(
        levels = c(
          "early_untreated",
          "early_treated",
          "late_untreated",
          "late_treated"
        )
      ),
    UCSC_RefGene_Group = UCSC_RefGene_Group %>%
      factor(
        levels = c(
          "TSS1500",
          "TSS200",
          "5'UTR",
          "1stExon",
          "Body",
          "ExonBnd",
          "3'UTR"
        )
      ),
    Relation_to_Island = Relation_to_Island %>%
      factor(
        levels = c(
          "N_Shelf",
          "N_Shore",
          "Island",
          "S_Shore",
          "S_Shelf",
          "OpenSea"
        )
      )
  )

# Plot
plot_gene_model <- df_beta_long %>%
  ggplot(
    aes(
      x = condition_id,
      y = beta,
    )
  ) +
  geom_jitter(
    aes(
      fill = condition_id,
      colour = condition_id
    ),
    width = 0.1,
    size = 0.4,
    alpha = 0.2,
    stroke = 0.02
  ) +
  geom_flat_violin(
    aes(
      fill = condition_id,
      colour = condition_id
    ),
    position = position_nudge(x = 0.2)
  ) +
  geom_boxplot(
    aes(
      colour = condition_id
    ),
    outliers = FALSE,
    width = 0.2
  ) +
  theme_bw() +
  scale_colour_manual(
    values = palette_msc,
    labels = c(
      "early_untreated" = "P+5 Ctrl",
      "early_treated" = "P+5 Hep",
      "late_untreated" = "P+13 Ctrl",
      "late_treated" = "P+13 Hep"
    )
  ) +
  scale_fill_manual(
    values = palette_msc %>%
      lighten(0.4),
    labels = c(
      "early_untreated" = "P+5 Ctrl",
      "early_treated" = "P+5 Hep",
      "late_untreated" = "P+13 Ctrl",
      "late_treated" = "P+13 Hep"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "early_untreated" = "P+5 Ctrl",
      "early_treated" = "P+5 Hep",
      "late_untreated" = "P+13 Ctrl",
      "late_treated" = "P+13 Hep"
    )
  ) +
  labs(
    x = NULL,
    y = "Beta values"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
  ) +
  guides(colour = guide_legend(ncol = 4)) +
  facet_wrap(
    ~UCSC_RefGene_Group,
    ncol = 2
  )

# Plot
plot_island_model <- df_beta_long %>%
  ggplot(
    aes(
      x = condition_id,
      y = beta,
    )
  ) +
  geom_jitter(
    aes(
      fill = condition_id,
      colour = condition_id
    ),
    width = 0.1,
    size = 0.4,
    alpha = 0.2,
    stroke = 0.02
  ) +
  geom_flat_violin(
    aes(
      fill = condition_id,
      colour = condition_id
    ),
    position = position_nudge(x = 0.2)
  ) +
  geom_boxplot(
    aes(
      colour = condition_id
    ),
    outliers = FALSE,
    width = 0.2
  ) +
  theme_bw() +
  scale_colour_manual(
    values = palette_msc,
    labels = c(
      "early_untreated" = "P+5 Ctrl",
      "early_treated" = "P+5 Hep",
      "late_untreated" = "P+13 Ctrl",
      "late_treated" = "P+13 Hep"
    )
  ) +
  scale_fill_manual(
    values = palette_msc %>%
      lighten(0.4),
    labels = c(
      "early_untreated" = "P+5 Ctrl",
      "early_treated" = "P+5 Hep",
      "late_untreated" = "P+13 Ctrl",
      "late_treated" = "P+13 Hep"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "early_untreated" = "P+5 Ctrl",
      "early_treated" = "P+5 Hep",
      "late_untreated" = "P+13 Ctrl",
      "late_treated" = "P+13 Hep"
    )
  ) +
  labs(
    x = NULL,
    y = "Beta values"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
  ) +
  guides(colour = guide_legend(ncol = 4)) +
  facet_wrap(
    ~Relation_to_Island,
    ncol = 2
  )

# Get summary statistics
summary_gene_model <- df_beta %>%
  group_by(UCSC_RefGene_Group) %>%
  filter() %>%
  summarise(n = n())

summary_island_model <- df_beta %>%
  group_by(Relation_to_Island) %>%
  filter() %>%
  summarise(n = n())

# Save plots
ggsave(
  filename = "violin_gene_model.png",
  plot = plot_gene_model,
  path = here::here(
    "output",
    "plots_explore"
  ),
  scale = 0.8,
  width = 8,
  height = 10,
  units = "in",
  dpi = 144
)

ggsave(
  filename = "violin_island_model.png",
  plot = plot_island_model,
  path = here::here(
    "output",
    "plots_explore"
  ),
  scale = 0.8,
  width = 8,
  height = 10,
  units = "in",
  dpi = 144
)
