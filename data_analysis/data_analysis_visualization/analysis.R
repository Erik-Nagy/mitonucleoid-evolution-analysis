# 11–14 plots for P-site characterisation in mitochondrial nucleoid proteins.
# Requires: psite_long, save_fig(), COL_AA, COL_DISORDER, COL_SASA  (00_load_data.R)

library(tidyverse)
library(ggpubr)

n_sites <- nrow(psite_long)

# 1. Total conservation — histogram
p01 <- psite_long %>%
  filter(!is.na(exact_cons)) %>%
  ggplot(aes(x = exact_cons)) +
  geom_histogram(binwidth = 12.5, fill = "#4393c3", colour = "white",
                 boundary = 0, na.rm = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 12.5), limits = c(-6, 106)) +
  labs(
    title    = "Evolutionary conservation of P-sites",
    subtitle = sprintf("n = %d P-sites, up to 8 homolog species", n_sites),
    x        = "Exact Conservation (%)",
    y        = "Number of P-sites"
  )
save_fig(p01, "01_conservation_total", width = 6, height = 5)


# 2. Conservation by residue type (pS / pT / pY)
p02 <- psite_long %>%
  filter(!is.na(residue_type), !is.na(exact_cons)) %>%
  ggplot(aes(x = residue_type, y = exact_cons, fill = residue_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.7, colour = "gray30", na.rm = TRUE) +
  scale_fill_manual(values = COL_AA, guide = "none") +
  labs(
    title = "Exact conservation by phosphorylated residue type",
    x     = "Residue",
    y     = "Exact Conservation (%)"
  )
save_fig(p02, "02_conservation_by_residue", width = 6, height = 5)


# 3. Total ordered / disordered — bar chart

dis_total <- psite_long %>%
  filter(!is.na(disorder_state)) %>%
  count(disorder_state, name = "n") %>%
  mutate(pct = 100 * n / sum(n))

p03 <- dis_total %>%
  ggplot(aes(x = disorder_state, y = n, fill = disorder_state)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = COL_DISORDER, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(
    title    = "Ordered vs disordered P-sites",
    subtitle = sprintf("n = %d P-sites", n_sites),
    x        = NULL,
    y        = "Number of P-sites"
  ) +
  theme(plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "mm"))
save_fig(p03, "03_disorder_total", width = 5, height = 5)


# 4. Disorder score by residue type (pS / pT / pY)
p04 <- psite_long %>%
  filter(!is.na(residue_type), !is.na(disorder_score)) %>%
  ggplot(aes(x = residue_type, y = disorder_score, fill = residue_type)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "gray50",
             linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.7, colour = "gray30", na.rm = TRUE) +
  scale_fill_manual(values = COL_AA, guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.5), clip = "off") +
  annotate("text", x = 0.55, y = 0.53, label = "threshold 0.5",
           hjust = 0, size = 3, colour = "gray50") +
  stat_compare_means(
    comparisons = list(c("pS", "pT"), c("pS", "pY"), c("pT", "pY")),
    method = "wilcox.test", label = "p.signif", tip.length = 0.01
  ) +
  labs(
    title = "Disorder score by phosphorylated residue type",
    x     = "Residue",
    y     = "Disorder score (0 = ordered, 1 = disordered)"
  )
save_fig(p04, "04_disorder_by_residue", width = 5, height = 5)


# 5. Total exposed / buried — bar chart
exp_total <- psite_long %>%
  filter(!is.na(sasa_location)) %>%
  count(sasa_location, name = "n") %>%
  mutate(pct = 100 * n / sum(n))

p05 <- exp_total %>%
  ggplot(aes(x = sasa_location, y = n, fill = sasa_location)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = COL_SASA, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.40))) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Exposed vs buried P-sites (SASA threshold: 20 Å²)",
    subtitle = sprintf("n = %d P-sites", n_sites),
    x        = NULL,
    y        = "Number of P-sites"
  ) +
  theme(plot.margin = margin(t = 15, r = 25, b = 5, l = 5, unit = "mm"))
save_fig(p05, "05_exposure_total", width = 6, height = 5)


# 6. Exposed / buried by residue type (pS / pT / pY)
exp_by_aa <- psite_long %>%
  filter(!is.na(residue_type), !is.na(sasa_location)) %>%
  count(residue_type, sasa_location, name = "n") %>%
  group_by(residue_type) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p06 <- exp_by_aa %>%
  ggplot(aes(x = residue_type, y = pct, fill = sasa_location)) +
  geom_col(colour = "white", linewidth = 0.4, width = 0.6) +
  geom_text(aes(label = ifelse(pct > 6, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5),
            size = 3.8, colour = "white", fontface = "bold") +
  scale_fill_manual(values = COL_SASA, name = "Location") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Surface exposure by phosphorylated residue type",
    x     = "Residue",
    y     = "Percentage of P-sites"
  ) +
  theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5, unit = "mm"))
save_fig(p06, "06_exposure_by_residue", width = 6, height = 5)


# 7. Total pLDDT — histogram with AF2 confidence zones
zone_df <- data.frame(
  xmin = c(0,  50, 70, 90),
  xmax = c(50, 70, 90, 100),
  fill = c("#d73027", "#fc8d59", "#91cf60", "#1a9850")
)

p07 <- ggplot() +
  geom_rect(data = zone_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                fill = fill),
            alpha = 0.12, inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_histogram(data = filter(psite_long, !is.na(plddt)), aes(x = plddt),
                 binwidth = 5, fill = "gray30", colour = "white", alpha = 0.85,
                 boundary = 0, na.rm = TRUE) +
  scale_x_continuous(breaks = c(0, 50, 70, 90, 100)) +
  labs(
    title    = "AlphaFold2 pLDDT scores at P-sites",
    subtitle = sprintf("n = %d P-sites", n_sites),
    x        = "pLDDT Score",
    y        = "Number of P-sites"
  ) +
  annotate("text", x = c(25, 60, 80, 95), y = Inf,
           label = c("Very low", "Low", "Confident", "Very\nhigh"),
           vjust = 1.4, size = 2.8, colour = "gray30")
save_fig(p07, "07_plddt_total", width = 6, height = 5)


# 8. pLDDT by residue type (pS / pT / pY)
p08 <- psite_long %>%
  filter(!is.na(residue_type), !is.na(plddt)) %>%
  ggplot(aes(x = residue_type, y = plddt, fill = residue_type)) +
  geom_hline(yintercept = c(50, 70, 90), linetype = "dashed",
             colour = "gray70", linewidth = 0.6) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.7, colour = "gray30", na.rm = TRUE) +
  scale_fill_manual(values = COL_AA, guide = "none") +
  scale_y_continuous(breaks = c(0, 50, 70, 90, 100)) +
  labs(
    title = "AlphaFold2 pLDDT score by phosphorylated residue type",
    x     = "Residue",
    y     = "pLDDT Score"
  )
save_fig(p08, "08_plddt_by_residue", width = 6, height = 5)


# 9. Conservation vs disorder — boxplot (ordered vs disordered)
p09 <- psite_long %>%
  filter(!is.na(disorder_state), !is.na(exact_cons)) %>%
  ggplot(aes(x = disorder_state, y = exact_cons, fill = disorder_state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.18, size = 2.2, alpha = 0.7, colour = "gray30", na.rm = TRUE) +
  scale_fill_manual(values = COL_DISORDER, guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.35, label.y = 105) +
  labs(
    title = "Conservation in ordered vs disordered P-sites",
    x     = NULL,
    y     = "Exact Conservation (%)"
  )
save_fig(p09, "09_conservation_vs_disorder", width = 5, height = 5)


# 10. Conservation vs surface exposure — boxplot
p10 <- psite_long %>%
  filter(!is.na(sasa_location), !is.na(exact_cons)) %>%
  ggplot(aes(x = sasa_location, y = exact_cons, fill = sasa_location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.18, size = 2.2, alpha = 0.7, colour = "gray30", na.rm = TRUE) +
  scale_fill_manual(values = COL_SASA, guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.35, label.y = 105) +
  labs(
    title = "Conservation vs surface exposure",
    x     = NULL,
    y     = "Exact Conservation (%)"
  )
save_fig(p10, "10_conservation_vs_exposure", width = 5, height = 5)


# 11. Disorder vs surface exposure — boxplot
p11 <- psite_long %>%
  filter(!is.na(sasa_location), !is.na(disorder_score)) %>%
  ggplot(aes(x = sasa_location, y = disorder_score, fill = sasa_location)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "gray60",
             linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.18, size = 2.2, alpha = 0.7, colour = "gray30", na.rm = TRUE) +
  scale_fill_manual(values = COL_SASA, guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.15), clip = "off") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.35, label.y = 1.08) +
  labs(
    title = "Disorder score vs surface exposure",
    x     = NULL,
    y     = "Disorder score"
  )
save_fig(p11, "11_disorder_vs_exposure", width = 5, height = 5)


# 12-14. pSTY vs background STY — disorder score, pLDDT, surface exposure
#         Requires all_sty_data.tsv from extract_all_sty_data.py

if (!is.null(all_sty)) {

  # Paired fill colours: background (light) vs P-site (dark) per residue class
  COL_STY <- c(
    "S"  = "#fcbba1", "pS" = "#cb181d",
    "T"  = "#9ecae1", "pT" = "#2171b5",
    "Y"  = "#fed99c", "pY" = "#d94701"
  )

  # Comparison pairs for stat_compare_means
  sty_pairs <- list(c("S", "pS"), c("T", "pT"), c("Y", "pY"))

  # 12. Disorder score: pSTY vs STY
  p12 <- all_sty %>%
    filter(!is.na(residue_type), !is.na(disorder_score)) %>%
    ggplot(aes(x = residue_type, y = disorder_score, fill = residue_type)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "gray50",
               linewidth = 0.8) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
    geom_jitter(width = 0.14, size = 1.6, alpha = 0.5, colour = "gray30",
                na.rm = TRUE) +
    scale_fill_manual(values = COL_STY, guide = "none") +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(ylim = c(0, 1.35), clip = "off") +
    annotate("text", x = 0.55, y = 0.53, label = "threshold 0.5",
             hjust = 0, size = 3, colour = "gray50") +
    stat_compare_means(
      comparisons = sty_pairs,
      method = "wilcox.test", label = "p.signif", tip.length = 0.01
    ) +
    labs(
      title    = "Disorder score: phosphorylated vs all S/T/Y residues",
      subtitle = sprintf("n = %d P-sites, %d background STY",
                         sum(all_sty$is_psite), sum(!all_sty$is_psite)),
      x        = "Residue type  (lighter = background, darker = P-site)",
      y        = "Disorder score (0 = ordered, 1 = disordered)"
    ) +
    theme(plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "mm"))
  save_fig(p12, "12_disorder_pSTY_vs_STY", width = 7, height = 6)

  # 13. pLDDT: pSTY vs STY
  p13 <- all_sty %>%
    filter(!is.na(residue_type), !is.na(plddt)) %>%
    ggplot(aes(x = residue_type, y = plddt, fill = residue_type)) +
    geom_hline(yintercept = c(50, 70, 90), linetype = "dashed",
               colour = "gray70", linewidth = 0.6) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
    geom_jitter(width = 0.14, size = 1.6, alpha = 0.5, colour = "gray30",
                na.rm = TRUE) +
    scale_fill_manual(values = COL_STY, guide = "none") +
    scale_y_continuous(breaks = c(0, 50, 70, 90, 100)) +
    coord_cartesian(ylim = c(0, 130), clip = "off") +
    stat_compare_means(
      comparisons = sty_pairs,
      method = "wilcox.test", label = "p.signif", tip.length = 0.01
    ) +
    labs(
      title    = "AlphaFold2 pLDDT: phosphorylated vs all S/T/Y residues",
      subtitle = sprintf("n = %d P-sites, %d background STY",
                         sum(all_sty$is_psite), sum(!all_sty$is_psite)),
      x        = "Residue type  (lighter = background, darker = P-site)",
      y        = "pLDDT Score"
    ) +
    theme(plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "mm"))
  save_fig(p13, "13_plddt_pSTY_vs_STY", width = 7, height = 6)

  # 14. Surface exposure: pSTY vs STY (stacked 100% bar)
  exp_sty <- all_sty %>%
    filter(!is.na(residue_type), !is.na(sasa_location)) %>%
    count(residue_type, sasa_location, name = "n") %>%
    group_by(residue_type) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  p14 <- exp_sty %>%
    ggplot(aes(x = residue_type, y = pct, fill = sasa_location)) +
    geom_col(colour = "white", linewidth = 0.4, width = 0.65) +
    geom_text(aes(label = ifelse(pct > 6, sprintf("%.0f%%", pct), "")),
              position = position_stack(vjust = 0.5),
              size = 3.5, colour = "white", fontface = "bold") +
    scale_fill_manual(values = COL_SASA, name = "Location") +
    scale_y_continuous(labels = scales::percent_format(scale = 1),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = "Surface exposure: phosphorylated vs all S/T/Y residues",
      subtitle = sprintf("n = %d P-sites, %d background STY (SASA threshold: 20 Å²)",
                         sum(all_sty$is_psite & !is.na(all_sty$sasa_location)),
                         sum(!all_sty$is_psite & !is.na(all_sty$sasa_location))),
      x        = "Residue type  (lighter shade = background, darker = P-site)",
      y        = "Percentage of residues"
    ) +
    theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5, unit = "mm"))
  save_fig(p14, "14_exposure_pSTY_vs_STY", width = 7, height = 6)

  message(sprintf("Figures 12-14 saved to %s/", FIG_DIR))
} else {
  message("Skipping plots 12-14 (all_sty_data.tsv not found — run 'make process').")
}

message(sprintf("All figures saved to %s/", FIG_DIR))
