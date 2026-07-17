# ── 07 · ATTITUDES ─────────────────────────────────────────────────────────────
# Two things: (1) transition-matrix figures for the attitude-stratified matrices
# built in 02 (pooled and three-cohort), and (2) descriptive % conservative
# trends by religious tradition and birth cohort (childhood and current). The
# trend panels are not transition matrices — they are computed directly from the
# cleaned data.
#
# Input:  data/derived/matrices.rds, data/derived/gss_clean.rds
# Output: output/figures/attitude/{pooled,3coh,cohort_trends,cohort_trends16}/*.png

library(dplyr)
library(tidyr)
library(patchwork)
source("code/utils.R")

matrices = readRDS("data/derived/matrices.rds")
clean    = readRDS("data/derived/gss_clean.rds")
data     = clean$data

P_list_att_pooled   = matrices$att_pooled$P
pi0_list_att_pooled = matrices$att_pooled$pi0
pooled_summary      = matrices$att_pooled$summary

# ── ATTITUDE FIGURES ──────────────────────────────────────────────────────────

att_labels = c(
  evolved_bin  = "Evolution",
  abany_bin    = "Abortion (any reason)",
  homosex_bin  = "Homosexuality",
  premarsx_bin = "Premarital sex",
  pornlaw_bin  = "Pornography law",
  cappun_bin   = "Capital punishment"
)

dir.create("output/figures/attitude/pooled", recursive = TRUE, showWarnings = FALSE)

for (v in att_vars) {
  for (grp_lbl in c("conservative", "liberal")) {
    key = paste(v, grp_lbl, sep = "_")
    if (is.null(P_list_att_pooled[[key]])) next
    P       = P_list_att_pooled[[key]]
    pi0     = pi0_list_att_pooled[[key]]
    pistar  = pi_star(P)
    med_coh = pooled_summary$median_cohort[pooled_summary$variable == v &
                                             pooled_summary$group    == grp_lbl]
    ttl     = paste0(att_labels[v], " — ", tools::toTitleCase(grp_lbl), " (pooled)")
    p       = make_combined(P, pi0, pistar, levels = rel_level_order, title_str = ttl) +
                patchwork::plot_annotation(subtitle = paste0("Median birth cohort: ", med_coh))
    ggsave(paste0("output/figures/attitude/pooled/", key, ".png"),
           p, width = 10, height = 7, dpi = 200)
  }
}

# ── ATTITUDE TRENDS BY RELTRAD AND BIRTH COHORT ──────────────────────────────
# Conservative positions: evolved_bin == 0, abany_bin == 0, homosex_bin == 0

dir.create("output/figures/attitude/cohort_trends", recursive = TRUE, showWarnings = FALSE)

reltrad_colors_att = c(
  "Catholic"    = "#0072B2",
  "Evangelical" = "#D55E00",
  "Mainline"    = "#009E73",
  "Other"       = "#CC79A7",
  "None"        = "#999999"
)

att_cohort_df = data |>
  filter(!is.na(reltrad_alt), cohort_5 >= 1942.5, cohort_5 <= 1982.5) |>
  pivot_longer(
    cols      = c(evolved_bin, abany_bin, homosex_bin, cappun_bin),
    names_to  = "attitude",
    values_to = "liberal"
  ) |>
  filter(!is.na(liberal)) |>
  group_by(attitude, reltrad_alt, cohort_5) |>
  summarise(
    n                = n(),
    pct_conservative = mean(liberal == 0L) * 100,
    se               = sqrt(pct_conservative / 100 * (1 - pct_conservative / 100) / n) * 100,
    .groups          = "drop"
  ) |>
  mutate(
    attitude = factor(attitude,
      levels = c("evolved_bin", "abany_bin", "homosex_bin", "cappun_bin"),
      labels = c("Evolution (% deny)", "Abortion (% oppose)", "Homosexuality (% morally wrong)",
                 "Capital punishment (% favor)")
    ),
    reltrad_alt = factor(reltrad_alt,
      levels = c("catholic", "evangelical", "mainline", "other", "none"),
      labels = c("Catholic", "Evangelical", "Mainline", "Other", "None")
    )
  )

make_att_cohort_plot = function(att_label) {
  att_cohort_df |>
    filter(attitude == att_label) |>
    ggplot(aes(x = cohort_5, y = pct_conservative, color = reltrad_alt, fill = reltrad_alt, group = reltrad_alt)) +
    geom_ribbon(aes(ymin = pct_conservative - 1.96 * se, ymax = pct_conservative + 1.96 * se),
                alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    scale_color_manual(values = reltrad_colors_att, name = NULL) +
    scale_fill_manual(values = reltrad_colors_att, name = NULL) +
    scale_x_continuous(breaks = seq(1940, 1980, by = 10)) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    labs(title = att_label, x = NULL, y = sub(".*\\((.+)\\)", "\\1", att_label)) +
    healy_theme +
    theme(plot.title = element_text(size = 12, face = "plain"))
}

p_evolved = make_att_cohort_plot("Evolution (% deny)")
p_abany   = make_att_cohort_plot("Abortion (% oppose)")
p_homosex = make_att_cohort_plot("Homosexuality (% morally wrong)")
p_cappun  = make_att_cohort_plot("Capital punishment (% favor)")

caption_str = "Source: GSS. 5-year birth cohort bins, 1940–1980. Ribbons show 95% CIs. Current religious affiliation (reltrad_alt)."

p_evolved = p_evolved + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)
p_abany   = p_abany   + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)
p_homosex = p_homosex + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)
p_cappun  = p_cappun  + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)

ggsave("output/figures/attitude/cohort_trends/att_trends_evolution.png",
       p_evolved, width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends/att_trends_abortion.png",
       p_abany,   width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends/att_trends_homosexuality.png",
       p_homosex, width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends/att_trends_cappun.png",
       p_cappun,  width = 7, height = 5, dpi = 200)

# ── ATTITUDE TRENDS BY RELTRAD16 (CHILDHOOD RELIGION) AND BIRTH COHORT ────────

dir.create("output/figures/attitude/cohort_trends16", recursive = TRUE, showWarnings = FALSE)

att_cohort_df16 = data |>
  filter(!is.na(reltrad16_alt), cohort_5 >= 1942.5, cohort_5 <= 1982.5) |>
  pivot_longer(
    cols      = c(evolved_bin, abany_bin, homosex_bin, cappun_bin),
    names_to  = "attitude",
    values_to = "liberal"
  ) |>
  filter(!is.na(liberal)) |>
  group_by(attitude, reltrad16_alt, cohort_5) |>
  summarise(
    n                = n(),
    pct_conservative = mean(liberal == 0L) * 100,
    se               = sqrt(pct_conservative / 100 * (1 - pct_conservative / 100) / n) * 100,
    .groups          = "drop"
  ) |>
  mutate(
    attitude = factor(attitude,
      levels = c("evolved_bin", "abany_bin", "homosex_bin", "cappun_bin"),
      labels = c("Evolution (% deny)", "Abortion (% oppose)", "Homosexuality (% morally wrong)",
                 "Capital punishment (% favor)")
    ),
    reltrad16_alt = factor(reltrad16_alt,
      levels = c("catholic", "evangelical", "mainline", "other", "none"),
      labels = c("Catholic", "Evangelical", "Mainline", "Other", "None")
    )
  )

make_att_cohort_plot16 = function(att_label) {
  att_cohort_df16 |>
    filter(attitude == att_label) |>
    ggplot(aes(x = cohort_5, y = pct_conservative, color = reltrad16_alt, fill = reltrad16_alt, group = reltrad16_alt)) +
    geom_ribbon(aes(ymin = pct_conservative - 1.96 * se, ymax = pct_conservative + 1.96 * se),
                alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    scale_color_manual(values = reltrad_colors_att, name = NULL) +
    scale_fill_manual(values = reltrad_colors_att, name = NULL) +
    scale_x_continuous(breaks = seq(1940, 1980, by = 10)) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    labs(title = att_label, x = NULL, y = sub(".*\\((.+)\\)", "\\1", att_label)) +
    healy_theme +
    theme(plot.title = element_text(size = 12, face = "plain"))
}

caption_str16 = "Source: GSS. 5-year birth cohort bins, 1940–1980. Ribbons show 95% CIs. Childhood religious affiliation (reltrad16_alt)."

p_evolved16 = make_att_cohort_plot16("Evolution (% deny)") +
  labs(x = "Birth cohort (5-year bin)") + plot_annotation(caption = caption_str16)
p_abany16   = make_att_cohort_plot16("Abortion (% oppose)") +
  labs(x = "Birth cohort (5-year bin)") + plot_annotation(caption = caption_str16)
p_homosex16 = make_att_cohort_plot16("Homosexuality (% morally wrong)") +
  labs(x = "Birth cohort (5-year bin)") + plot_annotation(caption = caption_str16)
p_cappun16  = make_att_cohort_plot16("Capital punishment (% favor)") +
  labs(x = "Birth cohort (5-year bin)") + plot_annotation(caption = caption_str16)

ggsave("output/figures/attitude/cohort_trends16/att_trends16_evolution.png",
       p_evolved16, width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends16/att_trends16_abortion.png",
       p_abany16,   width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends16/att_trends16_homosexuality.png",
       p_homosex16, width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends16/att_trends16_cappun.png",
       p_cappun16,  width = 7, height = 5, dpi = 200)
