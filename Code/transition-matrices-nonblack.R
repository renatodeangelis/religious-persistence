library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gssr)
source("Code/utils.R")

reltrad_labels  = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")

reltrad_colors = c(
  catholic    = "#0072B2",
  evangelical = "#D55E00",
  mainline    = "#009E73",
  other       = "#CC79A7",
  none        = "#999999"
)
reltrad_labels_tc = c(
  catholic    = "Catholic",
  evangelical = "Evangelical",
  mainline    = "Mainline",
  other       = "Other",
  none        = "None"
)

# ── DATA (non-Black sample: race != 2) ───────────────────────────────────────

data(gss_all)
data_nb = gss_all |>
  select(year, cohort, sex, reltrad, reltrad16, region, born,
         evolved, abany, homosex, premarsx, pornlaw, cappun, cappun2, race) |>
  filter(!(year %in% c(1972, 2021))) |>
  filter(as.numeric(race) != 2) |>
  mutate(across(c(reltrad, reltrad16),
                ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  mutate(across(c(reltrad, reltrad16),
                ~ case_when(
                  . == "jewish"           ~ "other",
                  . == "black protestant" ~ "evangelical",
                  TRUE                    ~ .
                ),
                .names = "{.col}_alt")) |>
  mutate(age = year - cohort) |>
  filter(age > 25, cohort >= 1900) |>
  mutate(
    cohort_10 = floor((cohort - 1900) / 10) * 10 + 1900
  )

states_alt = sort(unique(c(data_nb$reltrad_alt, data_nb$reltrad16_alt)))
states_alt = states_alt[!is.na(states_alt)]

# ── DECADAL TRANSITION MATRICES (1940–1989, non-Black) ───────────────────────

cohorts_nb = c(1940, 1950, 1960, 1970, 1980)

P_list_nb      = list()
pi0_list_nb    = list()
pistar_list_nb = list()
n_list_nb      = list()

for (coh in cohorts_nb) {
  sub = data_nb[
    !is.na(data_nb$cohort_10)     & data_nb$cohort_10     == coh &
    !is.na(data_nb$reltrad16_alt) & !is.na(data_nb$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_nb[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_nb[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_nb[[key]] = pi_star(P_list_nb[[key]])
  n_list_nb[[key]]      = nrow(sub)
}

# Print matrices
for (key in names(P_list_nb)) {
  cat("\n── Cohort", key, "–", as.integer(key) + 9,
      "  (N =", n_list_nb[[key]], ") ──\n")
  print(round(P_list_nb[[key]], 3))
}

# ── IM MEMORY CURVES (t = 0:6) ───────────────────────────────────────────────

im_rows_nb = vector("list", length(P_list_nb))
names(im_rows_nb) = names(P_list_nb)

for (key in names(P_list_nb)) {
  rows = lapply(0:6, function(t) {
    vals = im_from_P(P_list_nb[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_nb[[key]] = do.call(rbind, rows)
}

im_df_nb = do.call(rbind, im_rows_nb)
im_df_nb$origin = factor(im_df_nb$origin, levels = rel_level_order)

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/nonblack", recursive = TRUE, showWarnings = FALSE)

# Transition matrix heatmaps
for (key in names(P_list_nb)) {
  p = make_combined(
    P_list_nb[[key]], pi0_list_nb[[key]], pistar_list_nb[[key]],
    levels    = rel_level_order,
    title_str = paste0("Non-Black Sample — Cohort ", key, "–",
                       as.integer(key) + 9, "  (N = ", n_list_nb[[key]], ")")
  )
  ggsave(paste0("output/figures/nonblack/trans_", key, "_10yr_nb.png"),
         p, width = 10, height = 7, dpi = 200)
}

# IM memory curves (all cohorts, faceted)
p_im_nb = ggplot(im_df_nb, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL,
       title = "Individual Memory by Cohort — Non-Black Sample (10-year bins, t = 0–6)") +
  healy_theme

ggsave("output/figures/nonblack/im_memory_10yr_nb.png",
       p_im_nb, width = 14, height = 5, dpi = 200)

# IM memory curves: one panel per origin state, cohorts as lines
p_im_nb_byorigin = ggplot(im_df_nb, aes(x = t, y = im, color = factor(cohort), group = cohort)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ origin, nrow = 1,
             labeller = labeller(origin = c(
               catholic    = "Catholic",
               evangelical = "Evangelical",
               mainline    = "Mainline",
               other       = "Other",
               none        = "None"
             ))) +
  scale_color_brewer(palette = "Dark2", name = "Birth cohort") +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       title = "Individual Memory by Origin — Non-Black Sample (cohorts 1940–1980)") +
  healy_theme

ggsave("output/figures/nonblack/im_memory_byorigin_nb.png",
       p_im_nb_byorigin, width = 14, height = 5, dpi = 200)

# ── FULL-SAMPLE MATRICES (same cohorts, all races) ───────────────────────────
# Built here to enable cell-by-cell difference with P_list_nb.

data_all = gss_all |>
  select(year, cohort, reltrad, reltrad16, race) |>
  filter(!(year %in% c(1972, 2021))) |>
  mutate(across(c(reltrad, reltrad16),
                ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  mutate(across(c(reltrad, reltrad16),
                ~ case_when(
                  . == "jewish"           ~ "other",
                  . == "black protestant" ~ "evangelical",
                  TRUE                    ~ .
                ),
                .names = "{.col}_alt")) |>
  mutate(age = year - cohort) |>
  filter(age > 25, cohort >= 1900) |>
  mutate(cohort_10 = floor((cohort - 1900) / 10) * 10 + 1900)

P_list_all    = list()
pi0_list_all  = list()
pistar_list_all = list()
n_list_all    = list()

for (coh in cohorts_nb) {
  sub = data_all[
    !is.na(data_all$cohort_10)     & data_all$cohort_10     == coh &
    !is.na(data_all$reltrad16_alt) & !is.na(data_all$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)
  P_list_all[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_all[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_all[[key]] = pi_star(P_list_all[[key]])
  n_list_all[[key]]      = nrow(sub)
}

# ── DIFFERENCE MATRICES (P_overall − P_nonblack) ─────────────────────────────
# Positive = overall > non-Black (Black respondents raise this cell's probability)
# Negative = non-Black > overall

diff_theme = theme_bw(base_size = 12) +
  theme(
    axis.ticks      = element_blank(),
    axis.title      = element_text(size = 12),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    legend.title    = element_text(size = 11),
    plot.title      = element_text(hjust = 0.5, size = 13, face = "plain"),
    panel.grid      = element_blank(),
    panel.border    = element_blank()
  )

plot_diff_heatmap = function(D, title_str, lim, text_size = 4.5) {
  df = as.data.frame(as.table(D))
  names(df) = c("origin", "current", "diff")
  df$origin  = factor(df$origin,  levels = rel_level_order)
  df$current = factor(df$current, levels = rev(rel_level_order))

  ggplot(df, aes(x = current, y = origin, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.3f", diff)), size = text_size) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-lim, lim), name = "Δ prob."
    ) +
    labs(x = "Current religion", y = "Origin religion", title = title_str) +
    diff_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11))
}

# Δ column (π₀ or π*): diverging tile, same color scale as the main matrix.
plot_diff_col = function(delta_vec, title_str, lim, text_size = 4.5) {
  df = data.frame(
    origin = factor(names(delta_vec), levels = rel_level_order),
    diff   = as.numeric(delta_vec)
  )
  ggplot(df, aes(x = 1, y = origin, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.3f", diff)), size = text_size, color = "black") +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-lim, lim), guide = "none"
    ) +
    labs(title = title_str) +
    diff_theme +
    theme(axis.title.x  = element_blank(), axis.text.x  = element_blank(),
          axis.title.y  = element_blank(), axis.text.y  = element_blank())
}

dir.create("output/figures/nonblack", recursive = TRUE, showWarnings = FALSE)

diff_rows = list()

for (key in names(P_list_nb)) {
  if (is.null(P_list_all[[key]])) next

  D        = P_list_all[[key]] - P_list_nb[[key]]
  delta_pi0   = pi0_list_all[[key]]   - pi0_list_nb[[key]]
  delta_pistar = pistar_list_all[[key]] - pistar_list_nb[[key]]

  # Shared symmetric limit so matrix and columns use the same scale
  lim = max(abs(c(as.numeric(D), delta_pi0, delta_pistar)), na.rm = TRUE)

  coh_yr = as.integer(key)
  n_all  = n_list_all[[key]]
  n_nb   = n_list_nb[[key]]

  g_mat  = plot_diff_heatmap(D,
    title_str = paste0("P overall − P non-Black\nCohort ", key, "–", coh_yr + 9,
                       "  (N all = ", n_all, "; N nb = ", n_nb, ")"),
    lim = lim)
  g_pi0  = plot_diff_col(delta_pi0,    title_str = "Δπ₀",  lim = lim)
  g_pist = plot_diff_col(delta_pistar, title_str = "Δπ*",  lim = lim)

  p_diff = patchwork::wrap_plots(g_mat, g_pi0, g_pist, widths = c(6, 1, 1))

  ggsave(paste0("output/figures/nonblack/diff_", key, "_10yr.png"),
         p_diff, width = 9, height = 6, dpi = 200)

  # Collect for faceted summary plot
  df_d = as.data.frame(as.table(D))
  names(df_d) = c("origin", "current", "diff")
  df_d$cohort = coh_yr
  diff_rows[[key]] = df_d
}

# ── FACETED SUMMARY: all cohorts in one figure ───────────────────────────────

diff_all = do.call(rbind, diff_rows)
diff_all$origin  = factor(diff_all$origin,  levels = rel_level_order)
diff_all$current = factor(diff_all$current, levels = rev(rel_level_order))
lim_global = max(abs(diff_all$diff), na.rm = TRUE)

p_diff_facet = ggplot(diff_all, aes(x = current, y = origin, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.3f", diff)), size = 3) +
  facet_wrap(~ cohort, nrow = 1,
             labeller = labeller(cohort = function(x) paste0(x, "–", as.integer(x) + 9))) +
  scale_fill_gradient2(
    low      = "#4393C3",
    mid      = "white",
    high     = "#D6604D",
    midpoint = 0,
    limits   = c(-lim_global, lim_global),
    name     = "Δ prob. (overall − non-Black)"
  ) +
  labs(x = "Current religion", y = "Origin religion",
       title = "Difference Matrices: P overall − P non-Black (10-year cohorts, 1940–1980)") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y     = element_text(size = 9),
    axis.ticks      = element_blank(),
    axis.title      = element_text(size = 11),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text      = element_text(size = 10),
    panel.grid      = element_blank(),
    panel.border    = element_blank()
  )

ggsave("output/figures/nonblack/diff_all_cohorts_10yr.png",
       p_diff_facet, width = 22, height = 5, dpi = 200)
