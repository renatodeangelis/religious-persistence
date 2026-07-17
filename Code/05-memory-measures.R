# ── 05 · NATIONAL MEMORY MEASURES & FIGURES ────────────────────────────────────
# Individual memory curves (IM), overall/exchange/structural mobility, mean time
# to exit (MTE), and the national transition-matrix heatmaps. Figure blocks use
# the pre-built national matrices from 02; the mobility and MTE time series
# recompute at 1-year cohort resolution from the cleaned data.
#
# Input:  data/derived/matrices.rds, data/derived/gss_clean.rds
# Output: output/figures/*.png

library(tidyr)
source("code/utils.R")

matrices   = readRDS("data/derived/matrices.rds")
clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

P_list_5      = matrices$nat5$P
pi0_list_5    = matrices$nat5$pi0
pistar_list_5 = matrices$nat5$pistar
P_list_10      = matrices$nat10$P
pi0_list_10    = matrices$nat10$pi0
pistar_list_10 = matrices$nat10$pistar
P_list_20      = matrices$nat20$P
pi0_list_20    = matrices$nat20$pi0
pistar_list_20 = matrices$nat20$pistar

# ── NATIONAL IM COMPUTATION ──────────────────────────────────────────────────

# ── IM LOOP (10-year cohorts, t = 0:4) ──────────────────────────────────────
im_rows_10 = vector("list", length(P_list_10))
names(im_rows_10) = names(P_list_10)

for (key in names(P_list_10)) {
  rows = lapply(0:6, function(t) {
    vals = im_from_P(P_list_10[[key]], t = t)
    data.frame(cohort = as.numeric(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_10[[key]] = do.call(rbind, rows)
}

im_df_10 = do.call(rbind, im_rows_10)

# ── IM LOOP (20-year cohorts, t = 0:4) ──────────────────────────────────────
im_rows_20 = vector("list", length(P_list_20))
names(im_rows_20) = names(P_list_20)

for (key in names(P_list_20)) {
  rows = lapply(0:4, function(t) {
    vals = im_from_P(P_list_20[[key]], t = t)
    data.frame(cohort = as.numeric(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_20[[key]] = do.call(rbind, rows)
}

im_df_20 = do.call(rbind, im_rows_20)

# ── NATIONAL MOBILITY ────────────────────────────────────────────────────────

mob_rows = lapply(1920:1980, function(coh) {
  sub = data[!is.na(data$cohort) & data$cohort == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) return(NULL)
  P   = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0 = pi_0(sub, "reltrad16_alt")
  data.frame(cohort = coh, mobility = overall_mobility(P, pi0))
})
mob_df = do.call(rbind, Filter(Negate(is.null), mob_rows))

# ── EXCHANGE/STRUCTURAL MOBILITY (1-year cohorts, t = 0) ─────────────────────

em_sm_rows = lapply(1920:1985, function(coh) {
  sub = data[!is.na(data$cohort) & data$cohort == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) return(NULL)
  P   = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0 = pi_0(sub, "reltrad16_alt")
  om_v = overall_mobility(P, pi0)
  sm_v = sm(P, pi0, t = 0)
  data.frame(cohort = coh, EM = om_v - sm_v, SM = sm_v, OM = om_v)
})
em_sm_df = do.call(rbind, Filter(Negate(is.null), em_sm_rows))

# ── NATIONAL FIGURES ──────────────────────────────────────────────────────────

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# Keys are bin midpoints (edge + halfwidth); recover the integer left edge for
# titles/filenames so each figure still labels its actual cohort range.
for (key in names(P_list_5)) {
  edge = as.numeric(key) - 2.5
  p = make_combined(P_list_5[[key]], pi0_list_5[[key]], pistar_list_5[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", edge, "–", edge + 4))
  ggsave(paste0("output/figures/trans_", edge, "_5yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

for (key in names(P_list_10)) {
  edge = as.numeric(key) - 5
  p = make_combined(P_list_10[[key]], pi0_list_10[[key]], pistar_list_10[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", edge, "–", edge + 9))
  ggsave(paste0("output/figures/trans_", edge, "_10yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

for (key in names(P_list_20)) {
  edge = as.numeric(key) - 10
  p = make_combined(P_list_20[[key]], pi0_list_20[[key]], pistar_list_20[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", edge, "–", edge + 19))
  ggsave(paste0("output/figures/trans_", edge, "_20yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

# Okabe-Ito palette (Healy) — lowercase keys match rel_level_order
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

im_df_10 = im_df_10[im_df_10$cohort != 1925, ]   # drop earliest 10-yr bin (edge 1920, midpoint 1925)
im_df_10$origin = factor(im_df_10$origin, levels = rel_level_order)

p_im = ggplot(im_df_10, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL, title = "Individual Memory by Cohort (10-year bins, t = 0–4)") +
  healy_theme

ggsave("output/figures/im_memory_10yr.png", p_im, width = 8, height = 6, dpi = 200)

im_df_20$origin = factor(im_df_20$origin, levels = rel_level_order)

p_im_20 = ggplot(im_df_20, aes(x = t, y = im, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL, title = "Individual Memory by Cohort (20-year bins, t = 0–4)") +
  healy_theme

ggsave("output/figures/im_memory_20yr.png", p_im_20, width = 10, height = 5, dpi = 200)

p_mob = ggplot(mob_df[mob_df$cohort >= 1930 & mob_df$cohort <= 1985, ], aes(x = cohort, y = mobility)) +
  geom_point(size = 1.5, alpha = 0.6, color = "#0072B2") +
  geom_smooth(method = "loess", se = TRUE, span = 0.4, alpha = 0.2,
              color = "#0072B2", fill = "#0072B2") +
  scale_y_continuous(limits = c(0.15, 0.45)) +
  labs(x = "Birth cohort", y = "Probability to Move",
       title = "Overall Mobility by Birth Cohort") +
  healy_theme

ggsave("output/figures/overall_mobility.png", p_mob, width = 8, height = 5, dpi = 200)

em_sm_long = pivot_longer(
  em_sm_df[em_sm_df$cohort >= 1930 & em_sm_df$cohort <= 1985, ],
  cols = c(EM, SM),
  names_to = "measure", values_to = "value"
)

p_em_sm = ggplot(em_sm_long, aes(x = cohort, y = value, color = measure,
                                  fill = measure, group = measure)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, span = 0.4, alpha = 0.2) +
  scale_color_manual(name   = NULL,
                     values = c(EM = "#D55E00", SM = "#009E73"),
                     labels = c(EM = "Exchange mobility", SM = "Structural mobility")) +
  scale_fill_manual(name   = NULL,
                    values = c(EM = "#D55E00", SM = "#009E73"),
                    labels = c(EM = "Exchange mobility", SM = "Structural mobility")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Birth cohort", y = "Probability to Move",
       title = "Exchange and Structural Mobility by Birth Cohort") +
  healy_theme

ggsave("output/figures/em_sm_pooled.png", p_em_sm, width = 8, height = 5, dpi = 200)

# ── MTE TIME SERIES (5-year cohorts, 1930–1980) ───────────────────────────────

mte_rows_5 = lapply(names(P_list_5), function(key) {
  coh = as.numeric(key)                          # 5-yr bin midpoint (e.g. 1932.5)
  if (coh < 1932.5 | coh > 1982.5) return(NULL)  # edges 1930–1980
  vals = mte(P_list_5[[key]])
  data.frame(cohort = coh, origin = names(vals), mte = vals, row.names = NULL)
})
mte_df_5 = do.call(rbind, Filter(Negate(is.null), mte_rows_5))
mte_df_5$origin = factor(mte_df_5$origin, levels = rel_level_order)

p_mte = ggplot(mte_df_5, aes(x = cohort, y = mte, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = seq(1930, 1980, by = 10)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Birth cohort (5-year bins)", y = "Mean time to exit (steps)",
       color = NULL,
       title = "Mean Time to Exit by Religious Origin (5-year cohorts, 1930–1980)") +
  healy_theme

ggsave("output/figures/mte_5yr.png", p_mte, width = 8, height = 5, dpi = 200)

# ── MTE TIME SERIES (1-year cohorts, 1930–1980) ───────────────────────────────

mte_rows_1 = lapply(1930:1980, function(coh) {
  sub = data[!is.na(data$cohort) & data$cohort == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) return(NULL)
  P    = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  vals = mte(P)
  data.frame(cohort = coh, origin = names(vals), mte = vals, row.names = NULL)
})
mte_df_1 = do.call(rbind, Filter(Negate(is.null), mte_rows_1))
mte_df_1$origin = factor(mte_df_1$origin, levels = rel_level_order)

p_mte_1 = ggplot(mte_df_1, aes(x = cohort, y = mte, color = origin, group = origin)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, span = 0.4) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = seq(1930, 1980, by = 10)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Birth cohort", y = "Mean time to exit (steps)",
       color = NULL,
       title = "Mean Time to Exit by Religious Origin (1-year cohorts, 1930–1980)") +
  healy_theme

ggsave("output/figures/mte_1yr.png", p_mte_1, width = 8, height = 5, dpi = 200)
