# ── 10 · ROBUSTNESS: NON-BLACK SAMPLE ─────────────────────────────────────────
# Re-estimates the national 10-year cohort matrices on the non-Black sample
# (race != 2) and plots cell-by-cell difference heatmaps against the full-sample
# matrices (P overall − P non-Black). Both samples share the derived pipeline
# frame, so the only thing that varies is the race restriction.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_nonblack.rds
#         output/figures/nonblack/*.png

library(dplyr)
library(ggplot2)
library(patchwork)
source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")
reltrad_colors = c(
  catholic    = "#0072B2",
  evangelical = "#D55E00",
  mainline    = "#009E73",
  other       = "#CC79A7",
  none        = "#999999"
)
reltrad_labels_tc = c(
  catholic = "Catholic", evangelical = "Evangelical", mainline = "Mainline",
  other = "Other", none = "None"
)

# 10-year cohort midpoints (edges 1940–1980). cohort_10 in the pipeline is the
# bin midpoint (edge + 5); the edge is used for titles and filenames.
mids_10 = c(1940, 1950, 1960, 1970, 1980)

data_nb = data[as.numeric(data$race) != 2, ]

# ── HELPER: build P/pi0/pistar lists for a given data frame ───────────────────

build_lists = function(df) {
  P = pi0 = pistar = nn = list()
  for (mid in mids_10) {
    sub = df[!is.na(df$cohort_10) & df$cohort_10 == mid &
             !is.na(df$reltrad16_alt) & !is.na(df$reltrad_alt), ]
    if (nrow(sub) < 30) next
    key = as.character(mid - 5)   # edge label
    P[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0[[key]]    = pi_0(sub, "reltrad16_alt")
    pistar[[key]] = pi_star(P[[key]])
    nn[[key]]     = nrow(sub)
  }
  list(P = P, pi0 = pi0, pistar = pistar, n = nn)
}

nb  = build_lists(data_nb)   # non-Black
all = build_lists(data)      # full sample (same cohorts, all races)

# Print non-Black matrices
for (key in names(nb$P)) {
  cat("\n── Non-Black cohort", key, "–", as.integer(key) + 9,
      "  (N =", nb$n[[key]], ") ──\n")
  print(round(nb$P[[key]], 3))
}

# ── IM MEMORY CURVES (non-Black, t = 0:6) ────────────────────────────────────

im_rows = lapply(names(nb$P), function(key) {
  do.call(rbind, lapply(0:6, function(t) {
    vals = im_from_P(nb$P[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  }))
})
im_df_nb = do.call(rbind, im_rows)
im_df_nb$origin = factor(im_df_nb$origin, levels = rel_level_order)

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/nonblack", recursive = TRUE, showWarnings = FALSE)

# Transition-matrix heatmaps
for (key in names(nb$P)) {
  p = make_combined(
    nb$P[[key]], nb$pi0[[key]], nb$pistar[[key]],
    levels    = rel_level_order,
    title_str = paste0("Non-Black Sample — Cohort ", key, "–",
                       as.integer(key) + 9, "  (N = ", nb$n[[key]], ")")
  )
  ggsave(paste0("output/figures/nonblack/trans_", key, "_10yr_nb.png"),
         p, width = 10, height = 7, dpi = 200)
}

# IM memory curves (all cohorts, faceted by cohort)
p_im_nb = ggplot(im_df_nb, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)", color = NULL,
       title = "Individual Memory by Cohort — Non-Black Sample (10-year bins, t = 0–6)") +
  healy_theme

ggsave("output/figures/nonblack/im_memory_10yr_nb.png",
       p_im_nb, width = 14, height = 5, dpi = 200)

# IM memory curves: one panel per origin, cohorts as lines
p_im_nb_byorigin = ggplot(im_df_nb, aes(x = t, y = im, color = factor(cohort), group = cohort)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ origin, nrow = 1, labeller = labeller(origin = reltrad_labels_tc)) +
  scale_color_brewer(palette = "Dark2", name = "Birth cohort") +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       title = "Individual Memory by Origin — Non-Black Sample (cohorts 1940–1980)") +
  healy_theme

ggsave("output/figures/nonblack/im_memory_byorigin_nb.png",
       p_im_nb_byorigin, width = 14, height = 5, dpi = 200)

# ── DIFFERENCE MATRICES (P overall − P non-Black) ────────────────────────────
# Positive = overall > non-Black (Black respondents raise this cell's probability)

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
    scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                         midpoint = 0, limits = c(-lim, lim), name = "Δ prob.") +
    labs(x = "Current religion", y = "Origin religion", title = title_str) +
    diff_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11))
}

plot_diff_col = function(delta_vec, title_str, lim, text_size = 4.5) {
  df = data.frame(origin = factor(names(delta_vec), levels = rel_level_order),
                  diff   = as.numeric(delta_vec))
  ggplot(df, aes(x = 1, y = origin, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.3f", diff)), size = text_size, color = "black") +
    scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                         midpoint = 0, limits = c(-lim, lim), guide = "none") +
    labs(title = title_str) +
    diff_theme +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank())
}

diff_rows = list()

for (key in names(nb$P)) {
  if (is.null(all$P[[key]])) next

  D            = all$P[[key]]      - nb$P[[key]]
  delta_pi0    = all$pi0[[key]]    - nb$pi0[[key]]
  delta_pistar = all$pistar[[key]] - nb$pistar[[key]]

  lim = max(abs(c(as.numeric(D), delta_pi0, delta_pistar)), na.rm = TRUE)

  g_mat  = plot_diff_heatmap(D,
    title_str = paste0("P overall − P non-Black\nCohort ", key, "–", as.integer(key) + 9,
                       "  (N all = ", all$n[[key]], "; N nb = ", nb$n[[key]], ")"),
    lim = lim)
  g_pi0  = plot_diff_col(delta_pi0,    title_str = "Δπ₀", lim = lim)
  g_pist = plot_diff_col(delta_pistar, title_str = "Δπ*", lim = lim)

  p_diff = patchwork::wrap_plots(g_mat, g_pi0, g_pist, widths = c(6, 1, 1))
  ggsave(paste0("output/figures/nonblack/diff_", key, "_10yr.png"),
         p_diff, width = 9, height = 6, dpi = 200)

  df_d = as.data.frame(as.table(D))
  names(df_d) = c("origin", "current", "diff")
  df_d$cohort = as.integer(key)
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
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-lim_global, lim_global),
                       name = "Δ prob. (overall − non-Black)") +
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

saveRDS(
  list(
    nonblack = list(P = nb$P,  pi0 = nb$pi0,  pistar = nb$pistar,  n = nb$n),
    full      = list(P = all$P, pi0 = all$pi0, pistar = all$pistar, n = all$n)
  ),
  "data/derived/matrices_nonblack.rds"
)
cat("\nDone. Non-Black robustness figures in output/figures/nonblack/.\n")
cat("Wrote data/derived/matrices_nonblack.rds\n")
