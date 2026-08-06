# ── 05 · NATIONAL MEMORY MEASURES & FIGURES ────────────────────────────────────
# Individual memory curves (IM), overall mobility, mean time to exit (MTE), and
# the national transition-matrix heatmaps. Figure blocks use the pre-built
# national matrices from 02; the mobility and MTE time series recompute at
# 1-year cohort resolution from the cleaned data.
#
# Input:  data/derived/matrices.rds, data/derived/gss_clean.rds
# Output: output/figures/*.png

source("code/utils.R")

matrices   = readRDS("data/derived/matrices.rds")
clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

P_list_10      = matrices$nat10$P
pi0_list_10    = matrices$nat10$pi0
pistar_list_10 = matrices$nat10$pistar

# ── NATIONAL IM COMPUTATION ──────────────────────────────────────────────────

# ── IM LOOP (10-year cohorts, t = 0:4) ──────────────────────────────────────
im_rows_10 = vector("list", length(P_list_10))
names(im_rows_10) = names(P_list_10)

for (key in names(P_list_10)) {
  rows = lapply(0:4, function(t) {
    vals = im_from_P(P_list_10[[key]], t = t)
    data.frame(cohort = as.numeric(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_10[[key]] = do.call(rbind, rows)
}

im_df_10 = do.call(rbind, im_rows_10)


# ── NATIONAL FIGURES ──────────────────────────────────────────────────────────

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_10)) {
  edge = as.numeric(key) - 5
  p = make_combined(P_list_10[[key]], pi0_list_10[[key]], pistar_list_10[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", edge, "–", sprintf("%02d", (edge + 9) %% 100)))
  ggsave(paste0("output/figures/trans_", key, "_10yr.png"), p,
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

im_df_10$origin = factor(im_df_10$origin, levels = rel_level_order)

mids_im = c(1930, 1940, 1950, 1960, 1970, 1980)
im_df_10$cohort_label = factor(
  paste0(im_df_10$cohort - 5, "–", sprintf("%02d", (im_df_10$cohort + 4) %% 100)),
  levels = paste0(mids_im - 5, "–", sprintf("%02d", (mids_im + 4) %% 100))
)

p_im = ggplot(im_df_10 |> filter(cohort %in% mids_im), aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ cohort_label, nrow = 2) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL, title = "Individual Memory by Birth Cohort (t = 0–4)") +
  healy_theme

ggsave("output/figures/im_memory_10yr.png", p_im, width = 8, height = 6, dpi = 200)


# ── MTE TIME SERIES (10-year cohorts, 1925–1984) ─────────────────────────────

mte_rows_10 = lapply(names(P_list_10), function(key) {
  vals = mte(P_list_10[[key]])
  data.frame(cohort = as.numeric(key), origin = names(vals), mte = vals, row.names = NULL)
})
mte_df_10 = do.call(rbind, Filter(Negate(is.null), mte_rows_10))
mte_df_10$origin = factor(mte_df_10$origin, levels = rel_level_order)

p_mte = ggplot(mte_df_10, aes(x = cohort, y = mte, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = mids_im,
                     labels = paste0(mids_im - 5, "–", sprintf("%02d", (mids_im + 4) %% 100))) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Birth cohort (10-year bins)", y = "Mean time to exit (steps)",
       color = NULL,
       title = "Mean Time to Exit by Religious Origin (10-year cohorts, 1925–1984)") +
  healy_theme

ggsave("output/figures/mte_10yr.png", p_mte, width = 8, height = 5, dpi = 200)

# ── SHANNON ENTROPY BY COHORT (10-year bins) ─────────────────────────────────
# E(π₀): entropy of the childhood religion distribution (origin)
# E(π₁): entropy of the adult religion distribution (current); π₁ = π₀ P

entropy_df = do.call(rbind, lapply(names(P_list_10), function(key) {
  pi0 = pi0_list_10[[key]]
  pi1 = as.numeric(pi0 %*% P_list_10[[key]])
  names(pi1) = names(pi0)
  data.frame(cohort = as.numeric(key),
             e_pi0  = shannon_entropy(pi0),
             e_pi1  = shannon_entropy(pi1))
}))

entropy_long = rbind(
  data.frame(cohort = entropy_df$cohort, distribution = "Origin (π₀)", entropy = entropy_df$e_pi0),
  data.frame(cohort = entropy_df$cohort, distribution = "Current (π₁)", entropy = entropy_df$e_pi1)
)
entropy_long$distribution = factor(entropy_long$distribution,
                                   levels = c("Origin (π₀)", "Current (π₁)"))

p_entropy = ggplot(entropy_long, aes(x = cohort, y = entropy,
                                     color = distribution, shape = distribution,
                                     group = distribution)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Origin (π₀)" = "#5C4A8F", "Current (π₁)" = "#2D6A4F"),
                     name = NULL) +
  scale_shape_manual(values = c("Origin (π₀)" = 16, "Current (π₁)" = 15),
                     name = NULL) +
  scale_x_continuous(breaks = seq(1930, 1990, by = 10)) +
  labs(x = "Birth cohort (10-year bin midpoint)",
       y = expression(paste("Shannon entropy  ", E(mu) == -Sigma, mu[i], ln(mu[i]))),
       title = "Shannon Entropy of Origin and Current Religious Distributions by Cohort") +
  healy_theme

ggsave("output/figures/shannon_entropy_10yr.png", p_entropy, width = 8, height = 5, dpi = 200)

# ── TRANSITION MATRIX GRID (10-year cohorts, 1935–1984) ──────────────────────
# Five panels in one figure: one column per birth cohort. Rows = origin (RELIG16),
# columns = current (RELIG). Same layout as P_grid_U_vs_S from 08 but a single
# row — gives the reader a cross-cohort view of the full matrix at a glance.

mids_grid = c(1930, 1940, 1950, 1960, 1970, 1980)   # cohort_10 midpoints → edges 1925–1984

n_per_mid = sapply(mids_grid, function(m)
  format(sum(matrices$nat10$N[[as.character(m)]]), big.mark = ","))

grid_df = do.call(rbind, lapply(seq_along(mids_grid), function(i) {
  mid = mids_grid[i]
  key = as.character(mid)
  P   = matrices$nat10$P[[key]]
  if (is.null(P)) return(NULL)
  df  = as.data.frame(as.table(P))
  names(df) = c("origin", "current", "p")
  df$cohort = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100), "\n(N = ", n_per_mid[i], ")")
  df
}))

grid_df$origin  = factor(grid_df$origin,  levels = rel_level_order)
grid_df$current = factor(grid_df$current, levels = rev(rel_level_order))
grid_df$cohort  = factor(grid_df$cohort,
  levels = paste0(mids_grid - 5, "–", sprintf("%02d", (mids_grid + 4) %% 100), "\n(N = ", n_per_mid, ")"))

p_grid_national = ggplot(grid_df, aes(current, origin, fill = p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p)), size = 2.8) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1), name = "P[i→j]") +
  labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
       title = "Intergenerational transition matrices by birth cohort") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 9),
    plot.caption     = element_text(size = 8, color = "grey50")
  )

ggsave("output/figures/P_grid_national_10yr.png", p_grid_national,
       width = 12, height = 9, dpi = 200)

# ── π₀ AND π* DISTRIBUTION GRID (10-year cohorts, 1935–1994) ─────────────────
# Two companions to the matrix grid above. The grid figure (facet_grid) shows
# the bar-chart layout cohort-by-cohort; the line figure shows the trend across
# cohorts for each religion — more useful for reading the narrative.

dist_df = do.call(rbind, lapply(mids_grid, function(mid) {
  key    = as.character(mid)
  pi0    = matrices$nat10$pi0[[key]]
  pistar = matrices$nat10$pistar[[key]]
  if (is.null(pi0)) return(NULL)
  cohort_lbl = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100))
  rbind(
    data.frame(mid = mid, cohort = cohort_lbl, measure = "π₀  (origin)",
               religion = names(pi0),    value = as.numeric(pi0)),
    data.frame(mid = mid, cohort = cohort_lbl, measure = "π* (stationary)",
               religion = names(pistar), value = as.numeric(pistar))
  )
}))

dist_df$religion = factor(dist_df$religion, levels = rel_level_order)
dist_df$cohort   = factor(dist_df$cohort,
                           levels = paste0(mids_grid - 5, "–", sprintf("%02d", (mids_grid + 4) %% 100)))
dist_df$measure  = factor(dist_df$measure,
                           levels = c("π₀  (origin)", "π* (stationary)"))

# Grid version: facet_grid(measure × cohort), one bar chart per cell.
# coord_flip() keeps religion labels horizontal and readable at small size.
p_dist_grid = ggplot(dist_df, aes(y = religion, x = value, fill = religion)) +
  geom_col(width = 0.65) +
  facet_grid(measure ~ cohort) +
  scale_fill_manual(values = reltrad_colors, labels = reltrad_labels_tc, name = NULL) +
  scale_x_continuous(limits = c(0, 0.65), breaks = c(0, 0.25, 0.5),
                     labels = c("0", ".25", ".5")) +
  scale_y_discrete(labels = reltrad_labels_tc) +
  labs(x = "Share", y = NULL,
       title = "Origin (π₀) and stationary (π*) distributions by birth cohort",
) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position   = "none",
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(size = 8),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 7),
    plot.caption      = element_text(size = 8, color = "grey50")
  )

ggsave("output/figures/pi_dist_grid_10yr.png", p_dist_grid,
       width = 14, height = 5, dpi = 200)
