library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gssr)
source("Code/utils.R")

reltrad_labels = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

# 6-state order: black protestant retained, jewish → other
rel_level_order_6 = c("catholic", "evangelical", "black protestant", "mainline", "other", "none")

reltrad_colors_6 = c(
  catholic           = "#0072B2",
  evangelical        = "#D55E00",
  "black protestant" = "#E69F00",
  mainline           = "#009E73",
  other              = "#CC79A7",
  none               = "#999999"
)
reltrad_labels_6 = c(
  catholic           = "Catholic",
  evangelical        = "Evangelical",
  "black protestant" = "Black Protestant",
  mainline           = "Mainline",
  other              = "Other",
  none               = "None"
)

# ── DATA (full sample, 6-state space) ────────────────────────────────────────

data(gss_all)
data_bp = gss_all |>
  select(year, cohort, reltrad, reltrad16) |>
  filter(!(year %in% c(1972, 2021))) |>
  mutate(across(c(reltrad, reltrad16),
                ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  # Collapse only jewish → other; black protestant stays
  mutate(across(c(reltrad, reltrad16),
                ~ if_else(. == "jewish", "other", .),
                .names = "{.col}_6")) |>
  mutate(age = year - cohort) |>
  filter(age > 25, cohort >= 1900) |>
  mutate(
    cohort_10 = floor((cohort - 1900) / 10) * 10 + 1900,
    cohort_20 = floor((cohort - 1900) / 20) * 20 + 1900
  )

states_6 = rel_level_order_6  # fix ordering throughout

# ── CELL-COUNT DIAGNOSTIC ─────────────────────────────────────────────────────
# Black protestant rows are likely thin in early cohorts — inspect before trusting matrices.

cell_diag = lapply(sort(unique(data_bp$cohort_20[data_bp$cohort_20 >= 1920 &
                                                   data_bp$cohort_20 <= 1980])), function(coh) {
  sub = data_bp[!is.na(data_bp$cohort_20) & data_bp$cohort_20 == coh &
                !is.na(data_bp$reltrad16_6) & !is.na(data_bp$reltrad_6), ]
  tab = table(
    factor(sub$reltrad16_6, levels = states_6),
    factor(sub$reltrad_6,   levels = states_6)
  )
  data.frame(
    cohort    = coh,
    n         = nrow(sub),
    n_bp_orig = sum(tab["black protestant", ]),
    min_cell  = min(tab),
    cells_lt5 = sum(tab < 5),
    cells_0   = sum(tab == 0),
    row.names = NULL
  )
}) |> do.call(what = rbind)

cat("\n── Cell-count diagnostic (6-state, 20-year cohorts) ──\n")
print(cell_diag, row.names = FALSE)

# ── 20-YEAR COHORT TRANSITION MATRICES ───────────────────────────────────────

cohorts_20 = sort(unique(data_bp$cohort_20[!is.na(data_bp$cohort_20) &
                                             data_bp$cohort_20 >= 1920 &
                                             data_bp$cohort_20 <= 1980]))

P_list_bp      = list()
pi0_list_bp    = list()
pistar_list_bp = list()
n_list_bp      = list()

for (coh in cohorts_20) {
  sub = data_bp[!is.na(data_bp$cohort_20) & data_bp$cohort_20 == coh &
                !is.na(data_bp$reltrad16_6) & !is.na(data_bp$reltrad_6), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_bp[[key]]      = p_matrix(sub, "reltrad16_6", "reltrad_6", levels = states_6)
  pi0_list_bp[[key]]    = pi_0(sub, "reltrad16_6")
  pistar_list_bp[[key]] = pi_star(P_list_bp[[key]])
  n_list_bp[[key]]      = nrow(sub)
}

# Print matrices
for (key in names(P_list_bp)) {
  cat("\n── Cohort", key, "–", as.integer(key) + 19,
      "  (N =", n_list_bp[[key]], ") ──\n")
  print(round(P_list_bp[[key]], 3))
}

# ── IM MEMORY CURVES (t = 0:6) ───────────────────────────────────────────────

im_rows_bp = vector("list", length(P_list_bp))
names(im_rows_bp) = names(P_list_bp)

for (key in names(P_list_bp)) {
  rows = lapply(0:6, function(t) {
    vals = im_from_P(P_list_bp[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_bp[[key]] = do.call(rbind, rows)
}

im_df_bp = do.call(rbind, im_rows_bp)
im_df_bp$origin = factor(im_df_bp$origin, levels = rel_level_order_6)

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/bp", recursive = TRUE, showWarnings = FALSE)

# Transition matrix heatmaps (π₀ and π* columns included via make_combined)
for (key in names(P_list_bp)) {
  p = make_combined(
    P_list_bp[[key]], pi0_list_bp[[key]], pistar_list_bp[[key]],
    levels    = rel_level_order_6,
    title_str = paste0("6-State — Cohort ", key, "–", as.integer(key) + 19,
                       "  (N = ", n_list_bp[[key]], ")")
  )
  ggsave(paste0("output/figures/bp/trans_", key, "_20yr_6state.png"),
         p, width = 11, height = 7, dpi = 200)
}

# IM memory curves: faceted by cohort, colored by origin
p_im_bp = ggplot(im_df_bp, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ cohort, nrow = 1,
             labeller = labeller(cohort = function(x) paste0(x, "–", as.integer(x) + 19))) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL,
       title = "Individual Memory by Cohort — 6-State Space (20-year bins, t = 0–6)") +
  healy_theme

ggsave("output/figures/bp/im_memory_20yr_6state.png",
       p_im_bp, width = 14, height = 5, dpi = 200)

# IM memory curves: faceted by origin, cohorts as lines
p_im_bp_byorigin = ggplot(im_df_bp,
    aes(x = t, y = im, color = factor(cohort), group = cohort)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ origin, nrow = 1,
             labeller = labeller(origin = reltrad_labels_6)) +
  scale_color_brewer(palette = "Dark2", name = "Birth cohort") +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       title = "Individual Memory by Origin — 6-State Space (cohorts 1920–1980)") +
  healy_theme

ggsave("output/figures/bp/im_memory_byorigin_20yr_6state.png",
       p_im_bp_byorigin, width = 16, height = 5, dpi = 200)

# Diagonal persistence over cohorts
diag_df = do.call(rbind, lapply(names(P_list_bp), function(key) {
  P = P_list_bp[[key]]
  data.frame(cohort = as.integer(key), origin = rownames(P),
             persistence = diag(P), row.names = NULL)
}))
diag_df$origin = factor(diag_df$origin, levels = rel_level_order_6)

p_diag = ggplot(diag_df, aes(x = cohort, y = persistence, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = seq(1920, 1980, by = 20)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Birth cohort (20-year bins)", y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by Origin — 6-State Space (20-year cohorts)") +
  healy_theme

ggsave("output/figures/bp/diagonal_persistence_20yr_6state.png",
       p_diag, width = 8, height = 5, dpi = 200)

# ── 10-YEAR COHORTS, 1940–1989 (6-state) ─────────────────────────────────────

cohorts_10_bp = c(1940, 1950, 1960, 1970, 1980)

# Cell-count diagnostic — BP row likely thin for some decades
cell_diag_10 = lapply(cohorts_10_bp, function(coh) {
  sub = data_bp[!is.na(data_bp$cohort_10) & data_bp$cohort_10 == coh &
                !is.na(data_bp$reltrad16_6) & !is.na(data_bp$reltrad_6), ]
  tab = table(
    factor(sub$reltrad16_6, levels = states_6),
    factor(sub$reltrad_6,   levels = states_6)
  )
  data.frame(
    cohort    = coh,
    n         = nrow(sub),
    n_bp_orig = sum(tab["black protestant", ]),
    min_cell  = min(tab),
    cells_lt5 = sum(tab < 5),
    cells_0   = sum(tab == 0),
    row.names = NULL
  )
}) |> do.call(what = rbind)

cat("\n── Cell-count diagnostic (6-state, 10-year cohorts, 1940–1989) ──\n")
print(cell_diag_10, row.names = FALSE)

P_list_bp10      = list()
pi0_list_bp10    = list()
pistar_list_bp10 = list()
n_list_bp10      = list()

for (coh in cohorts_10_bp) {
  sub = data_bp[!is.na(data_bp$cohort_10) & data_bp$cohort_10 == coh &
                !is.na(data_bp$reltrad16_6) & !is.na(data_bp$reltrad_6), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_bp10[[key]]      = p_matrix(sub, "reltrad16_6", "reltrad_6", levels = states_6)
  pi0_list_bp10[[key]]    = pi_0(sub, "reltrad16_6")
  pistar_list_bp10[[key]] = pi_star(P_list_bp10[[key]])
  n_list_bp10[[key]]      = nrow(sub)
}

# Print matrices
for (key in names(P_list_bp10)) {
  cat("\n── Cohort", key, "–", as.integer(key) + 9,
      "  (N =", n_list_bp10[[key]], ") ──\n")
  print(round(P_list_bp10[[key]], 3))
}

# IM memory curves
im_rows_bp10 = vector("list", length(P_list_bp10))
names(im_rows_bp10) = names(P_list_bp10)

for (key in names(P_list_bp10)) {
  rows = lapply(0:6, function(t) {
    vals = im_from_P(P_list_bp10[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_bp10[[key]] = do.call(rbind, rows)
}

im_df_bp10 = do.call(rbind, im_rows_bp10)
im_df_bp10$origin = factor(im_df_bp10$origin, levels = rel_level_order_6)

# Heatmaps
for (key in names(P_list_bp10)) {
  p = make_combined(
    P_list_bp10[[key]], pi0_list_bp10[[key]], pistar_list_bp10[[key]],
    levels    = rel_level_order_6,
    title_str = paste0("6-State — Cohort ", key, "–", as.integer(key) + 9,
                       "  (N = ", n_list_bp10[[key]], ")")
  )
  ggsave(paste0("output/figures/bp/trans_", key, "_10yr_6state.png"),
         p, width = 11, height = 7, dpi = 200)
}

# IM curves faceted by cohort
p_im_bp10 = ggplot(im_df_bp10, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ cohort, nrow = 1,
             labeller = labeller(cohort = function(x) paste0(x, "–", as.integer(x) + 9))) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL,
       title = "Individual Memory by Cohort — 6-State Space (10-year bins, 1940–1989)") +
  healy_theme

ggsave("output/figures/bp/im_memory_10yr_6state.png",
       p_im_bp10, width = 14, height = 5, dpi = 200)

# Diagonal persistence over 10-year cohorts
diag_df10 = do.call(rbind, lapply(names(P_list_bp10), function(key) {
  P = P_list_bp10[[key]]
  data.frame(cohort = as.integer(key), origin = rownames(P),
             persistence = diag(P), row.names = NULL)
}))
diag_df10$origin = factor(diag_df10$origin, levels = rel_level_order_6)

p_diag10 = ggplot(diag_df10, aes(x = cohort, y = persistence, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = cohorts_10_bp) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Birth cohort (10-year bins)", y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by Origin — 6-State Space (10-year cohorts, 1940–1989)") +
  healy_theme

ggsave("output/figures/bp/diagonal_persistence_10yr_6state.png",
       p_diag10, width = 8, height = 5, dpi = 200)
