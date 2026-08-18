# ── 22 · ROBUSTNESS: SEVEN-STATE RELTRAD SCHEME (2022 & 2024 DROPPED) ─────────
# Identical to 21-robustness-7state.R except that GSS survey years 2022 and
# 2024 are excluded before estimation — a sensitivity check for the possibility
# that the two most-recent panels drive the youngest-cohort estimates.
# Also produces a difference grid (P_full − P_drop2224) as Figure 4.
#
# Input:  data/derived/gss_clean.rds
# Output: output/figures/7state-drop2224/*.png

library(dplyr)
library(ggplot2)
source("code/utils.R")

clean     = readRDS("data/derived/gss_clean.rds")
data_full = clean$data

# Filtered sample: drop survey years 2022 and 2024
data = data_full[!(data_full$year %in% c(2022, 2024)), ]

# ── STATE SPACE ───────────────────────────────────────────────────────────────

states_7 = sort(unique(c(data$reltrad, data$reltrad16)))
states_7  = states_7[!is.na(states_7)]

rel_level_order_7 = c("catholic", "evangelical", "black protestant",
                       "mainline", "jewish", "other", "none")

reltrad_colors_7 = c(
  "catholic"          = "#0072B2",
  "evangelical"       = "#D55E00",
  "black protestant"  = "#E69F00",
  "mainline"          = "#009E73",
  "jewish"            = "#56B4E9",
  "other"             = "#CC79A7",
  "none"              = "#999999"
)

reltrad_labels_7 = c(
  "catholic"          = "Catholic",
  "evangelical"       = "Evangelical",
  "black protestant"  = "Black Prot.",
  "mainline"          = "Mainline",
  "jewish"            = "Jewish",
  "other"             = "Other",
  "none"              = "None"
)

mids_10 = c(1930, 1940, 1950, 1960, 1970, 1980)

# ── THIN-CELL CHECK ───────────────────────────────────────────────────────────

thin_check = do.call(rbind, lapply(mids_10, function(mid) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == mid &
             !is.na(data$reltrad16), ]
  tab = table(sub$reltrad16)
  data.frame(cohort = mid, state = names(tab), n = as.numeric(tab),
             row.names = NULL)
}))

thin_at_risk = thin_check[
  thin_check$state %in% c("jewish", "black protestant") & thin_check$n < 30, ]

cat("Thin origin cells (N < 30) for Jewish or Black Protestant:\n")
if (nrow(thin_at_risk) > 0) print(thin_at_risk) else cat("  None found.\n")

thin_cohorts = unique(thin_at_risk$cohort)

# ── MATRIX ESTIMATION ─────────────────────────────────────────────────────────

P_list_7      = list()
pi0_list_7    = list()
pistar_list_7 = list()
N_list_7      = list()

for (coh in mids_10) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == coh &
             !is.na(data$reltrad16) & !is.na(data$reltrad), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_7[[key]]      = suppressWarnings(
    p_matrix(sub, "reltrad16", "reltrad", levels = states_7)
  )
  pi0_list_7[[key]]    = pi_0(sub, "reltrad16")
  pistar_list_7[[key]] = pi_star(P_list_7[[key]])
  N_list_7[[key]]      = suppressWarnings(
    count_matrix(sub, "reltrad16", "reltrad", levels = states_7)
  )
}

dir.create("output/figures/7state-drop2224", recursive = TRUE, showWarnings = FALSE)

# ── FIGURE 1: TRANSITION MATRIX GRID ─────────────────────────────────────────

cohort_labels_7 = sapply(seq_along(mids_10), function(i) {
  mid = mids_10[i]
  key = as.character(mid)
  n_total = format(sum(N_list_7[[key]]), big.mark = ",")
  paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100),
         "\n(N = ", n_total, ")",
         if (mid %in% thin_cohorts) "*" else "")
})

grid_df_7 = do.call(rbind, lapply(seq_along(mids_10), function(i) {
  mid = mids_10[i]
  key = as.character(mid)
  P   = P_list_7[[key]]
  if (is.null(P)) return(NULL)
  df  = as.data.frame(as.table(P))
  names(df) = c("origin", "current", "p")
  df$cohort  = cohort_labels_7[i]
  df
}))

grid_df_7$origin  = factor(grid_df_7$origin,  levels = rel_level_order_7)
grid_df_7$current = factor(grid_df_7$current, levels = rev(rel_level_order_7))
grid_df_7$cohort  = factor(grid_df_7$cohort,  levels = cohort_labels_7)

p_grid_7 = ggplot(grid_df_7, aes(current, origin, fill = p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p)), size = 2.4) +
  facet_wrap(~ cohort, nrow = 3) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                       name = "P[i→j]") +
  scale_x_discrete(labels = reltrad_labels_7) +
  scale_y_discrete(labels = reltrad_labels_7) +
  labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
       title = "Intergenerational transition matrices by birth cohort (7-state; 2022 & 2024 excluded)",
       caption = "* Cohort contains a Jewish or Black Protestant origin cell with N < 30.") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y      = element_text(size = 7),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 8),
    plot.caption     = element_text(size = 7, color = "grey50")
  )

ggsave("output/figures/7state-drop2224/P_grid_7state_10yr.png", p_grid_7,
       width = 10, height = 14, dpi = 200)

# ── FIGURE 2: λ₂ TREND WITH BOOTSTRAPPED CIs ─────────────────────────────────

lambda2 = function(P) sort(Mod(eigen(P)$values), decreasing = TRUE)[2]

n_boot = 250
set.seed(42)

boot_rows_7 = lapply(mids_10, function(mid) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == mid &
             !is.na(data$reltrad16) & !is.na(data$reltrad), ]

  P_obs  = suppressWarnings(
    p_matrix(sub, "reltrad16", "reltrad", levels = states_7)
  )
  l2_obs = lambda2(P_obs)

  l2_boot = replicate(n_boot, {
    idx = sample(nrow(sub), replace = TRUE)
    P_b = suppressWarnings(
      p_matrix(sub[idx, ], "reltrad16", "reltrad", levels = states_7)
    )
    lambda2(P_b)
  })

  data.frame(
    mid    = mid,
    cohort = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100)),
    l2     = l2_obs,
    lo     = quantile(l2_boot, 0.025),
    hi     = quantile(l2_boot, 0.975),
    n      = nrow(sub),
    row.names = NULL
  )
})

boot_df_7 = do.call(rbind, boot_rows_7)
boot_df_7$cohort = factor(boot_df_7$cohort, levels = boot_df_7$cohort)

cat("\nλ₂ estimates (7-state, 2022 & 2024 excluded) with 95% bootstrap CIs (n_boot = 250, seed = 42):\n")
print(boot_df_7[, c("cohort", "l2", "lo", "hi", "n")], row.names = FALSE, digits = 3)

p_l2_7 = ggplot(boot_df_7, aes(x = cohort, y = l2, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.15) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  color     = "#0072B2",
                  linewidth = 0.7,
                  size      = 0.6) +
  scale_y_continuous(limits = c(0.5, 0.9), breaks = seq(0.5, 0.9, 0.1)) +
  labs(
    x     = "Birth cohort (10-year window)",
    y     = expression(lambda[2]),
    title = expression(lambda[2] ~ "(7-state; 2022 & 2024 excluded) by birth cohort with 95% bootstrap CIs")
  ) +
  healy_theme

ggsave("output/figures/7state-drop2224/lambda2_trend_7state_10yr.png", p_l2_7,
       width = 8, height = 5, dpi = 200)

# ── FIGURE 3: π₀ AND π∞ DISTRIBUTION GRID ────────────────────────────────────

dist_df_7 = do.call(rbind, lapply(seq_along(mids_10), function(i) {
  mid    = mids_10[i]
  key    = as.character(mid)
  pi0    = pi0_list_7[[key]]
  pistar = pistar_list_7[[key]]
  if (is.null(pi0)) return(NULL)
  cohort_lbl = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100))
  rbind(
    data.frame(mid = mid, cohort = cohort_lbl, measure = "π₀  (origin)",
               religion = names(pi0),    value = as.numeric(pi0)),
    data.frame(mid = mid, cohort = cohort_lbl, measure = "π∞ (stationary)",
               religion = names(pistar), value = as.numeric(pistar))
  )
}))

dist_df_7$religion = factor(dist_df_7$religion, levels = rel_level_order_7)
dist_df_7$cohort   = factor(dist_df_7$cohort,
  levels = paste0(mids_10 - 5, "–", sprintf("%02d", (mids_10 + 4) %% 100)))
dist_df_7$measure  = factor(dist_df_7$measure,
  levels = c("π₀  (origin)", "π∞ (stationary)"))

p_dist_grid_7 = ggplot(dist_df_7, aes(y = religion, x = value, fill = religion)) +
  geom_col(width = 0.65) +
  facet_grid(measure ~ cohort) +
  scale_fill_manual(values = reltrad_colors_7, labels = reltrad_labels_7, name = NULL) +
  scale_x_continuous(limits = c(0, 0.65), breaks = c(0, 0.25, 0.5),
                     labels = c("0", ".25", ".5")) +
  scale_y_discrete(labels = reltrad_labels_7) +
  labs(x = "Share", y = NULL,
       title = "Origin (π₀) and stationary (π∞) distributions by birth cohort (7-state; 2022 & 2024 excluded)") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "none",
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(size = 8),
    axis.text.y        = element_text(size = 7),
    axis.text.x        = element_text(size = 7),
    plot.caption       = element_text(size = 8, color = "grey50")
  )

ggsave("output/figures/7state-drop2224/pi_dist_grid_7state_10yr.png", p_dist_grid_7,
       width = 14, height = 5, dpi = 200)

cat("\nSaved all 7-state (2022 & 2024 excluded) figures to output/figures/7state-drop2224/\n")

# ── FIGURE 4: DIFFERENCE GRID (P_full − P_drop2224) ──────────────────────────
# Re-estimate matrices on the full sample (all survey years) so the difference
# is computed within this script without depending on script 21's workspace.

P_list_full = list()

for (coh in mids_10) {
  sub = data_full[!is.na(data_full$cohort_10) & data_full$cohort_10 == coh &
                  !is.na(data_full$reltrad16) & !is.na(data_full$reltrad), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)
  P_list_full[[key]] = suppressWarnings(
    p_matrix(sub, "reltrad16", "reltrad", levels = states_7)
  )
}

diff_df_7 = do.call(rbind, lapply(seq_along(mids_10), function(i) {
  mid = mids_10[i]
  key = as.character(mid)
  P_f = P_list_full[[key]]
  P_d = P_list_7[[key]]
  if (is.null(P_f) || is.null(P_d)) return(NULL)
  diff_mat = P_f - P_d
  df = as.data.frame(as.table(diff_mat))
  names(df) = c("origin", "current", "diff")
  df$cohort = cohort_labels_7[i]
  df
}))

diff_df_7$origin  = factor(diff_df_7$origin,  levels = rel_level_order_7)
diff_df_7$current = factor(diff_df_7$current, levels = rev(rel_level_order_7))
diff_df_7$cohort  = factor(diff_df_7$cohort,  levels = cohort_labels_7)

diff_limit = max(abs(diff_df_7$diff), na.rm = TRUE)

p_diff_7 = ggplot(diff_df_7, aes(current, origin, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.3f", diff)), size = 2.2) +
  facet_wrap(~ cohort, nrow = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-diff_limit, diff_limit),
                       name = "P_full − P_drop") +
  scale_x_discrete(labels = reltrad_labels_7) +
  scale_y_discrete(labels = reltrad_labels_7) +
  labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
       title = "Difference in transition matrices: full sample minus 2022 & 2024 excluded (7-state)",
       caption = "Positive = cell probability is higher in the full sample (including 2022 & 2024).") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y      = element_text(size = 7),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 8),
    plot.caption     = element_text(size = 7, color = "grey50")
  )

ggsave("output/figures/7state-drop2224/P_diff_7state_10yr.png", p_diff_7,
       width = 10, height = 14, dpi = 200)

cat("Saved difference grid to output/figures/7state-drop2224/P_diff_7state_10yr.png\n")
