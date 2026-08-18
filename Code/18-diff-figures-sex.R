# ── 16 · DIFFERENCE FIGURES: SEX ──────────────────────────────────────────────
# Figure 1: Difference matrix grid (P_male − P_female) across 6 cohort windows,
#           laid out like P_grid_national (facet_wrap, nrow = 2).
# Figure 2: π₀ and π∞ dot-plot comparison, male vs female, facet_grid(measure ~ cohort).
#
# Input:  data/derived/matrices_sex.rds
# Output: output/figures/sex/diff_grid_sex.png
#         output/figures/sex/pi_comparison_sex.png

library(ggplot2)
library(dplyr)
library(tidyr)
source("code/utils.R")

mat        = readRDS("data/derived/matrices_sex.rds")
P_sex      = mat$P
pi0_sex    = mat$pi0
pistar_sex = mat$pistar
n_sex      = mat$n

mids = c(1930, 1940, 1950, 1960, 1970, 1980)    # cohort_10 midpoints (edges 1925–1984)

cohort_lbl = function(mid)
  paste0(mid - 5, "–9", sprintf("%d", ((mid - 5) %/% 100) * 100 + (mid + 4) %% 100))

# Use the cleaner format used throughout the project
cohort_lbl = function(mid) {
  edge = mid - 5
  paste0(edge, "–", sprintf("%02d", (edge + 9) %% 100))
}

rel_level_order   = c("catholic", "evangelical", "mainline", "other", "none")
reltrad_labels_tc = c(
  catholic = "Catholic", evangelical = "Evangelical", mainline = "Mainline",
  other = "Other", none = "None"
)

diff_theme = theme_bw(base_size = 10) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 9)
  )

# ── FIGURE 1: DIFFERENCE MATRIX GRID ─────────────────────────────────────────

diff_rows = list()
flag_log  = character(0)

for (mid in mids) {
  key_m = paste0("male_",   mid)
  key_f = paste0("female_", mid)
  if (is.null(P_sex[[key_m]]) || is.null(P_sex[[key_f]])) next

  D = P_sex[[key_m]] - P_sex[[key_f]]

  # Flag cells where |ΔP| > 0.10 and the smaller group N < 100
  n_min = min(n_sex[[key_m]], n_sex[[key_f]])
  idx   = which(abs(D) > 0.10 & n_min < 100, arr.ind = TRUE)
  if (nrow(idx) > 0) {
    for (k in seq_len(nrow(idx))) {
      flag_log = c(flag_log, sprintf(
        "Cohort %s  %s→%s  ΔP=%+.3f  N_male=%d  N_female=%d",
        cohort_lbl(mid),
        rownames(D)[idx[k, 1]], colnames(D)[idx[k, 2]],
        D[idx[k, 1], idx[k, 2]],
        n_sex[[key_m]], n_sex[[key_f]]
      ))
    }
  }

  df_d           = as.data.frame(as.table(D))
  names(df_d)    = c("origin", "current", "diff")
  df_d$cohort    = cohort_lbl(mid)
  df_d$mid       = mid
  diff_rows[[as.character(mid)]] = df_d
}

if (length(flag_log) > 0) {
  cat("\n── Reliability flags (|ΔP| > 0.10 and N < 100) ──\n")
  cat(paste(flag_log, collapse = "\n"), "\n\n")
} else {
  cat("No reliability flags.\n")
}

diff_all         = do.call(rbind, diff_rows)
diff_all$origin  = factor(diff_all$origin,  levels = rel_level_order)
diff_all$current = factor(diff_all$current, levels = rev(rel_level_order))
diff_all$cohort  = factor(diff_all$cohort,  levels = sapply(mids, cohort_lbl))

lim_global = max(abs(diff_all$diff), na.rm = TRUE)

p_diff_grid = ggplot(diff_all, aes(x = current, y = origin, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.3f", diff)), size = 2.8) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-lim_global, lim_global),
    name = "ΔP (male − female)"
  ) +
  labs(
    x = "Current religion", y = "Origin religion",
    title = "Difference matrices: P male − P female (1925–1984)"
  ) +
  diff_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

dir.create("output/figures/sex", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/sex/diff_grid_sex.png", p_diff_grid,
       width = 12, height = 9, dpi = 200)
cat("Saved output/figures/sex/diff_grid_sex.png\n")

# ── FIGURE 2: π₀ AND π∞ DOT-PLOT COMPARISON ──────────────────────────────────

pi_rows = list()

for (mid in mids) {
  key_m = paste0("male_",   mid)
  key_f = paste0("female_", mid)
  if (is.null(pi0_sex[[key_m]]) || is.null(pi0_sex[[key_f]])) next

  clbl = cohort_lbl(mid)

  for (measure in c("π₀  (origin)", "π∞ (stationary)")) {
    vec_m = if (measure == "π₀  (origin)") pi0_sex[[key_m]] else pistar_sex[[key_m]]
    vec_f = if (measure == "π₀  (origin)") pi0_sex[[key_f]] else pistar_sex[[key_f]]

    pi_rows[[length(pi_rows) + 1]] = data.frame(
      cohort  = clbl,
      mid     = mid,
      measure = measure,
      group   = rep(c("Male", "Female"), each = length(vec_m)),
      origin  = rep(names(vec_m), 2),
      value   = c(as.numeric(vec_m), as.numeric(vec_f)),
      stringsAsFactors = FALSE
    )
  }
}

pi_df         = do.call(rbind, pi_rows)
pi_df$origin  = factor(pi_df$origin,  levels = rel_level_order)
pi_df$cohort  = factor(pi_df$cohort,  levels = sapply(mids, cohort_lbl))
pi_df$measure = factor(pi_df$measure, levels = c("π₀  (origin)", "π∞ (stationary)"))
pi_df$group   = factor(pi_df$group,   levels = c("Male", "Female"))

# Wide form for segments connecting Male and Female points
pi_wide = pi_df |>
  pivot_wider(names_from = group, values_from = value)

p_pi_comp = ggplot() +
  geom_segment(
    data = pi_wide,
    aes(x = Male, xend = Female, y = origin, yend = origin),
    color = "grey70", linewidth = 0.5
  ) +
  geom_point(
    data = pi_df,
    aes(x = value, y = origin, color = group, shape = group),
    size = 2.2
  ) +
  facet_grid(measure ~ cohort) +
  scale_color_manual(values = c(Male = "#0072B2", Female = "#D55E00"), name = NULL) +
  scale_shape_manual(values = c(Male = 16, Female = 17), name = NULL) +
  scale_x_continuous(limits = c(0, NA), breaks = c(0, 0.25, 0.5),
                     labels = c("0", ".25", ".5")) +
  scale_y_discrete(labels = reltrad_labels_tc) +
  labs(
    x = "Share", y = NULL,
    title = "π₀ and π∞ by sex and birth cohort (1925–1984)"
  ) +
  healy_theme +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text.x       = element_text(size = 8),
    strip.text.y       = element_text(size = 9),
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 7)
  )

ggsave("output/figures/sex/pi_comparison_sex.png", p_pi_comp,
       width = 16, height = 5, dpi = 200)
cat("Saved output/figures/sex/pi_comparison_sex.png\n")
