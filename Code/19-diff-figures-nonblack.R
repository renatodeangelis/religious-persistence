# ── 17 · DIFFERENCE FIGURES: NON-BLACK ────────────────────────────────────────
# Figure 1: Difference matrix grid (P_nonblack − P_full) across 6 cohort windows,
#           laid out like P_grid_national (facet_wrap, nrow = 2).
# Figure 2: π₀ and π∞ dot-plot comparison, non-Black vs full sample,
#           facet_grid(measure ~ cohort).
#
# Input:  data/derived/matrices_nonblack.rds
# Output: output/figures/nonblack/diff_grid_nonblack.png
#         output/figures/nonblack/pi_comparison_nonblack.png

library(ggplot2)
library(dplyr)
library(tidyr)
source("code/utils.R")

mat    = readRDS("data/derived/matrices_nonblack.rds")
nb     = mat$nonblack   # keys are edge labels: "1925", "1935", ..., "1975"
full   = mat$full

# Edge labels (keys) and corresponding midpoints for cohort_10 binning
edges  = c("1925", "1935", "1945", "1955", "1965", "1975")

cohort_lbl = function(edge)
  paste0(edge, "–", sprintf("%02d", (as.integer(edge) + 9) %% 100))

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

for (edge in edges) {
  if (is.null(nb$P[[edge]]) || is.null(full$P[[edge]])) next

  D = nb$P[[edge]] - full$P[[edge]]

  # Flag cells where |ΔP| > 0.10 and non-Black N < 100 (the smaller, limiting group)
  n_nb = nb$n[[edge]]
  idx  = which(abs(D) > 0.10 & n_nb < 100, arr.ind = TRUE)
  if (nrow(idx) > 0) {
    for (k in seq_len(nrow(idx))) {
      flag_log = c(flag_log, sprintf(
        "Cohort %s  %s→%s  ΔP=%+.3f  N_nonblack=%d  N_full=%d",
        cohort_lbl(edge),
        rownames(D)[idx[k, 1]], colnames(D)[idx[k, 2]],
        D[idx[k, 1], idx[k, 2]],
        nb$n[[edge]], full$n[[edge]]
      ))
    }
  }

  df_d        = as.data.frame(as.table(D))
  names(df_d) = c("origin", "current", "diff")
  df_d$cohort = cohort_lbl(edge)
  df_d$edge   = edge
  diff_rows[[edge]] = df_d
}

if (length(flag_log) > 0) {
  cat("\n── Reliability flags (|ΔP| > 0.10 and N_nonblack < 100) ──\n")
  cat(paste(flag_log, collapse = "\n"), "\n\n")
} else {
  cat("No reliability flags.\n")
}

diff_all         = do.call(rbind, diff_rows)
diff_all$origin  = factor(diff_all$origin,  levels = rel_level_order)
diff_all$current = factor(diff_all$current, levels = rev(rel_level_order))
diff_all$cohort  = factor(diff_all$cohort,  levels = sapply(edges, cohort_lbl))

lim_global = max(abs(diff_all$diff), na.rm = TRUE)

p_diff_grid = ggplot(diff_all, aes(x = current, y = origin, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.3f", diff)), size = 2.8) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-lim_global, lim_global),
    name = "ΔP (non-Black − full)"
  ) +
  labs(
    x = "Current religion", y = "Origin religion",
    title = "Difference matrices: P non-Black − P full sample (1925–1984)"
  ) +
  diff_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

dir.create("output/figures/nonblack", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/nonblack/diff_grid_nonblack.png", p_diff_grid,
       width = 12, height = 9, dpi = 200)
cat("Saved output/figures/nonblack/diff_grid_nonblack.png\n")

# ── FIGURE 2: π₀ AND π∞ DOT-PLOT COMPARISON ──────────────────────────────────

pi_rows = list()

for (edge in edges) {
  if (is.null(nb$pi0[[edge]]) || is.null(full$pi0[[edge]])) next

  clbl = cohort_lbl(edge)

  for (measure in c("π₀  (origin)", "π∞ (stationary)")) {
    vec_nb  = if (measure == "π₀  (origin)") nb$pi0[[edge]]   else nb$pistar[[edge]]
    vec_all = if (measure == "π₀  (origin)") full$pi0[[edge]] else full$pistar[[edge]]

    pi_rows[[length(pi_rows) + 1]] = data.frame(
      cohort  = clbl,
      edge    = edge,
      measure = measure,
      group   = rep(c("Non-Black", "Full"), each = length(vec_nb)),
      origin  = rep(names(vec_nb), 2),
      value   = c(as.numeric(vec_nb), as.numeric(vec_all)),
      stringsAsFactors = FALSE
    )
  }
}

pi_df         = do.call(rbind, pi_rows)
pi_df$origin  = factor(pi_df$origin,  levels = rel_level_order)
pi_df$cohort  = factor(pi_df$cohort,  levels = sapply(edges, cohort_lbl))
pi_df$measure = factor(pi_df$measure, levels = c("π₀  (origin)", "π∞ (stationary)"))
pi_df$group   = factor(pi_df$group,   levels = c("Non-Black", "Full"))

pi_wide = pi_df |>
  pivot_wider(names_from = group, values_from = value)

p_pi_comp = ggplot() +
  geom_segment(
    data = pi_wide,
    aes(x = `Non-Black`, xend = Full, y = origin, yend = origin),
    color = "grey70", linewidth = 0.5
  ) +
  geom_point(
    data = pi_df,
    aes(x = value, y = origin, color = group, shape = group),
    size = 2.2
  ) +
  facet_grid(measure ~ cohort) +
  scale_color_manual(values = c("Non-Black" = "#0072B2", "Full" = "#999999"), name = NULL) +
  scale_shape_manual(values  = c("Non-Black" = 16,        "Full" = 1),         name = NULL) +
  scale_x_continuous(limits = c(0, NA), breaks = c(0, 0.25, 0.5),
                     labels = c("0", ".25", ".5")) +
  scale_y_discrete(labels = reltrad_labels_tc) +
  labs(
    x = "Share", y = NULL,
    title = "π₀ and π∞ by sample (non-Black vs full) and birth cohort (1925–1984)"
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

ggsave("output/figures/nonblack/pi_comparison_nonblack.png", p_pi_comp,
       width = 16, height = 5, dpi = 200)
cat("Saved output/figures/nonblack/pi_comparison_nonblack.png\n")
