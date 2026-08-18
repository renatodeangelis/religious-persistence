# ── 18 · DIFFERENCE FIGURES: POLITICAL ────────────────────────────────────────
# Figure 1 (×2): Binary difference matrix grid — 6 cohort windows, nrow = 2
#   partyid_broad:  P_dem − P_rep   (excludes "other" from this figure only)
#   polviews_broad: P_liberal − P_conservative  (excludes "moderate")
# Figure 2 (×2): π₀ and π∞ dot-plot — all three groups per variable
#   partyid_broad:  dem / rep / other
#   polviews_broad: liberal / moderate / conservative
#
# Input:  data/derived/matrices_political.rds
# Output: output/figures/political/diff_grid_partyid_broad.png
#         output/figures/political/diff_grid_polviews_broad.png
#         output/figures/political/pi_comparison_partyid_broad.png
#         output/figures/political/pi_comparison_polviews_broad.png

library(ggplot2)
library(dplyr)
library(tidyr)
source("code/utils.R")

mat        = readRDS("data/derived/matrices_political.rds")
P_pol      = mat$P
pi0_pol    = mat$pi0
pistar_pol = mat$pistar
n_pol      = mat$n

mids = c(1930, 1940, 1950, 1960, 1970, 1980)

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

dir.create("output/figures/political", recursive = TRUE, showWarnings = FALSE)

# ── FIGURE 1 HELPER: binary difference matrix grid ────────────────────────────

make_diff_grid = function(vname, grp_a, grp_b, title_str, legend_label, outfile) {
  diff_rows = list()
  flag_log  = character(0)

  for (mid in mids) {
    key_a = paste(vname, grp_a, mid, sep = "_")
    key_b = paste(vname, grp_b, mid, sep = "_")
    if (is.null(P_pol[[key_a]]) || is.null(P_pol[[key_b]])) next

    D     = P_pol[[key_a]] - P_pol[[key_b]]
    n_min = min(n_pol[[key_a]], n_pol[[key_b]])
    idx   = which(abs(D) > 0.10 & n_min < 100, arr.ind = TRUE)
    if (nrow(idx) > 0) {
      for (k in seq_len(nrow(idx))) {
        flag_log = c(flag_log, sprintf(
          "Cohort %s  %s→%s  ΔP=%+.3f  N_%s=%d  N_%s=%d",
          cohort_lbl(mid),
          rownames(D)[idx[k, 1]], colnames(D)[idx[k, 2]],
          D[idx[k, 1], idx[k, 2]],
          grp_a, n_pol[[key_a]], grp_b, n_pol[[key_b]]
        ))
      }
    }

    df_d        = as.data.frame(as.table(D))
    names(df_d) = c("origin", "current", "diff")
    df_d$cohort = cohort_lbl(mid)
    df_d$mid    = mid
    diff_rows[[as.character(mid)]] = df_d
  }

  if (length(flag_log) > 0) {
    cat("\n── Reliability flags for", outfile, "(|ΔP| > 0.10 and N < 100) ──\n")
    cat(paste(flag_log, collapse = "\n"), "\n\n")
  } else {
    cat("No reliability flags for", basename(outfile), "\n")
  }

  diff_all         = do.call(rbind, diff_rows)
  diff_all$origin  = factor(diff_all$origin,  levels = rel_level_order)
  diff_all$current = factor(diff_all$current, levels = rev(rel_level_order))
  diff_all$cohort  = factor(diff_all$cohort,  levels = sapply(mids, cohort_lbl))

  lim_global = max(abs(diff_all$diff), na.rm = TRUE)

  p = ggplot(diff_all, aes(x = current, y = origin, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%+.3f", diff)), size = 2.8) +
    facet_wrap(~ cohort, nrow = 2) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-lim_global, lim_global),
      name = legend_label
    ) +
    labs(x = "Current religion", y = "Origin religion", title = title_str) +
    diff_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8)
    )

  ggsave(outfile, p, width = 12, height = 9, dpi = 200)
  cat("Saved", outfile, "\n")
}

# ── FIGURE 2 HELPER: π₀ and π∞ dot-plot, all groups ──────────────────────────
# Segments span min-to-max across groups for each origin × cohort × measure cell,
# showing the full range without implying a specific pairwise comparison.

make_pi_plot = function(vname, groups, group_colors, group_shapes,
                        group_labels, title_str, outfile) {
  pi_rows = list()

  for (mid in mids) {
    clbl = cohort_lbl(mid)
    for (grp in groups) {
      key = paste(vname, grp, mid, sep = "_")
      if (is.null(pi0_pol[[key]])) next

      for (measure in c("π₀  (origin)", "π∞ (stationary)")) {
        vec = if (measure == "π₀  (origin)") pi0_pol[[key]] else pistar_pol[[key]]
        pi_rows[[length(pi_rows) + 1]] = data.frame(
          cohort  = clbl,
          mid     = mid,
          measure = measure,
          group   = grp,
          origin  = names(vec),
          value   = as.numeric(vec),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  pi_df         = do.call(rbind, pi_rows)
  pi_df$origin  = factor(pi_df$origin,  levels = rel_level_order)
  pi_df$cohort  = factor(pi_df$cohort,  levels = sapply(mids, cohort_lbl))
  pi_df$measure = factor(pi_df$measure, levels = c("π₀  (origin)", "π∞ (stationary)"))
  pi_df$group   = factor(pi_df$group,   levels = groups)

  # Segments: min-to-max range across all groups per origin × cohort × measure
  seg_df = pi_df |>
    group_by(cohort, measure, origin) |>
    summarise(x_min = min(value), x_max = max(value), .groups = "drop")

  p = ggplot() +
    geom_segment(
      data = seg_df,
      aes(x = x_min, xend = x_max, y = origin, yend = origin),
      color = "grey70", linewidth = 0.5
    ) +
    geom_point(
      data = pi_df,
      aes(x = value, y = origin, color = group, shape = group),
      size = 2.2
    ) +
    facet_grid(measure ~ cohort) +
    scale_color_manual(values = group_colors, labels = group_labels, name = NULL) +
    scale_shape_manual(values = group_shapes, labels = group_labels, name = NULL) +
    scale_x_continuous(limits = c(0, NA), breaks = c(0, 0.25, 0.5),
                       labels = c("0", ".25", ".5")) +
    scale_y_discrete(labels = reltrad_labels_tc) +
    labs(x = "Share", y = NULL, title = title_str) +
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

  ggsave(outfile, p, width = 16, height = 5, dpi = 200)
  cat("Saved", outfile, "\n")
}

# ── FIGURE 1A: partyid_broad — P_dem − P_rep ──────────────────────────────────

make_diff_grid(
  vname        = "partyid_broad",
  grp_a        = "rep",
  grp_b        = "dem",
  title_str    = "Difference matrices: P Republican − P Democrat (broad party ID, 1925–1984)",
  legend_label = "ΔP (rep − dem)",
  outfile      = "output/figures/political/diff_grid_partyid_broad.png"
)

# ── FIGURE 1B: polviews_broad — P_liberal − P_conservative ────────────────────

make_diff_grid(
  vname        = "polviews_broad",
  grp_a        = "conservative",
  grp_b        = "liberal",
  title_str    = "Difference matrices: P Conservative − P Liberal (broad political views, 1925–1984)",
  legend_label = "ΔP (conservative − liberal)",
  outfile      = "output/figures/political/diff_grid_polviews_broad.png"
)

# ── FIGURE 2A: partyid_broad — dem / rep / other ──────────────────────────────

make_pi_plot(
  vname        = "partyid_broad",
  groups       = c("dem", "rep", "other"),
  group_colors = c(dem = "#0072B2", rep = "#D55E00", other = "#999999"),
  group_shapes = c(dem = 16, rep = 17, other = 15),
  group_labels = c(dem = "Democrat", rep = "Republican", other = "Other/Ind."),
  title_str    = "π₀ and π∞ by party ID (broad) and birth cohort (1925–1984)",
  outfile      = "output/figures/political/pi_comparison_partyid_broad.png"
)

# ── FIGURE 2B: polviews_broad — liberal / moderate / conservative ──────────────

make_pi_plot(
  vname        = "polviews_broad",
  groups       = c("liberal", "moderate", "conservative"),
  group_colors = c(liberal = "#0072B2", moderate = "#999999", conservative = "#D55E00"),
  group_shapes = c(liberal = 16, moderate = 15, conservative = 17),
  group_labels = c(liberal = "Liberal", moderate = "Moderate", conservative = "Conservative"),
  title_str    = "π₀ and π∞ by political views (broad) and birth cohort (1925–1984)",
  outfile      = "output/figures/political/pi_comparison_polviews_broad.png"
)

# ── FIGURE 1C: partyid_narrow — P_dem − P_rep ─────────────────────────────────

make_diff_grid(
  vname        = "partyid_narrow",
  grp_a        = "rep",
  grp_b        = "dem",
  title_str    = "Difference matrices: P Republican − P Democrat (narrow party ID, 1925–1984)",
  legend_label = "ΔP (rep − dem)",
  outfile      = "output/figures/political/diff_grid_partyid_narrow.png"
)

# ── FIGURE 1D: polviews_narrow — P_liberal − P_conservative ───────────────────

make_diff_grid(
  vname        = "polviews_narrow",
  grp_a        = "conservative",
  grp_b        = "liberal",
  title_str    = "Difference matrices: P Conservative − P Liberal (narrow political views, 1925–1984)",
  legend_label = "ΔP (conservative − liberal)",
  outfile      = "output/figures/political/diff_grid_polviews_narrow.png"
)

# ── FIGURE 2C: partyid_narrow — dem / rep / other ─────────────────────────────

make_pi_plot(
  vname        = "partyid_narrow",
  groups       = c("dem", "rep", "other"),
  group_colors = c(dem = "#0072B2", rep = "#D55E00", other = "#999999"),
  group_shapes = c(dem = 16, rep = 17, other = 15),
  group_labels = c(dem = "Democrat", rep = "Republican", other = "Other/Ind."),
  title_str    = "π₀ and π∞ by party ID (narrow) and birth cohort (1925–1984)",
  outfile      = "output/figures/political/pi_comparison_partyid_narrow.png"
)

# ── FIGURE 2D: polviews_narrow — liberal / moderate / conservative ─────────────

make_pi_plot(
  vname        = "polviews_narrow",
  groups       = c("liberal", "moderate", "conservative"),
  group_colors = c(liberal = "#0072B2", moderate = "#999999", conservative = "#D55E00"),
  group_shapes = c(liberal = 16, moderate = 15, conservative = 17),
  group_labels = c(liberal = "Liberal", moderate = "Moderate", conservative = "Conservative"),
  title_str    = "π₀ and π∞ by political views (narrow) and birth cohort (1925–1984)",
  outfile      = "output/figures/political/pi_comparison_polviews_narrow.png"
)
