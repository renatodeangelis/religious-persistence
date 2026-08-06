# ── 04 · HOMOGENEITY TESTS ─────────────────────────────────────────────────────
# Anderson-Goodman (1957) chi-square tests of whether the RELIG16 -> RELIG
# transition matrix is common across birth cohorts: omnibus, row-wise,
# pairwise-adjacent, and a rolling-window sensitivity. Emits three figures and
# one LaTeX table. Operates on RAW COUNTS (cell frequencies), consuming
# nat10$N for all tests (six 10-year cohort windows: 1925–1984).
#
# Input:  data/derived/matrices.rds
# Output: output/figures/homogeneity/*.png, output/tables/homogeneity_tests.tex

source("code/utils.R")

matrices  = readRDS("data/derived/matrices.rds")
N_list_10 = matrices$nat10$N

# ── HOMOGENEITY TEST ACROSS COHORT DECADES (Anderson-Goodman 1957) ───────────
# H0: the RELIG16 -> RELIG transition matrix is common across the 10-year birth
# cohorts. Because each cohort matrix is an independent GSS cross-section, this
# is a test of HOMOGENEITY across cohorts, not stationarity over time — a
# rejection means the one-step matrix has genuinely shifted across cohorts,
# beyond sampling noise. (The chi2_joint console message inherited from the
# occupation-mobility port says "not stationary"; read it as "not homogeneous".)

# Omnibus test: are all cohort-decade matrices drawn from a common P?
cat("\n════ Omnibus homogeneity test across 10-year cohorts ════\n")
homog_joint_10 = chi2_joint(N_list_10, alpha = 0.05)

# Row-wise decomposition: which origin states drive non-homogeneity? Needs no
# ordinality — respects the nominal state space.
homog_row_10 = chi2_row(N_list_10, alpha = 0.05)
cat("\n── Row-wise homogeneity by origin state (10-year cohorts) ──\n")
print(homog_row_10, row.names = FALSE)

# Pairwise-adjacent tests: localize WHICH cohort-to-cohort step carries the
# change. Each test compares two consecutive cohort decades.
cohort_keys_10 = names(N_list_10)
pairwise_10 = do.call(rbind, lapply(seq_len(length(cohort_keys_10) - 1), function(i) {
  k1  = cohort_keys_10[i]
  k2  = cohort_keys_10[i + 1]
  res = chi2_joint(N_list_10[c(k1, k2)], alpha = 0.05)
  data.frame(
    from_cohort = as.integer(k1),
    to_cohort   = as.integer(k2),
    chi2        = res$chi2,
    df          = res$df,
    p_value     = res$p_value,
    significant = res$significant,
    row.names   = NULL
  )
}))
cat("\n── Pairwise-adjacent homogeneity tests (10-year cohorts) ──\n")
print(pairwise_10, row.names = FALSE)

# ── HOMOGENEITY FIGURES ──────────────────────────────────────────────────────

dir.create("output/figures/homogeneity", recursive = TRUE, showWarnings = FALSE)

# ── ROLLING-WINDOW SENSITIVITY (10-year cohort bins) ─────────────────────────
# Slide a window of k consecutive 10-year cohort matrices across N_list_10 and
# run the joint homogeneity test inside each window. With six 10-year bins,
# k = 2, 3, 4 are the useful widths (k = 5 or 6 collapses toward the omnibus).
# A dip below p = 0.05 marks a neighborhood where the transition matrix is
# changing fast; comparing widths shows whether that localization is robust.

joint_stat = function(mats) {
  rw = chi2_row(mats)
  jc = sum(rw$chi2); jd = sum(rw$df)
  c(chi2 = jc, df = jd, p = pchisq(jc, jd, lower.tail = FALSE))
}

roll_years = as.numeric(names(N_list_10))   # ordered midpoints: 1930, 1940, ..., 1980
bin_widths = c(2, 3, 4)

roll_results = do.call(rbind, lapply(bin_widths, function(k) {
  half      = floor((k - 1) / 2)
  valid_idx = (half + 1):(length(roll_years) - (k - half - 1))
  do.call(rbind, lapply(valid_idx, function(i) {
    window_idx = (i - half):(i - half + k - 1)
    mats = N_list_10[window_idx]
    if (length(mats) < 2) return(NULL)
    jt = joint_stat(mats)
    data.frame(
      k          = k,
      center     = roll_years[i],
      year_min   = roll_years[min(window_idx)],
      year_max   = roll_years[max(window_idx)],
      n_matrices = length(mats),
      chi2       = jt[["chi2"]],
      p_value    = jt[["p"]],
      row.names  = NULL
    )
  }))
}))

roll_results$k_label = factor(paste0("k = ", roll_results$k, " bins"),
                              levels = paste0("k = ", bin_widths, " bins"))

# x-axis: center midpoint labeled as its bin range ("1935–44" for mid = 1940)
roll_x_breaks = c(1930, 1940, 1950, 1960, 1970, 1980)
roll_x_labels = paste0(roll_x_breaks - 5, "–",
                       sprintf("%02d", (roll_x_breaks + 4) %% 100))

p_roll = ggplot(roll_results, aes(x = center, y = p_value,
                                  color = k_label, shape = k_label, group = k_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.65, alpha = 0.85) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("indianred4", "goldenrod2", "royalblue4"),
    name   = "Window width\n(10-yr bins)") +
  scale_shape_manual(values = c(15, 16, 17),
                     name   = "Window width\n(10-yr bins)") +
  scale_x_continuous(breaks = roll_x_breaks, labels = roll_x_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    x = "Center cohort bin",
    y = "Joint homogeneity test p-value") +
  healy_theme +
  theme(legend.key.width = unit(1.2, "cm"))

# ── ROLLING-WINDOW SENSITIVITY (5-year cohort bins) ───────────────────────────
# Builds N_list_5 inline from gss_clean.rds (not saved in matrices.rds) and
# reruns the same rolling-window test at k = 2, 3, 4. Provides finer resolution
# at the cost of smaller N per bin; compare with the 10-year figure above.

clean_5yr  = readRDS("data/derived/gss_clean.rds")
data_5yr   = clean_5yr$data
states_5yr = clean_5yr$states_alt

data_5yr$cohort_5 = (floor((data_5yr$cohort - 1925) / 5) * 5 + 1925) + 2
mids_5 = sort(unique(data_5yr$cohort_5[!is.na(data_5yr$cohort_5)]))

N_list_5 = list()
for (mid in mids_5) {
  sub = data_5yr[!is.na(data_5yr$cohort_5)    & data_5yr$cohort_5      == mid &
                   !is.na(data_5yr$reltrad16_alt) & !is.na(data_5yr$reltrad_alt), ]
  if (nrow(sub) < 30) next
  N_list_5[[as.character(mid)]] = count_matrix(sub, "reltrad16_alt", "reltrad_alt",
                                               levels = states_5yr)
}

roll_years_5 = as.numeric(names(N_list_5))

roll_results_5 = do.call(rbind, lapply(bin_widths, function(k) {
  half      = floor((k - 1) / 2)
  valid_idx = (half + 1):(length(roll_years_5) - (k - half - 1))
  do.call(rbind, lapply(valid_idx, function(i) {
    window_idx = (i - half):(i - half + k - 1)
    mats = N_list_5[window_idx]
    if (length(mats) < 2) return(NULL)
    jt = joint_stat(mats)
    data.frame(k       = k,
               center  = roll_years_5[i],
               chi2    = jt[["chi2"]],
               p_value = jt[["p"]],
               row.names = NULL)
  }))
}))

roll_results_5$k_label = factor(paste0("k = ", roll_results_5$k, " bins"),
                                levels = paste0("k = ", bin_widths, " bins"))

roll_x_breaks_5 = as.numeric(names(N_list_5))
roll_x_labels_5 = paste0(roll_x_breaks_5 - 2, "–",
                          sprintf("%02d", (roll_x_breaks_5 + 2) %% 100))

p_roll_5 = ggplot(roll_results_5, aes(x = center, y = p_value,
                                       color = k_label, shape = k_label, group = k_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.65, alpha = 0.85) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c("indianred4", "goldenrod2", "royalblue4"),
                     name   = "Window width\n(5-yr bins)") +
  scale_shape_manual(values = c(15, 16, 17),
                     name   = "Window width\n(5-yr bins)") +
  scale_x_continuous(breaks = roll_x_breaks_5, labels = roll_x_labels_5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Center cohort bin",
       y = "Joint homogeneity test p-value") +
  healy_theme +
  theme(legend.key.width = unit(1.2, "cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 9))

p_roll_combined = patchwork::wrap_plots(p_roll, p_roll_5, ncol = 1) +
  patchwork::plot_annotation(
    tag_levels = "A",
    caption    = paste("Joint homogeneity tests (Anderson-Goodman 1957) on k consecutive cohort",
                       "matrices. Panel A: 10-year bins (1925–1984); Panel B: 5-year bins (1925–1984).",
                       "k = window width in number of matrices. Dashed line: p = 0.05."),
    theme = theme(plot.caption = element_text(size = 9, color = "grey40", hjust = 0))
  )

ggsave("output/figures/homogeneity/rolling_window_sensitivity.png", p_roll_combined,
       width = 14, height = 12, dpi = 300)

# ── HOMOGENEITY LATEX TABLE (grid: origins x cohort transitions) ─────────────
# Cells show chi2 with significance stars; the All column is the omnibus across
# every cohort, each remaining column a successive-decade pair (Anderson-Goodman 1957).

# Significance stars as a LaTeX superscript (cells show chi2; stars encode p)
star_only = function(p) ifelse(p < 0.001, "$^{***}$",
                        ifelse(p < 0.01,  "$^{**}$",
                        ifelse(p < 0.05,  "$^{*}$", "")))

# per-origin chi2/p (+ joint) for one count-list, origins in rel_level_order
grid_col_stats = function(N_list) {
  rw = chi2_row(N_list)
  rw = rw[match(rel_level_order, rw$state), ]
  jc = sum(rw$chi2); jd = sum(rw$df)
  list(chi2       = setNames(rw$chi2,    rel_level_order),
       p          = setNames(rw$p_value, rel_level_order),
       joint_chi2 = jc,
       joint_p    = pchisq(jc, jd, lower.tail = FALSE))
}

# Columns: all-cohort omnibus, then each successive 10-year-cohort pair
cohort_keys_10 = names(N_list_10)
pair_from = as.integer(cohort_keys_10[-length(cohort_keys_10)])
pair_to   = as.integer(cohort_keys_10[-1])
pair_lbls = sprintf("%d--%02d", pair_from, pair_to %% 100)

grid_cols = c(list(N_list_10),
              lapply(seq_along(pair_from), function(i)
                N_list_10[cohort_keys_10[c(i, i + 1)]]))
col_lbls  = c("All", pair_lbls)
col_stats = lapply(grid_cols, grid_col_stats)

# Body rows (one per origin) + bold Joint row
origin_rows = vapply(rel_level_order, function(st) {
  cells = vapply(col_stats, function(cs)
    sprintf("%.1f%s", cs$chi2[[st]], star_only(cs$p[[st]])), character(1))
  paste0(tools::toTitleCase(st), " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1))

joint_cells = vapply(col_stats, function(cs)
  sprintf("\\textbf{%.1f%s}", cs$joint_chi2, star_only(cs$joint_p)), character(1))
joint_row = paste0("\\textbf{Joint} & ", paste(joint_cells, collapse = " & "), " \\\\")

# Dynamic column spec / header (width adapts to the number of cohort pairs)
col_spec   = paste0("l", strrep("r", length(col_lbls)))
header_row = paste0("Origin & ", paste(col_lbls, collapse = " & "), " \\\\")
df_all     = (length(rel_level_order) - 1) * (length(N_list_10) - 1)

homog_tex = c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Chi-Square Tests of Homogeneity Across Birth Cohorts}",
  "\\label{tab:homogeneity}",
  paste0("\\begin{tabular}{", col_spec, "}"),
  "\\toprule",
  header_row,
  "\\midrule",
  origin_rows,
  "\\midrule",
  joint_row,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{0.5em}\\footnotesize",
  "\\textit{Notes:} Cells report $\\chi^2$ statistics testing whether the RELIG16 $\\to$ RELIG",
  "transition probabilities out of each origin state are constant across birth cohorts",
  "(Anderson and Goodman 1957). The \\textit{All} column pools every 10-year cohort; each",
  "remaining column compares two successive cohort decades, localizing when the matrix",
  paste0("shifts. Each origin-state test has df $= (s-1)(T-1)$ --- 4 for the pairwise columns and ",
         df_all, " for \\textit{All} --- reduced where a pooled destination cell is empty; the"),
  "Joint row sums $\\chi^2$ and df across origins. Because each cohort matrix is an independent",
  "GSS cross-section, these are tests of homogeneity across cohorts. The pairwise columns",
  "involve multiple comparisons; read a lone $^{*}$ cautiously.",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)


# ── STATIONARITY TABLE: ROW-SPECIFIC AND JOINT, TWO PANELS ───────────────────
# Panel A: all six decadal cohorts (N_list_10, midpoints 1930–1980).
# Panel B: central three cohorts only (midpoints 1940, 1950, 1960 = bins 1935–64).
# Columns: Origin State | chi2 | df | p-value (stars on p-value column only).
# Reuses star_only() and chi2_row() already defined above; avoids chi2_joint()
# to suppress its console output.

make_panel_rows = function(N_sub) {
  rw = chi2_row(N_sub)
  rw = rw[match(rel_level_order, rw$state), ]
  jc = sum(rw$chi2); jd = sum(rw$df)
  jp = pchisq(jc, jd, lower.tail = FALSE)

  origin_rows = vapply(seq_len(nrow(rw)), function(i) {
    sprintf("%s & %.2f & %d & %.3f%s \\\\",
            tools::toTitleCase(rw$state[i]),
            rw$chi2[i], rw$df[i],
            rw$p_value[i], star_only(rw$p_value[i]))
  }, character(1))

  joint_row = sprintf(
    "\\textbf{Joint} & \\textbf{%.2f} & \\textbf{%d} & \\textbf{%.3f}%s \\\\",
    jc, jd, jp, star_only(jp))

  c(origin_rows, joint_row)
}

panel_A = make_panel_rows(N_list_10)
panel_B = make_panel_rows(N_list_10[c("1940", "1950", "1960")])

stat_tex = c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Chi-Square Tests of Homogeneity: Row-Specific and Joint Statistics}",
  "\\label{tab:stationarity}",
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  "Origin State & $\\chi^2$ & df & $p$-value \\\\",
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Panel A: All cohorts (1925--1984)}} \\\\",
  "\\midrule",
  panel_A,
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Panel B: Central cohorts (1935--1964)}} \\\\",
  "\\midrule",
  panel_B,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{0.5em}\\footnotesize",
  "\\textit{Notes:} Each row reports the Anderson-Goodman (1957) chi-square statistic",
  "testing whether the RELIG16 $\\to$ RELIG transition probabilities out of origin",
  "state $i$ are constant across the indicated birth-cohort windows. Panel A includes",
  "all six 10-year cohort bins (1925--1984, $T=6$); Panel B restricts to the three",
  "central bins (1935--1964, $T=3$). df $= (s-1)(T-1)$, reduced where a pooled",
  "destination cell is empty. The Joint row sums $\\chi^2$ and df across all origin",
  "states. $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)

writeLines(stat_tex, "output/figures/homogeneity/homogeneity_stationarity.tex")
cat("\nWrote output/figures/homogeneity/homogeneity_stationarity.tex\n")

# ── PAIRWISE COMPARISON TABLE ─────────────────────────────────────────────────
# Same grid format as the old homogeneity_tests.tex (chi2 + stars, no separate
# df/p columns) but with the five adjacent-decade pairs first, all-cohorts last.
# Reuses col_stats / pair_lbls already computed above.

pw_col_stats = col_stats[c(2, 3, 4, 5, 6, 1)]   # pairwise first, all-cohorts last
pw_col_lbls  = c(pair_lbls, "All")

pw_col_spec  = paste0("l", strrep("r", length(pw_col_lbls)))
pw_hdr       = paste0("Origin & ", paste(pw_col_lbls, collapse = " & "), " \\\\")

pw_origin_rows = vapply(rel_level_order, function(st) {
  cells = vapply(pw_col_stats, function(cs)
    sprintf("%.1f%s", cs$chi2[[st]], star_only(cs$p[[st]])), character(1))
  paste0(tools::toTitleCase(st), " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1))

pw_joint_cells = vapply(pw_col_stats, function(cs)
  sprintf("\\textbf{%.1f%s}", cs$joint_chi2, star_only(cs$joint_p)), character(1))
pw_joint_row = paste0("\\textbf{Joint} & ", paste(pw_joint_cells, collapse = " & "), " \\\\")

pw_df_pair = (length(rel_level_order) - 1) * (2 - 1)
pw_df_all  = (length(rel_level_order) - 1) * (length(N_list_10) - 1)

pw_tex = c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Chi-Square Tests of Homogeneity: Adjacent-Cohort Comparisons}",
  "\\label{tab:pairwise}",
  paste0("\\begin{tabular}{", pw_col_spec, "}"),
  "\\toprule",
  pw_hdr,
  "\\midrule",
  pw_origin_rows,
  "\\midrule",
  pw_joint_row,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{0.5em}\\footnotesize",
  "\\textit{Notes:} Cells report $\\chi^2$ statistics (Anderson-Goodman 1957).",
  paste0("Each pairwise column tests two consecutive 10-year cohort matrices;",
         " df $= (s-1)(T-1) = ", pw_df_pair, "$ per origin (Joint df $= ",
         pw_df_pair * length(rel_level_order), "$)."),
  paste0("The \\textit{All} column pools all six cohorts; df $= ", pw_df_all,
         "$ per origin (Joint df $= ", pw_df_all * length(rel_level_order), "$)."),
  "Both are reduced where a pooled destination cell is empty.",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)

writeLines(pw_tex, "output/figures/homogeneity/homogeneity_pairwise.tex")
cat("\nWrote output/figures/homogeneity/homogeneity_pairwise.tex\n")
