# ── 04 · HOMOGENEITY TESTS ─────────────────────────────────────────────────────
# Anderson-Goodman (1957) chi-square tests of whether the RELIG16 -> RELIG
# transition matrix is common across birth cohorts: omnibus, row-wise, pairwise-
# adjacent, and a rolling-window sensitivity. Emits two figures and one LaTeX
# table. Operates on RAW COUNTS (cell frequencies), so it consumes the N slots
# persisted in 02 — nat10$N for the decadal tests, nat5$N for the rolling window.
#
# Input:  data/derived/matrices.rds
# Output: output/figures/homogeneity/*.png, output/tables/homogeneity_tests.tex

source("code/utils.R")

matrices  = readRDS("data/derived/matrices.rds")
N_list_10 = matrices$nat10$N
N_list_5  = matrices$nat5$N

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

# Row-wise: chi-square magnitude per origin. The p-value axis is uninformative
# here — every origin rejects at 10-year resolution, so all bars would pin to
# zero — so plot the chi2 statistic, which shows how FAR each origin's transitions
# moved. Dashed line = chi2 critical value at p = .05.
homog_row_10$origin = factor(homog_row_10$state, levels = rel_level_order)
crit_05 = qchisq(0.95, df = max(homog_row_10$df))

p_homog_row = ggplot(homog_row_10, aes(x = origin, y = chi2, fill = significant)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = crit_05, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  scale_fill_manual(
    values = c(`TRUE` = "#D55E00", `FALSE` = "#999999"),
    labels = c(`TRUE` = "Reject H0 (p < .05)", `FALSE` = "Fail to reject"),
    name   = NULL) +
  scale_x_discrete(labels = tools::toTitleCase) +
  labs(x = "Origin religion (RELIG16)", y = expression(chi^2 ~ "statistic"),
       subtitle = paste0("Dashed line: χ² critical value at p = .05 (df = ",
                         max(homog_row_10$df), ")"),
       title = "Which Origins Change Most Across Cohorts? Row-wise Homogeneity (10-year cohorts)") +
  healy_theme

ggsave("output/figures/homogeneity/rowwise_10yr.png", p_homog_row,
       width = 8, height = 5, dpi = 200)

# Pairwise-adjacent: p-value across cohort transitions localizes the break
pairwise_10$midpoint = (pairwise_10$from_cohort + pairwise_10$to_cohort) / 2

p_homog_pairwise = ggplot(pairwise_10, aes(x = midpoint, y = p_value)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.7, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Cohort transition (midpoint of adjacent 10-year bins)",
       y = "Joint homogeneity test p-value",
       title = "Where Does the Transition Matrix Change? Adjacent-Cohort Tests",
       caption = paste("Each point tests H0 that two consecutive 10-year cohort",
                       "matrices are equal (Anderson-Goodman 1957). Dashed line: p = 0.05.")) +
  healy_theme

ggsave("output/figures/homogeneity/pairwise_adjacent_10yr.png", p_homog_pairwise,
       width = 8, height = 5, dpi = 200)

# ── ROLLING-WINDOW SENSITIVITY (adapted from chi-square-test.R) ───────────────
# For several window widths, slide across birth-cohort time and run the joint
# homogeneity test on the 5-year-cohort matrices inside each window. A dip below
# p = 0.05 marks a neighborhood of cohort-time where the transition matrix is
# changing fast; comparing widths shows whether that localization is robust to
# the choice of window size. 5-year bins (not 1-year) keep each matrix well
# populated enough for the chi-square approximation to hold.

# Quiet joint statistic (chi2_joint without the console print)
joint_stat = function(mats) {
  rw = chi2_row(mats)
  jc = sum(rw$chi2); jd = sum(rw$df)
  c(chi2 = jc, df = jd, p = pchisq(jc, jd, lower.tail = FALSE))
}

roll_years = as.numeric(names(N_list_5))   # ordered 5-year cohort midpoints (non-integer, e.g. 1922.5)
bin_widths = c(2, 3, 4, 5)                 # window width in number of 5-year bins

roll_results = do.call(rbind, lapply(bin_widths, function(k) {
  half = floor(k / 2)
  # indices where a full window of k consecutive bins is available
  valid_idx = (half + 1):(length(roll_years) - (k - half - 1))
  do.call(rbind, lapply(valid_idx, function(i) {
    window_idx = (i - half):(i - half + k - 1)
    mats = N_list_5[window_idx]
    if (length(mats) < 2) return(NULL)
    jt = joint_stat(mats)
    data.frame(
      k          = k,
      center     = roll_years[i],
      year_min   = min(roll_years[window_idx]),
      year_max   = max(roll_years[window_idx]),
      n_matrices = length(mats),
      chi2       = jt[["chi2"]],
      p_value    = jt[["p"]],
      row.names  = NULL
    )
  }))
}))

roll_results$k_label = factor(paste0("k = ", roll_results$k, " bins"),
                              levels = paste0("k = ", bin_widths, " bins"))

p_roll = ggplot(roll_results, aes(x = center, y = p_value,
                                  color = k_label, shape = k_label, group = k_label)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.65, alpha = 0.85) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("indianred4", "goldenrod2", "darkgreen", "royalblue4"),
    name   = "Adjacent 5-yr\ncohort matrices") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Adjacent 5-yr\ncohort matrices") +
  scale_x_continuous(breaks = seq(1925, 1995, by = 5), guide = guide_axis(angle = 45)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    x       = "Center birth cohort",
    y       = "Joint homogeneity test p-value",
    caption = paste("Each point is a joint homogeneity test (Anderson-Goodman 1957) on the k",
                    "adjacent 5-year cohort matrices centered on that cohort (k = number of",
                    "matrices, not\ncalendar span). Dashed line: p = 0.05.")) +
  healy_theme +
  theme(legend.key.width = unit(1.2, "cm"),
        plot.caption     = element_text(size = 8, color = "grey40", hjust = 0))

ggsave("output/figures/homogeneity/rolling_window_sensitivity.png", p_roll,
       width = 12, height = 7, dpi = 300)

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

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
writeLines(homog_tex, "output/tables/homogeneity_tests.tex")
cat("\nWrote output/tables/homogeneity_tests.tex\n")
