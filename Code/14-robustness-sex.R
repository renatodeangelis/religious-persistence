# ── 14 · ROBUSTNESS: SEX-STRATIFIED MATRICES ──────────────────────────────────
# Decadal transition matrices estimated separately for men and women.
# Five cohort windows: 1940–1949 through 1980–1989.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_sex.rds
#         output/figures/sex/*.png

source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")
reltrad_colors  = c(
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

# 10-year bin midpoints (edges 1940–1980)
mids_sex   = c(1940, 1950, 1960, 1970, 1980)
sex_labels = c("1" = "male", "2" = "female")

# ── BUILD MATRICES ────────────────────────────────────────────────────────────

P_sex = pi0_sex = pistar_sex = n_sex = list()

for (sx in c(1, 2)) {
  for (mid in mids_sex) {
    sub = data[!is.na(data$cohort_10)     & data$cohort_10         == mid  &
               !is.na(data$sex)           & as.numeric(data$sex)   == sx   &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
    if (nrow(sub) < 30) next
    key = paste(sex_labels[as.character(sx)], mid, sep = "_")
    P_sex[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0_sex[[key]]    = pi_0(sub, "reltrad16_alt")
    pistar_sex[[key]] = pi_star(P_sex[[key]])
    n_sex[[key]]      = nrow(sub)
  }
}

# ── CONSOLE OUTPUT ───────────────────────────────────────────────────────────

for (key in names(P_sex)) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  sex_lbl = tools::toTitleCase(sub("_(\\d{4})$", "", key))
  cat("\n──", sex_lbl, "| Cohort", edge, "–", edge + 9,
      " (N =", n_sex[[key]], ") ──\n")
  print(round(P_sex[[key]], 3))
}

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/sex", recursive = TRUE, showWarnings = FALSE)

# Individual heatmaps
for (key in names(P_sex)) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  sex_lbl = tools::toTitleCase(sub("_(\\d{4})$", "", key))
  p = make_combined(
    P_sex[[key]], pi0_sex[[key]], pistar_sex[[key]],
    levels    = rel_level_order,
    title_str = paste0(sex_lbl, " – Cohort ", edge, "–", edge + 9,
                       "  (N = ", n_sex[[key]], ")")
  )
  ggsave(paste0("output/figures/sex/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}

# Diagonal persistence by sex and cohort
diag_sex = do.call(rbind, lapply(names(P_sex), function(key) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  sex_lbl = tools::toTitleCase(sub("_(\\d{4})$", "", key))
  data.frame(cohort  = edge,
             sex     = sex_lbl,
             origin  = rel_level_order,
             persist = diag(P_sex[[key]])[rel_level_order],
             row.names = NULL)
}))
diag_sex$origin = factor(diag_sex$origin, levels = rel_level_order)

p_diag = ggplot(diag_sex, aes(x = cohort, y = persist, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ sex) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  labs(x = "Cohort (left edge of 10-year bin)",
       y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by Sex and Birth Cohort") +
  healy_theme

ggsave("output/figures/sex/diagonal_persistence_sex.png",
       p_diag, width = 10, height = 5, dpi = 200)

saveRDS(
  list(P = P_sex, pi0 = pi0_sex, pistar = pistar_sex, n = n_sex),
  "data/derived/matrices_sex.rds"
)
cat("Wrote output/figures/sex/\n")
cat("Wrote data/derived/matrices_sex.rds\n")
