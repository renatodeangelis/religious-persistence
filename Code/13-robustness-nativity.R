# ── 13 · ROBUSTNESS: NATIVITY-STRATIFIED MATRICES ─────────────────────────────
# Decadal transition matrices estimated separately for US-born and foreign-born
# respondents. Three cohort windows: 1950–1959, 1960–1969, 1970–1979.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_nativity.rds
#         output/figures/nativity/*.png

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

# 10-year bin midpoints (edges 1950/1960/1970)
mids_nat        = c(1950, 1960, 1970)
nativity_groups = c("Born in US", "Born abroad")

# ── BUILD MATRICES ────────────────────────────────────────────────────────────

P_nat = pi0_nat = pistar_nat = n_nat = list()

for (nat in nativity_groups) {
  for (mid in mids_nat) {
    sub = data[!is.na(data$cohort_10) & data$cohort_10 == mid &
               !is.na(data$nativity)  & data$nativity   == nat &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
    if (nrow(sub) < 30) next
    key = paste0(gsub(" ", "_", nat), "_", mid)
    P_nat[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0_nat[[key]]    = pi_0(sub, "reltrad16_alt")
    pistar_nat[[key]] = pi_star(P_nat[[key]])
    n_nat[[key]]      = nrow(sub)
  }
}

# ── CONSOLE OUTPUT ───────────────────────────────────────────────────────────

for (key in names(P_nat)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  nat  = gsub("_", " ", sub("_(\\d{4})$", "", key))
  cat("\n──", nat, "| Cohort", edge, "–", edge + 9,
      "  (N =", n_nat[[key]], ") ──\n")
  print(round(P_nat[[key]], 3))
}

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/nativity", recursive = TRUE, showWarnings = FALSE)

# Individual heatmaps
for (key in names(P_nat)) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  nat_lbl = gsub("_", " ", sub("_(\\d{4})$", "", key))
  p = make_combined(
    P_nat[[key]], pi0_nat[[key]], pistar_nat[[key]],
    levels    = rel_level_order,
    title_str = paste0(nat_lbl, " – Cohort ", edge, "–", edge + 9,
                       "  (N = ", n_nat[[key]], ")")
  )
  ggsave(paste0("output/figures/nativity/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}

# Diagonal persistence by nativity group and cohort
diag_nat = do.call(rbind, lapply(names(P_nat), function(key) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  nat_lbl = gsub("_", " ", sub("_(\\d{4})$", "", key))
  data.frame(cohort   = edge,
             nativity = nat_lbl,
             origin   = rel_level_order,
             persist  = diag(P_nat[[key]])[rel_level_order],
             row.names = NULL)
}))
diag_nat$origin = factor(diag_nat$origin, levels = rel_level_order)

p_diag = ggplot(diag_nat, aes(x = cohort, y = persist, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ nativity) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  labs(x = "Cohort (left edge of 10-year bin)",
       y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by Nativity and Birth Cohort") +
  healy_theme

ggsave("output/figures/nativity/diagonal_persistence_nativity.png",
       p_diag, width = 10, height = 5, dpi = 200)

saveRDS(
  list(P = P_nat, pi0 = pi0_nat, pistar = pistar_nat, n = n_nat),
  "data/derived/matrices_nativity.rds"
)
cat("Wrote output/figures/nativity/\n")
cat("Wrote data/derived/matrices_nativity.rds\n")
