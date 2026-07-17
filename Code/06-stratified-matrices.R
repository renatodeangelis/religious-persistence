# ── 06 · STRATIFIED TRANSITION MATRICES ────────────────────────────────────────
# Console tables and figures for every stratified matrix set built in 02:
# binary affiliated/unaffiliated 2×2, nativity, sex, and political (party ID /
# polviews).
#
# Input:  data/derived/matrices.rds, data/derived/gss_clean.rds
# Output: output/figures/{binary,nativity,sex,political}/*.png

library(patchwork)
source("code/utils.R")

matrices   = readRDS("data/derived/matrices.rds")
clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

P_list_2x2      = matrices$binary$P
pi0_list_2x2    = matrices$binary$pi0
pistar_list_2x2 = matrices$binary$pistar
n_list_2x2      = matrices$binary$n
states_2x2      = matrices$binary$states

P_list_nat      = matrices$nativity$P
pi0_list_nat    = matrices$nativity$pi0
pistar_list_nat = matrices$nativity$pistar
n_list_nat      = matrices$nativity$n

P_list_sex      = matrices$sex$P
pi0_list_sex    = matrices$sex$pi0
pistar_list_sex = matrices$sex$pistar
n_list_sex      = matrices$sex$n

P_list_pol      = matrices$political$P
pi0_list_pol    = matrices$political$pi0
pistar_list_pol = matrices$political$pistar
n_list_pol      = matrices$political$n

# ── BINARY (AFFILIATED / UNAFFILIATED) 2×2 MATRICES ──────────────────────────

# Print matrices with N
for (key in names(P_list_2x2)) {
  edge = as.numeric(key) - 5   # 10-yr bin midpoint → left edge
  cat("\n── Cohort", edge, "–", edge + 9, "  (N =", n_list_2x2[[key]], ") ──\n")
  print(round(P_list_2x2[[key]], 3))
}

# ── FIGURES: individual make_combined heatmaps per cohort ────────────────────

dir.create("output/figures/binary", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_2x2)) {
  p = make_combined(
    P_list_2x2[[key]], pi0_list_2x2[[key]], pistar_list_2x2[[key]],
    levels    = states_2x2,
    title_str = paste0("Cohort ", as.numeric(key) - 5, "–", as.numeric(key) + 4,
                       "  (N = ", n_list_2x2[[key]], ")")
  )
  ggsave(paste0("output/figures/binary/trans_", as.numeric(key) - 5, "_10yr_2x2.png"),
         p, width = 7, height = 5, dpi = 200)
}

# ── FIGURE: all 4 cells over cohorts ─────────────────────────────────────────

cells_2x2 = do.call(rbind, lapply(names(P_list_2x2), function(key) {
  P   = P_list_2x2[[key]]
  df  = as.data.frame(as.table(P), stringsAsFactors = FALSE)
  names(df) = c("origin", "current", "prob")
  df$cohort = as.numeric(key)   # 10-yr bin midpoint, used as continuous x
  df
}))

p_cells_2x2 = ggplot(cells_2x2,
    aes(x = cohort, y = prob, color = current, group = current)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ origin, nrow = 1,
             labeller = labeller(origin = c(
               affiliated   = "Origin: Affiliated",
               unaffiliated = "Origin: Unaffiliated"))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    values = c(affiliated = "#0072B2", unaffiliated = "#999999"),
    labels = c(affiliated = "→ Affiliated", unaffiliated = "→ Unaffiliated")) +
  labs(x = "Birth cohort (10-year bins)", y = "Transition probability",
       color = NULL,
       title = "2×2 Transition Probabilities by Birth Cohort") +
  healy_theme

ggsave("output/figures/binary/cells_2x2_10yr.png", p_cells_2x2,
       width = 10, height = 5, dpi = 200)

# ── FIGURE: diagonal persistence over cohorts ─────────────────────────────────

persistence_2x2 = do.call(rbind, lapply(names(P_list_2x2), function(key) {
  P = P_list_2x2[[key]]
  data.frame(cohort = as.numeric(key), state = rownames(P),
             persistence = diag(P), row.names = NULL)
}))

p_persist_2x2 = ggplot(persistence_2x2,
    aes(x = cohort, y = persistence, color = state, group = state)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_color_manual(
    values = c(affiliated = "#0072B2", unaffiliated = "#999999"),
    labels = c(affiliated = "Affiliated", unaffiliated = "Unaffiliated")) +
  labs(x = "Birth cohort (10-year bins)",
       y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence: Affiliated vs. Unaffiliated (10-year cohorts)") +
  healy_theme

ggsave("output/figures/binary/persistence_2x2_10yr.png", p_persist_2x2,
       width = 8, height = 5, dpi = 200)

# ── NATIVITY-SPLIT MATRICES (10-year cohorts: 1950, 1960, 1970) ──────────────

for (key in names(P_list_nat)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5   # midpoint → left edge
  cat("\n── ", gsub("_", " ", sub("_(\\d{4})$", "", key)),
      " | Cohort", edge, "–", edge + 9,
      "  (N =", n_list_nat[[key]], ") ──\n")
  print(round(P_list_nat[[key]], 3))
}

dir.create("output/figures/nativity", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_nat)) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5   # midpoint → left edge
  nat_lbl = gsub("_", " ", sub("_(\\d{4})$", "", key))
  p = make_combined(
    P_list_nat[[key]], pi0_list_nat[[key]], pistar_list_nat[[key]],
    levels    = rel_level_order,
    title_str = paste0(nat_lbl, " – ", edge, "–", edge + 9,
                       "  (N = ", n_list_nat[[key]], ")")
  )
  ggsave(paste0("output/figures/nativity/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}

# ── SEX-STRATIFIED DECADAL MATRICES (10-year cohorts, 1940–1980) ─────────────

for (key in names(P_list_sex)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5   # midpoint → left edge
  cat("\n──", tools::toTitleCase(sub("_(\\d{4})$", "", key)),
      "| Cohort", edge, "–", edge + 9,
      " (N =", n_list_sex[[key]], ") ──\n")
  print(round(P_list_sex[[key]], 3))
}

dir.create("output/figures/sex", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_sex)) {
  edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5   # midpoint → left edge
  sex_lbl = tools::toTitleCase(sub("_(\\d{4})$", "", key))
  p = make_combined(
    P_list_sex[[key]], pi0_list_sex[[key]], pistar_list_sex[[key]],
    levels    = rel_level_order,
    title_str = paste0(sex_lbl, " – ", edge, "–", edge + 9,
                       "  (N = ", n_list_sex[[key]], ")")
  )
  ggsave(paste0("output/figures/sex/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}

# ── POLITICAL STRATIFICATION DECADAL MATRICES (10-year cohorts, 1940–1989) ───

for (key in names(P_list_pol)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5   # midpoint → left edge
  lbl  = gsub("_", " ", sub("_\\d{4}$", "", key))
  cat("\n──", lbl, "| Cohort", edge, "–", edge + 9,
      " (N =", n_list_pol[[key]], ") ──\n")
  print(round(P_list_pol[[key]], 3))
}

dir.create("output/figures/political", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_pol)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5   # midpoint → left edge
  lbl  = gsub("_", " ", sub("_\\d{4}$", "", key))
  p = make_combined(
    P_list_pol[[key]], pi0_list_pol[[key]], pistar_list_pol[[key]],
    levels    = rel_level_order,
    title_str = paste0(lbl, " – ", edge, "–", edge + 9,
                       "  (N = ", n_list_pol[[key]], ")")
  )
  ggsave(paste0("output/figures/political/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}
