# ── 13 · ROBUSTNESS: NATIVITY-STRATIFIED MATRICES ─────────────────────────────
# Decadal transition matrices estimated separately for US-born and foreign-born
# respondents. Up to six cohort windows: 1925–1984. Early windows (esp. 1925–34,
# 1935–44) may be skipped for Born abroad due to thin cells (n < 30 guard).
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_nativity.rds
#         output/figures/nativity/*.png

source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt


# 10-year bin midpoints (edges 1925–1975); early windows may be skipped for Born abroad
mids_nat        = c(1930, 1940, 1950, 1960, 1970, 1980)
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
    title_str = paste0(nat_lbl, " – Cohort ", edge, "–",
                       sprintf("%02d", (edge + 9) %% 100), "  (N = ", n_nat[[key]], ")")
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
  scale_x_continuous(breaks = c(1925, 1935, 1945, 1955, 1965, 1975),
                     labels = c("1925–34", "1935–44", "1945–54", "1955–64", "1965–74", "1975–84")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  labs(x = "Birth cohort (10-year bin)",
       y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by Nativity and Birth Cohort (1925–1984)") +
  healy_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/nativity/diagonal_persistence_nativity.png",
       p_diag, width = 10, height = 5, dpi = 200)

saveRDS(
  list(P = P_nat, pi0 = pi0_nat, pistar = pistar_nat, n = n_nat),
  "data/derived/matrices_nativity.rds"
)
cat("Wrote output/figures/nativity/\n")
cat("Wrote data/derived/matrices_nativity.rds\n")

# ── π₀ AND π∞ GRID: 20-YEAR COHORTS × NATIVITY (6-STATE BP SCHEME) ──────────
# Three 20-year cohort windows (1925–44, 1945–64, 1965–84).
# Uses the 6-state BP scheme. Foreign-born Black Protestant cells may be thin
# (BP is a predominantly US-born category). If any foreign-born window drops
# due to n < 30, we stop — no fallback to the 5-state scheme.

states_bp = clean$states_bp
mids_20   = c(1935, 1955, 1975)
lbl_20    = c("1925–44", "1945–64", "1965–84")

data$cohort_20 = floor((data$cohort - 1925) / 20) * 20 + 1925 + 10
data$cohort_20[data$cohort < 1925 | data$cohort > 1984] = NA

P_nat20 = pi0_nat20 = pistar_nat20 = n_nat20 = list()

for (nat in nativity_groups) {
  for (mid in mids_20) {
    sub = data[!is.na(data$cohort_20) & data$cohort_20 == mid &
               !is.na(data$nativity)  & data$nativity  == nat &
               !is.na(data$reltrad16_bp) & !is.na(data$reltrad_bp), ]
    if (nrow(sub) < 30) next
    key                 = paste0(gsub(" ", "_", nat), "_", mid)
    P_nat20[[key]]      = p_matrix(sub, "reltrad16_bp", "reltrad_bp", levels = states_bp)
    pi0_nat20[[key]]    = pi_0(sub, "reltrad16_bp")
    pistar_nat20[[key]] = pi_star(P_nat20[[key]])
    n_nat20[[key]]      = nrow(sub)
  }
}

# Thin-cell guard for foreign-born — do not fall back to 5-state
fb_expected = paste0("Born_abroad_", mids_20)
fb_missing  = setdiff(fb_expected, names(P_nat20))
if (length(fb_missing) > 0) {
  # NOTE: Foreign-born 20-year cohort(s) dropped due to n < 30 in the BP scheme.
  # Missing windows stored in fb_missing. Per design, no fallback to 5-state.
  stop("Thin cells in foreign-born × BP scheme. Missing windows: ",
       paste(fb_missing, collapse = ", "))
}

# Long-format data for paired bar figure
pi_nat20_df = do.call(rbind, lapply(names(P_nat20), function(key) {
  mid_val    = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  nat_lbl    = gsub("_", " ", sub("_(\\d{4})$", "", key))
  cohort_lbl = lbl_20[match(mid_val, mids_20)]
  pi0v       = pi0_nat20[[key]]
  piv        = pistar_nat20[[key]]
  rbind(
    data.frame(nativity = nat_lbl, cohort = cohort_lbl,
               religion = names(pi0v), value = as.numeric(pi0v), measure = "π₀"),
    data.frame(nativity = nat_lbl, cohort = cohort_lbl,
               religion = names(piv),  value = as.numeric(piv),  measure = "π∞")
  )
}))

pi_nat20_df$religion = factor(pi_nat20_df$religion, levels = rel_level_order)
pi_nat20_df$cohort   = factor(pi_nat20_df$cohort,   levels = lbl_20)
pi_nat20_df$measure  = factor(pi_nat20_df$measure,  levels = c("π₀", "π∞"))
pi_nat20_df$nativity = factor(pi_nat20_df$nativity, levels = nativity_groups)

p_pi_nat20 = ggplot(pi_nat20_df, aes(x = religion, y = value)) +
  ggpattern::geom_col_pattern(
    aes(fill = religion, pattern = measure),
    position             = position_dodge(width = 0.8), width = 0.8,
    color                = "grey25", linewidth = 0.25,
    pattern_fill         = "white", pattern_colour = "white",
    pattern_angle        = 45,      pattern_density = 0.08,
    pattern_spacing      = 0.03,    pattern_key_scale_factor = 0.6
  ) +
  scale_fill_manual(values = reltrad_colors, labels = reltrad_labels_tc, guide = "none") +
  ggpattern::scale_pattern_manual(
    values = c("π₀" = "stripe", "π∞" = "none"),
    labels = c("π₀" = expression(Origin~(pi[0])),
               "π∞" = expression(Stationary~(pi[infinity])))
  ) +
  guides(pattern = guide_legend(override.aes = list(fill = "grey55", color = "grey25"))) +
  scale_x_discrete(labels = reltrad_labels_tc) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  facet_grid(nativity ~ cohort) +
  labs(x = NULL, y = "Share") +
  theme_bc(base_size = 12, x_angle = 45) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("output/figures/nativity/pi_dist_20yr_nativity.png", p_pi_nat20,
       width = 10, height = 6, dpi = 200)
cat("Wrote output/figures/nativity/pi_dist_20yr_nativity.png\n")
