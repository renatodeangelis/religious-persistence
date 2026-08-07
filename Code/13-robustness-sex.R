# ── 14 · ROBUSTNESS: SEX-STRATIFIED MATRICES ──────────────────────────────────
# Decadal transition matrices estimated separately for men and women.
# Six cohort windows: 1925–1984.
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

# 10-year bin midpoints (edges 1925–1975)
mids_sex   = c(1930, 1940, 1950, 1960, 1970, 1980)
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
    title_str = paste0(sex_lbl, " – Cohort ", edge, "–",
                       sprintf("%02d", (edge + 9) %% 100), "  (N = ", n_sex[[key]], ")")
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
  scale_x_continuous(breaks = c(1925, 1935, 1945, 1955, 1965, 1975),
                     labels = c("1925–34", "1935–44", "1945–54", "1955–64", "1965–74", "1975–84")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  labs(x = "Birth cohort (10-year bin)",
       y = "Probability to Stay",
       color = NULL,
       title = "Retention by Sex and Birth Cohort (1925–1984)") +
  healy_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/sex/diagonal_persistence_sex.png",
       p_diag, width = 10, height = 5, dpi = 200)

# ── TRANSITION MATRIX GRIDS (sex-stratified, 10-year cohorts) ────────────────
# One 6-panel grid per sex, saved separately for side-by-side appendix use.

for (sx in c("male", "female")) {
  grid_df_sx = do.call(rbind, lapply(mids_sex, function(mid) {
    key = paste(sx, mid, sep = "_")
    P   = P_sex[[key]]
    if (is.null(P)) return(NULL)
    df  = as.data.frame(as.table(P))
    names(df) = c("origin", "current", "p")
    df$cohort = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100),
                       "\n(N = ", format(n_sex[[key]], big.mark = ","), ")")
    df
  }))

  grid_df_sx$origin  = factor(grid_df_sx$origin,  levels = rel_level_order)
  grid_df_sx$current = factor(grid_df_sx$current, levels = rev(rel_level_order))
  cohort_lvls_sx = sapply(mids_sex, function(mid) {
    key = paste(sx, mid, sep = "_")
    if (is.null(P_sex[[key]])) return(NULL)
    paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100),
           "\n(N = ", format(n_sex[[key]], big.mark = ","), ")")
  })
  grid_df_sx$cohort = factor(grid_df_sx$cohort, levels = cohort_lvls_sx)

  p_grid_sx = ggplot(grid_df_sx, aes(current, origin, fill = p)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", p)), size = 2.8) +
    facet_wrap(~ cohort, nrow = 2) +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                         name = "P[i→j]") +
    labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
         title = paste0("Transition matrices by birth cohort — ",
                        tools::toTitleCase(sx))) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y      = element_text(size = 8),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = 9)
    )

  ggsave(paste0("output/figures/sex/P_grid_", sx, "_10yr.png"), p_grid_sx,
         width = 12, height = 9, dpi = 200)
}

saveRDS(
  list(P = P_sex, pi0 = pi0_sex, pistar = pistar_sex, n = n_sex),
  "data/derived/matrices_sex.rds"
)
cat("Wrote output/figures/sex/\n")
cat("Wrote data/derived/matrices_sex.rds\n")
