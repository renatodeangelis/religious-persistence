# ── 15 · ROBUSTNESS: POLITICAL STRATIFICATION MATRICES ────────────────────────
# Decadal transition matrices stratified by party ID and political views,
# each in narrow and broad coding. Six cohort windows: 1925–1984.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_political.rds
#         output/figures/political/*.png

source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt


pol_vars = list(
  partyid_narrow  = c("dem", "rep", "other"),
  partyid_broad   = c("dem", "rep", "other"),
  polviews_narrow = c("liberal", "moderate", "conservative"),
  polviews_broad  = c("liberal", "moderate", "conservative")
)

# 10-year bin midpoints (edges 1925–1975)
mids_pol = c(1930, 1940, 1950, 1960, 1970, 1980)

# ── BUILD MATRICES ────────────────────────────────────────────────────────────

P_pol = pi0_pol = pistar_pol = n_pol = list()

for (vname in names(pol_vars)) {
  for (grp in pol_vars[[vname]]) {
    for (mid in mids_pol) {
      sub = data[!is.na(data$cohort_10)     & data$cohort_10 == mid  &
                 !is.na(data[[vname]])      & data[[vname]]  == grp  &
                 !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
      if (nrow(sub) < 30) next
      key = paste(vname, grp, mid, sep = "_")
      P_pol[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
      pi0_pol[[key]]    = pi_0(sub, "reltrad16_alt")
      pistar_pol[[key]] = pi_star(P_pol[[key]])
      n_pol[[key]]      = nrow(sub)
    }
  }
}

# ── CONSOLE OUTPUT ───────────────────────────────────────────────────────────

for (key in names(P_pol)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  lbl  = gsub("_", " ", sub("_\\d{4}$", "", key))
  cat("\n──", lbl, "| Cohort", edge, "–", edge + 9,
      " (N =", n_pol[[key]], ") ──\n")
  print(round(P_pol[[key]], 3))
}

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/political", recursive = TRUE, showWarnings = FALSE)

# Individual heatmaps
for (key in names(P_pol)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  lbl  = gsub("_", " ", sub("_\\d{4}$", "", key))
  p = make_combined(
    P_pol[[key]], pi0_pol[[key]], pistar_pol[[key]],
    levels    = rel_level_order,
    title_str = paste0(lbl, " – Cohort ", edge, "–",
                       sprintf("%02d", (edge + 9) %% 100), "  (N = ", n_pol[[key]], ")")
  )
  ggsave(paste0("output/figures/political/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}

# Diagonal persistence — one figure per variable (faceted by group)
for (vname in names(pol_vars)) {
  keys_v = grep(paste0("^", vname, "_"), names(P_pol), value = TRUE)
  if (length(keys_v) == 0) next

  diag_v = do.call(rbind, lapply(keys_v, function(key) {
    edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
    grp  = sub(paste0("^", vname, "_(.*)_\\d{4}$"), "\\1", key)
    data.frame(cohort  = edge,
               group   = tools::toTitleCase(grp),
               origin  = rel_level_order,
               persist = diag(P_pol[[key]])[rel_level_order],
               row.names = NULL)
  }))
  diag_v$origin = factor(diag_v$origin, levels = rel_level_order)

  p_diag = ggplot(diag_v, aes(x = cohort, y = persist, color = origin, group = origin)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ group) +
    scale_x_continuous(breaks = c(1925, 1935, 1945, 1955, 1965, 1975),
                       labels = c("1925–34", "1935–44", "1945–54", "1955–64", "1965–74", "1975–84")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
    labs(x = "Birth cohort (10-year bin)",
         y = "Probability to Stay",
         color = NULL,
         title = paste0("Retention: ",
                        gsub("_", " ", tools::toTitleCase(vname)), " (1925–1984)")) +
    healy_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(paste0("output/figures/political/diagonal_persistence_", vname, ".png"),
         p_diag, width = 12, height = 5, dpi = 200)
}

# ── TRANSITION MATRIX GRIDS (polviews_broad, 10-year cohorts) ────────────────
# One 6-panel grid per stratum (liberal / moderate / conservative).
# Flag: add matching loop for partyid_broad if wanted — swap "polviews_broad"
#       and pol_vars[["partyid_broad"]] below and update the output filename.

for (grp in pol_vars[["polviews_broad"]]) {
  grid_df_pol = do.call(rbind, lapply(mids_pol, function(mid) {
    key = paste("polviews_broad", grp, mid, sep = "_")
    P   = P_pol[[key]]
    if (is.null(P)) return(NULL)
    df  = as.data.frame(as.table(P))
    names(df) = c("origin", "current", "p")
    df$cohort = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100),
                       "\n(N = ", format(n_pol[[key]], big.mark = ","), ")")
    df
  }))

  grid_df_pol$origin  = factor(grid_df_pol$origin,  levels = rel_level_order)
  grid_df_pol$current = factor(grid_df_pol$current, levels = rev(rel_level_order))
  cohort_lvls_pol = sapply(mids_pol, function(mid) {
    key = paste("polviews_broad", grp, mid, sep = "_")
    if (is.null(P_pol[[key]])) return(NULL)
    paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100),
           "\n(N = ", format(n_pol[[key]], big.mark = ","), ")")
  })
  grid_df_pol$cohort = factor(grid_df_pol$cohort, levels = cohort_lvls_pol)

  p_grid_pol = ggplot(grid_df_pol, aes(current, origin, fill = p)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", p)), size = 2.8) +
    facet_wrap(~ cohort, nrow = 2) +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                         name = "P[i→j]") +
    labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
         title = paste0("Transition matrices by birth cohort — ",
                        tools::toTitleCase(grp), " (polviews broad)")) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y      = element_text(size = 8),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = 9)
    )

  ggsave(paste0("output/figures/political/P_grid_polviews_broad_", grp, "_10yr.png"),
         p_grid_pol, width = 12, height = 9, dpi = 200)
}

saveRDS(
  list(P = P_pol, pi0 = pi0_pol, pistar = pistar_pol, n = n_pol),
  "data/derived/matrices_political.rds"
)
cat("Wrote output/figures/political/\n")
cat("Wrote data/derived/matrices_political.rds\n")
