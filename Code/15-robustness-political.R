# ── 15 · ROBUSTNESS: POLITICAL STRATIFICATION MATRICES ────────────────────────
# Decadal transition matrices stratified by party ID and political views,
# each in narrow and broad coding. Five cohort windows: 1940–1989.
#
# Input:  data/derived/gss_clean.rds
# Output: output/figures/political/*.png

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

pol_vars = list(
  partyid_narrow  = c("dem", "rep", "other"),
  partyid_broad   = c("dem", "rep", "other"),
  polviews_narrow = c("liberal", "moderate", "conservative"),
  polviews_broad  = c("liberal", "moderate", "conservative")
)

# 10-year bin midpoints (edges 1940–1980)
mids_pol = c(1945, 1955, 1965, 1975, 1985)

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
    title_str = paste0(lbl, " – Cohort ", edge, "–", edge + 9,
                       "  (N = ", n_pol[[key]], ")")
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
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
    labs(x = "Cohort (left edge of 10-year bin)",
         y = "Diagonal persistence P[i → i]",
         color = NULL,
         title = paste0("Diagonal Persistence: ",
                        gsub("_", " ", tools::toTitleCase(vname)))) +
    healy_theme

  ggsave(paste0("output/figures/political/diagonal_persistence_", vname, ".png"),
         p_diag, width = 12, height = 5, dpi = 200)
}

cat("Wrote output/figures/political/\n")
