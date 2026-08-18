# ── 12 · ROBUSTNESS: 6-STATE SPACE (BLACK PROTESTANT SEPARATE) ───────────────
# Re-estimates cohort transition matrices on the 6-state space that keeps Black
# Protestant as its own category (only jewish → other collapsed), at both 20-
# and 10-year resolution. Reports transition heatmaps, IM memory curves, and
# diagonal persistence. Uses the reltrad_bp/reltrad16_bp columns built in 01.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_bp.rds
#         output/figures/bp/*.png

library(dplyr)
library(ggplot2)
source("code/utils.R")

data(gss_all)
data = gss_all |>
  select(year, cohort, sex, reltrad, reltrad16, region, born,
         race, polviews, partyid, sibs_7222, childs) |>
  filter(!(year %in% c(1972, 2021))) |>
  mutate(across(c(reltrad, reltrad16),
                ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  # 6-state variant for the Black-Protestant robustness stage (12): collapse
  # only jewish → other, keeping "black protestant" as its own state.
  mutate(across(c(reltrad, reltrad16),
                ~ if_else(. == "jewish", "other", .),
                .names = "{.col}_bp")) |>
  # cohort arrives from gss_all as a haven_labelled vector; strip to plain
  # numeric so downstream median()/binning behave (median.haven_labelled errors)
  mutate(cohort = as.numeric(cohort)) |>
  mutate(age = year - cohort) |>
  filter(age >= 30, age <= 75, cohort >= 1925, cohort <= 1984) |>
  mutate(
    cohort_10  = (floor((cohort - 1925) / 10) * 10 + 1925) + 5,
    cohort_5 = (floor((cohort - 1925) / 5) * 5 + 1925) + 2.5,
    region_broad = case_when(
      as.numeric(region) == 1 ~ "Northeast",
      as.numeric(region) == 2 ~ "Midwest",
      as.numeric(region) == 3 ~ "South",
      as.numeric(region) == 4 ~ "West",
      TRUE ~ NA_character_
    ),
    nativity = case_when(
      as.numeric(born) == 1 ~ "Born in US",
      as.numeric(born) == 2 ~ "Born abroad"
    )
  )

# ── PARTY ID AND POLITICAL VIEWS RECODES ─────────────────────────────────────
# partyid: 0 = strong dem … 6 = strong rep, 7 = other party
# polviews: 1 = extremely liberal … 7 = extremely conservative

data = data |>
  mutate(
    partyid_narrow = case_when(
      as.numeric(partyid) %in% 0:1            ~ "dem",
      as.numeric(partyid) %in% 5:6            ~ "rep",
      as.numeric(partyid) %in% c(2, 3, 4, 7) ~ "other"
    ),
    partyid_broad = case_when(
      as.numeric(partyid) %in% 0:2        ~ "dem",
      as.numeric(partyid) %in% 4:6        ~ "rep",
      as.numeric(partyid) %in% c(3, 7)    ~ "other"
    ),
    polviews_narrow = case_when(
      as.numeric(polviews) %in% 1:3 ~ "liberal",
      as.numeric(polviews) == 4      ~ "moderate",
      as.numeric(polviews) %in% 5:7 ~ "conservative"
    ),
    polviews_broad = case_when(
      as.numeric(polviews) %in% 1:2 ~ "liberal",
      as.numeric(polviews) %in% 3:5 ~ "moderate",
      as.numeric(polviews) %in% 6:7 ~ "conservative"
    )
  )

# ── STATE SPACE ───────────────────────────────────────────────────────────────

# 6-state space for the Black-Protestant robustness stage (12)
states_bp = sort(unique(c(data$reltrad_bp, data$reltrad16_bp)))
states_bp = states_bp[!is.na(states_bp)]

# ── SAVE ──────────────────────────────────────────────────────────────────────
# Strip haven value labels so the saved frame is plain numeric/character. The
# gssr columns arrive as haven_labelled, and as.numeric() on those only works
# while haven's S3 methods are attached (they are here, via library(gssr), but
# not in the downstream scripts that read this file).

data = haven::zap_labels(data)


# Fixed 6-state display order (Black Protestant retained; jewish already → other)
rel_level_order_6 = c("catholic", "evangelical", "black protestant", "mainline", "other", "none")
reltrad_colors_6 = c(
  catholic           = "#0072B2",
  evangelical        = "#D55E00",
  "black protestant" = "#E69F00",
  mainline           = "#009E73",
  other              = "#CC79A7",
  none               = "#999999"
)
reltrad_labels_6 = c(
  catholic           = "Catholic",
  evangelical        = "Evangelical",
  "black protestant" = "Black Protestant",
  mainline           = "Mainline",
  other              = "Other",
  none               = "None"
)

stopifnot(setequal(states_bp, rel_level_order_6))

# ── HELPER: build P/pi0/pistar lists for a set of bin midpoints ───────────────
# bin_col is cohort_10 (pipeline midpoints); `offset` (5) converts the midpoint
# to its bin edge. Matrices are stored under the edge label so titles/filenames
# match the edge convention.

build_bp = function(bin_col, mids, offset) {
  P = pi0 = pistar = nn = list()
  for (mid in mids) {
    sub = data[!is.na(data[[bin_col]]) & data[[bin_col]] == mid &
               !is.na(data$reltrad16_bp) & !is.na(data$reltrad_bp), ]
    if (nrow(sub) < 30) next
    key = as.character(mid - offset)
    P[[key]]      = p_matrix(sub, "reltrad16_bp", "reltrad_bp", levels = rel_level_order_6)
    pi0[[key]]    = pi_0(sub, "reltrad16_bp")
    pistar[[key]] = pi_star(P[[key]])
    nn[[key]]     = nrow(sub)
  }
  list(P = P, pi0 = pi0, pistar = pistar, n = nn)
}

# ── CELL-COUNT DIAGNOSTIC ─────────────────────────────────────────────────────
# Black Protestant rows are thin in early cohorts — inspect before trusting.

cell_diag = function(bin_col, mids, offset, label) {
  cat("\n──", label, "──\n")
  out = lapply(mids, function(mid) {
    sub = data[!is.na(data[[bin_col]]) & data[[bin_col]] == mid &
               !is.na(data$reltrad16_bp) & !is.na(data$reltrad_bp), ]
    if (nrow(sub) == 0) return(NULL)
    tab = table(factor(sub$reltrad16_bp, levels = rel_level_order_6),
                factor(sub$reltrad_bp,   levels = rel_level_order_6))
    data.frame(cohort = mid - offset, n = nrow(sub),
               n_bp_orig = sum(tab["black protestant", ]),
               min_cell = min(tab), cells_lt5 = sum(tab < 5),
               cells_0 = sum(tab == 0), row.names = NULL)
  })
  print(do.call(rbind, out), row.names = FALSE)
}

# ── IM MEMORY CURVE DATA ──────────────────────────────────────────────────────

im_data = function(lst) {
  df = do.call(rbind, lapply(names(lst$P), function(key) {
    do.call(rbind, lapply(0:4, function(t) {
      vals = im_from_P(lst$P[[key]], t = t)
      data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
                 row.names = NULL)
    }))
  }))
  df$origin = factor(df$origin, levels = rel_level_order_6)
  df
}

# ── DIAGONAL PERSISTENCE DATA ─────────────────────────────────────────────────

diag_data = function(lst) {
  df = do.call(rbind, lapply(names(lst$P), function(key) {
    P = lst$P[[key]]
    data.frame(cohort = as.integer(key), origin = rownames(P),
               persistence = diag(P), row.names = NULL)
  }))
  df$origin = factor(df$origin, levels = rel_level_order_6)
  df
}

dir.create("output/figures/bp", recursive = TRUE, showWarnings = FALSE)

# ── 10-YEAR COHORTS (edges 1925–1975) ────────────────────────────────────────

mids_10 = c(1930, 1940, 1950, 1960, 1970, 1980)   # edges 1925–1975
cell_diag("cohort_10", mids_10, 5, "Cell-count diagnostic (6-state, 10-year cohorts, 1925–1984)")
bp10 = build_bp("cohort_10", mids_10, 5)

for (key in names(bp10$P)) {
  cat("\n── Cohort", key, "–", as.integer(key) + 9, "  (N =", bp10$n[[key]], ") ──\n")
  print(round(bp10$P[[key]], 3))
}

for (key in names(bp10$P)) {
  p = make_combined(
    bp10$P[[key]], bp10$pi0[[key]], bp10$pistar[[key]],
    levels    = rel_level_order_6,
    title_str = paste0("6-State — Cohort ", key, "–",
                       sprintf("%02d", (as.integer(key) + 9) %% 100),
                       "  (N = ", bp10$n[[key]], ")")
  )
  ggsave(paste0("output/figures/bp/trans_", key, "_10yr_6state.png"),
         p, width = 11, height = 7, dpi = 200)
}

im_df_bp10 = im_data(bp10)

p_im_bp10 = ggplot(im_df_bp10, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ cohort, nrow = 2,
             labeller = labeller(cohort = function(x) paste0(x, "–", sprintf("%02d", (as.integer(x) + 9) %% 100)))) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Step (t)", y = "log(TV distance from π∞)", color = NULL,
       title = "Individual Memory by Cohort — 6-State Space (10-year bins, 1925–1984)") +
  healy_theme

ggsave("output/figures/bp/im_memory_10yr_6state.png", p_im_bp10, width = 10, height = 8, dpi = 200)

diag_df10 = diag_data(bp10)

p_diag10 = ggplot(diag_df10, aes(x = cohort, y = persistence, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = c(1925, 1935, 1945, 1955, 1965, 1975),
                     labels = c("1925–34", "1935–44", "1945–54", "1955–64", "1965–74", "1975–84")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Birth cohort (10-year bins)", y = "Diagonal persistence P[i → i]", color = NULL,
       title = "Diagonal Persistence by Origin — 6-State Space (10-year cohorts, 1925–1984)") +
  healy_theme

ggsave("output/figures/bp/diagonal_persistence_10yr_6state.png",
       p_diag10, width = 8, height = 5, dpi = 200)

# ── TRANSITION MATRIX GRID (6-state, 10-year cohorts, 1925–1984) ──────────────

edges_bp   = names(bp10$P)   # edge labels: "1925", "1935", ..., "1975"
n_per_edge = sapply(edges_bp, function(k) format(bp10$n[[k]], big.mark = ","))

grid_df_bp = do.call(rbind, lapply(seq_along(edges_bp), function(i) {
  key = edges_bp[i]
  P   = bp10$P[[key]]
  if (is.null(P)) return(NULL)
  df  = as.data.frame(as.table(P))
  names(df) = c("origin", "current", "p")
  df$cohort = paste0(key, "–", sprintf("%02d", (as.integer(key) + 9) %% 100),
                     "\n(N = ", n_per_edge[i], ")")
  df
}))

# Match plot_pmat_heatmap orientation: origin = levels (catholic at bottom), current = rev(levels)
grid_df_bp$origin  = factor(grid_df_bp$origin,  levels = rel_level_order_6)
grid_df_bp$current = factor(grid_df_bp$current, levels = rev(rel_level_order_6))
grid_df_bp$cohort  = factor(grid_df_bp$cohort,
  levels = paste0(edges_bp, "–", sprintf("%02d", (as.integer(edges_bp) + 9) %% 100),
                  "\n(N = ", n_per_edge, ")"))

p_grid_bp = ggplot(grid_df_bp, aes(current, origin, fill = p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p)), size = 2.5) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1), name = "P[i→j]") +
  labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
       title = "Intergenerational transition matrices by birth cohort — 6-State Space") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 9),
    plot.caption     = element_text(size = 8, color = "grey50")
  )

ggsave("output/figures/bp/P_grid_bp_10yr.png", p_grid_bp,
       width = 14, height = 10, dpi = 200)

# ── π₀ AND π∞ DISTRIBUTION GRID (6-state, 10-year cohorts, 1925–1984) ─────────

dist_df_bp = do.call(rbind, lapply(edges_bp, function(key) {
  pi0    = bp10$pi0[[key]]
  pistar = bp10$pistar[[key]]
  if (is.null(pi0)) return(NULL)
  cohort_lbl = paste0(key, "–", sprintf("%02d", (as.integer(key) + 9) %% 100))
  rbind(
    data.frame(cohort = cohort_lbl, measure = "π₀  (origin)",
               religion = names(pi0),    value = as.numeric(pi0)),
    data.frame(cohort = cohort_lbl, measure = "π∞ (stationary)",
               religion = names(pistar), value = as.numeric(pistar))
  )
}))

dist_df_bp$religion = factor(dist_df_bp$religion, levels = rel_level_order_6)
dist_df_bp$cohort   = factor(dist_df_bp$cohort,
  levels = paste0(edges_bp, "–", sprintf("%02d", (as.integer(edges_bp) + 9) %% 100)))
dist_df_bp$measure  = factor(dist_df_bp$measure,
  levels = c("π₀  (origin)", "π∞ (stationary)"))

p_dist_grid_bp = ggplot(dist_df_bp, aes(y = religion, x = value, fill = religion)) +
  geom_col(width = 0.65) +
  facet_grid(measure ~ cohort) +
  scale_fill_manual(values = reltrad_colors_6, labels = reltrad_labels_6, name = NULL) +
  scale_x_continuous(limits = c(0, 0.65), breaks = c(0, 0.25, 0.5),
                     labels = c("0", ".25", ".5")) +
  scale_y_discrete(labels = reltrad_labels_6) +
  labs(x = "Share", y = NULL,
       title = "Origin (π₀) and stationary (π∞) distributions by birth cohort — 6-State Space") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "none",
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(size = 8),
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 7),
    plot.caption       = element_text(size = 8, color = "grey50")
  )

ggsave("output/figures/bp/pi_dist_grid_bp_10yr.png", p_dist_grid_bp,
       width = 14, height = 6, dpi = 200)

saveRDS(
  list(P = bp10$P, pi0 = bp10$pi0, pistar = bp10$pistar, n = bp10$n),
  "data/derived/matrices_bp.rds"
)
cat("\nDone. 6-state robustness figures in output/figures/bp/.\n")
cat("Wrote data/derived/matrices_bp.rds\n")
