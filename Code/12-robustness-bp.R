# ── 12 · ROBUSTNESS: 6-STATE SPACE (BLACK PROTESTANT SEPARATE) ───────────────
# Re-estimates cohort transition matrices on the 6-state space that keeps Black
# Protestant as its own category (only jewish → other collapsed), at both 20-
# and 10-year resolution. Reports transition heatmaps, IM memory curves, and
# diagonal persistence. Uses the reltrad_bp/reltrad16_bp columns built in 01.
#
# Input:  data/derived/gss_clean.rds
# Output: output/figures/bp/*.png

library(dplyr)
library(ggplot2)
source("code/utils.R")

clean     = readRDS("data/derived/gss_clean.rds")
data      = clean$data
states_bp = clean$states_bp

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
    do.call(rbind, lapply(0:6, function(t) {
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

# ── 10-YEAR COHORTS (edges 1940–1980) ────────────────────────────────────────

mids_10 = c(1945, 1955, 1965, 1975, 1985)   # edges 1940–1980
cell_diag("cohort_10", mids_10, 5, "Cell-count diagnostic (6-state, 10-year cohorts, 1940–1989)")
bp10 = build_bp("cohort_10", mids_10, 5)

for (key in names(bp10$P)) {
  cat("\n── Cohort", key, "–", as.integer(key) + 9, "  (N =", bp10$n[[key]], ") ──\n")
  print(round(bp10$P[[key]], 3))
}

for (key in names(bp10$P)) {
  p = make_combined(
    bp10$P[[key]], bp10$pi0[[key]], bp10$pistar[[key]],
    levels    = rel_level_order_6,
    title_str = paste0("6-State — Cohort ", key, "–", as.integer(key) + 9,
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
  facet_wrap(~ cohort, nrow = 1,
             labeller = labeller(cohort = function(x) paste0(x, "–", as.integer(x) + 9))) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)", color = NULL,
       title = "Individual Memory by Cohort — 6-State Space (10-year bins, 1940–1989)") +
  healy_theme

ggsave("output/figures/bp/im_memory_10yr_6state.png", p_im_bp10, width = 14, height = 5, dpi = 200)

diag_df10 = diag_data(bp10)

p_diag10 = ggplot(diag_df10, aes(x = cohort, y = persistence, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = reltrad_colors_6, labels = reltrad_labels_6) +
  scale_x_continuous(breaks = c(1940, 1950, 1960, 1970, 1980)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Birth cohort (10-year bins)", y = "Diagonal persistence P[i → i]", color = NULL,
       title = "Diagonal Persistence by Origin — 6-State Space (10-year cohorts, 1940–1989)") +
  healy_theme

ggsave("output/figures/bp/diagonal_persistence_10yr_6state.png",
       p_diag10, width = 8, height = 5, dpi = 200)

cat("\nDone. 6-state robustness figures in output/figures/bp/.\n")
