library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gssr)

reltrad_labels = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

data(gss_all)
data = gss_all |>
  select(year, cohort, reltrad, reltrad16) |>
  filter(!(year %in% c(1972, 2021))) |>
  mutate(across(c(reltrad, reltrad16),
                ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  mutate(across(c(reltrad, reltrad16),
                ~ case_when(
                  . == "jewish"           ~ "other",
                  . == "black protestant" ~ "evangelical",
                  TRUE                    ~ .
                ),
                .names = "{.col}_alt")) |>
  mutate(age = year - cohort) |>
  filter(age > 25, cohort >= 1900) |>
  mutate(
    cohort_5  = floor((cohort - 1900) / 5)  * 5  + 1900,
    cohort_10 = floor((cohort - 1900) / 10) * 10 + 1900
  )

# ── FUNCTIONS ────────────────────────────────────────────────────────────────────

# Matrix power operator (replaces expm::`%^%`)
`%^%` = function(M, k) {
  if (k == 0) {
    I = diag(nrow(M))
    dimnames(I) = dimnames(M)
    return(I)
  }
  result = M
  for (i in seq_len(k - 1)) result = result %*% M
  result
}

tv_norm = function(mu, nu) {
  0.5 * sum(abs(mu - nu))
}

# Initial distribution (π₀): unweighted share of each origin state.
pi_0 = function(data, origin) {
  tab = table(data[[origin]])
  v = as.numeric(tab / sum(tab))
  names(v) = names(tab)
  v
}

# Stationary distribution (π*): left eigenvector of P corresponding to eigenvalue 1.
pi_star = function(P) {
  eig = eigen(t(as.matrix(P)))
  v = Re(eig$vectors[, which.min(abs(eig$values - 1))])
  if (any(v < 0)) v = abs(v)
  setNames(v / sum(v), rownames(P))
}

# Unweighted row-stochastic transition matrix.
# levels: optional character vector to fix state ordering across cohort subsets.
p_matrix = function(data, origin, current, levels = NULL) {
  f   = if (!is.null(levels)) function(x) factor(x, levels = levels) else factor
  tab = table(f(data[[origin]]), f(data[[current]]))

  # Warn on true zero cells (sampling or structural)
  zero_idx = which(tab == 0, arr.ind = TRUE)
  if (nrow(zero_idx) > 0) {
    zero_labels = paste(rownames(tab)[zero_idx[, 1]], "->",
                        colnames(tab)[zero_idx[, 2]], collapse = "; ")
    warning("Zero cells in P: ", zero_labels)
  }

  P = tab / rowSums(tab)
  class(P) = "matrix"
  P
}

# Individual memory (IM): log TV distance from π* for each origin state at step t.
im = function(data, origin, current, t = 1) {
  P_mat = p_matrix(data, origin, current)
  pi_s  = pi_star(P_mat)
  P_t   = P_mat %^% t
  im_i  = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  log(im_i)
}

# IM from a pre-built matrix (reuses P_list_* to avoid recomputing P).
im_from_P = function(P, t = 1) {
  pi_s = pi_star(P)
  P_t  = P %^% t
  im_i = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  setNames(log(im_i), rownames(P))
}

# ── TRANSITION MATRICES ──────────────────────────────────────────────────────

states_alt = sort(unique(c(data$reltrad_alt, data$reltrad16_alt)))
states_alt = states_alt[!is.na(states_alt)]

# ── 5-year cohort loop ───────────────────────────────────────────────────────
cohorts_5 = sort(unique(data$cohort_5[!is.na(data$cohort_5) & data$cohort_5 >= 1920 & data$cohort_5 <= 1980]))

P_list_5      = list()
pi0_list_5    = list()
pistar_list_5 = list()

for (coh in cohorts_5) {
  sub = data[!is.na(data$cohort_5) & data$cohort_5 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_5[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_5[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_5[[key]] = pi_star(P_list_5[[key]])
}

# ── 10-year cohort loop ──────────────────────────────────────────────────────
cohorts_10 = sort(unique(data$cohort_10[!is.na(data$cohort_10) & data$cohort_10 >= 1920 & data$cohort_10 <= 1980]))

P_list_10      = list()
pi0_list_10    = list()
pistar_list_10 = list()

for (coh in cohorts_10) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_10[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_10[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_10[[key]] = pi_star(P_list_10[[key]])
}

# ── IM LOOP (10-year cohorts, t = 0:4) ──────────────────────────────────────
im_rows_10 = vector("list", length(P_list_10))
names(im_rows_10) = names(P_list_10)

for (key in names(P_list_10)) {
  rows = lapply(0:4, function(t) {
    vals = im_from_P(P_list_10[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_10[[key]] = do.call(rbind, rows)
}

im_df_10 = do.call(rbind, im_rows_10)

# ── FIGURES ───────────────────────────────────────────────────────────────────

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")

plot_pmat_heatmap = function(P, levels = NULL, text_size = 5, title_str = "P") {
  df = as.data.frame(as.table(P))
  names(df) = c("origin", "current", "est")
  if (!is.null(levels)) {
    df$origin  = factor(df$origin,  levels = levels)
    df$current = factor(df$current, levels = rev(levels))
  }
  ggplot(df, aes(x = current, y = origin, fill = est)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", est)), size = text_size) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(x = "Current religion", y = "Origin religion",
         fill = "Transition Prob.", title = title_str) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
          axis.text.y     = element_text(angle = 45, hjust = 1, size = 13),
          axis.ticks.x    = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.title      = element_text(size = 13),
          legend.position = "bottom",
          legend.text     = element_text(size = 11),
          legend.title    = element_text(size = 12),
          plot.title      = element_text(hjust = 0.5, size = 16),
          panel.grid      = element_blank())
}

plot_pi_column = function(vec, title_str, levels = NULL) {
  df = data.frame(origin = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df$origin = factor(df$origin, levels = levels)
  } else {
    df$origin = factor(df$origin, levels = rev(unique(df$origin)))
  }
  ggplot(df, aes(x = 1, y = origin, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5, color = "black") +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(title = title_str) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none",
          plot.title   = element_text(hjust = 0.5, size = 16),
          panel.grid   = element_blank())
}

make_combined = function(P, pi0, pistar, levels = NULL, title_str = "P") {
  g      = plot_pmat_heatmap(P,      levels = levels, title_str = title_str)
  g0     = plot_pi_column(pi0,    title_str = "π₀", levels = levels)
  g_star = plot_pi_column(pistar, title_str = "π*",      levels = levels)
  patchwork::wrap_plots(g, g0, g_star, widths = c(6, 1, 1))
}

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_10)) {
  p = make_combined(P_list_10[[key]], pi0_list_10[[key]], pistar_list_10[[key]],
                    levels = rel_level_order, title_str = paste("Cohort", key))
  ggsave(paste0("output/figures/trans_", key, "_10yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

for (key in names(P_list_5)) {
  p = make_combined(P_list_5[[key]], pi0_list_5[[key]], pistar_list_5[[key]],
                    levels = rel_level_order, title_str = paste("Cohort", key))
  ggsave(paste0("output/figures/trans_", key, "_5yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

# ── IM FIGURE ────────────────────────────────────────────────────────────────
im_df_10$origin = factor(im_df_10$origin, levels = rel_level_order)

p_im = ggplot(im_df_10, aes(x = t, y = im, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = "Origin", title = "Individual Memory by Cohort (10-year bins, t = 0–4)") +
  theme_minimal() +
  theme(
    strip.text      = element_text(size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

ggsave("output/figures/im_memory_10yr.png", p_im, width = 12, height = 7, dpi = 200)

