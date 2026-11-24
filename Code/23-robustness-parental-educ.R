# ── 23 · ROBUSTNESS: PARENTAL-EDUCATION-STRATIFIED MATRICES ───────────────────
# 20-year cohort-window transition matrices (6-state scheme), estimated
# separately for respondents whose highest-educated parent has a college degree
# (>= 16 years) vs. not. Parental education = higher of paeduc and maeduc.
#
# Input:  data/derived/gss_clean.rds  (needs paeduc, maeduc — see 01)
# Output: data/derived/matrices_parented.rds
#         output/figures/parental-educ/*.png

source("code/utils.R")
library(tidyr)

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_bp  = clean$states_bp

if (!all(c("paeduc", "maeduc") %in% names(data)))
  stop("paeduc/maeduc not in gss_clean.rds: rerun 01-prepare-data.R.")

# 6-state display order (Black Protestant separate), from utils.R
rel_order_bp = rel_level_order

# 20-year window midpoints (edges 1925–44, 1945–64, 1965–84)
mids20     = c(1935, 1955, 1975)
cohort_lbl = function(mid) paste0(mid - 10, "–", sprintf("%02d", (mid + 9) %% 100))
grp_levels = c("no_college", "college")
grp_labels = c(no_college = "No college", college = "College")
# Line-plot styling (same Okabe-Ito colors/shapes as the age-band figures)
grp_colors = c(no_college = "#0072B2", college = "#D55E00")
grp_shapes = c(no_college = 16,        college = 17)

# Mixing-time tolerance (convergence criterion from the paper's footnote)
mix_eps = 0.01

# ── MIXING TIME ───────────────────────────────────────────────────────────────
# t_mix = min{ t : d'(t) <= ε }, where
#   d'(t) = 0.5 * max_{i,k} Σ_j |(P^t)_ij − (P^t)_kj|
# is the maximum total-variation distance between any two rows of P^t. Once it is
# below ε, every row is (numerically) equal to π∞: destination is independent of
# origin. Returns NA if not reached within t_max generations.
mixing_time = function(P, eps = 0.01, t_max = 1000) {
  P     = as.matrix(P)
  Pt    = diag(nrow(P))
  pairs = combn(nrow(P), 2)                  # all origin-row pairs (i, k)
  for (t in seq_len(t_max)) {
    Pt = Pt %*% P
    d  = max(apply(pairs, 2, function(ix) 0.5 * sum(abs(Pt[ix[1], ] - Pt[ix[2], ]))))
    if (d <= eps) return(t)
  }
  NA_real_
}

# ── PARENTAL EDUCATION RECODE ─────────────────────────────────────────────────
# pared = max(paeduc, maeduc); if only one parent is reported, that parent's
# value is used. NA only when both are missing.

data = data |>
  mutate(
    paeduc = as.numeric(paeduc),
    maeduc = as.numeric(maeduc),
    pared  = if_else(is.na(paeduc) & is.na(maeduc), NA_real_,
                     pmax(paeduc, maeduc, na.rm = TRUE)),
    pared_source = case_when(
      is.na(paeduc) & is.na(maeduc) ~ "neither",
      is.na(paeduc) | is.na(maeduc) ~ "one parent",
      TRUE                          ~ "both parents"
    ),
    pared_grp = case_when(
      pared >= 16 ~ "college",
      pared <  16 ~ "no_college"
    ),
    cohort_20 = (floor((cohort - 1925) / 20) * 20 + 1925) + 10
  )

# ── MISSINGNESS REPORT ────────────────────────────────────────────────────────
# "one parent" cases can understate pared (the unreported parent may hold a degree).

miss_tab = data |>
  filter(!is.na(reltrad16_bp), !is.na(reltrad_bp), cohort_20 %in% mids20) |>
  group_by(cohort = sapply(cohort_20, cohort_lbl)) |>
  summarise(n_total      = n(),
            n_neither    = sum(pared_source == "neither"),
            n_one_parent = sum(pared_source == "one parent"),
            pct_dropped  = round(100 * mean(pared_source == "neither"), 1),
            pct_one_par  = round(100 * mean(pared_source == "one parent"), 1),
            .groups = "drop")
cat("\n── Parental-education coverage by cohort window ──\n")
print(as.data.frame(miss_tab), row.names = FALSE)

# ── BUILD MATRICES ────────────────────────────────────────────────────────────

P_pe = pi0_pe = pistar_pe = n_pe = nrow_pe = list()

for (g in grp_levels) {
  for (mid in mids20) {
    sub = data[!is.na(data$pared_grp) & data$pared_grp == g &
               !is.na(data$cohort_20) & data$cohort_20 == mid &
               !is.na(data$reltrad16_bp) & !is.na(data$reltrad_bp), ]
    if (nrow(sub) < 30) next
    key = paste(g, mid, sep = "_")
    P_pe[[key]]      = p_matrix(sub, "reltrad16_bp", "reltrad_bp", levels = states_bp)
    pi0_pe[[key]]    = pi_0(sub, "reltrad16_bp")
    pistar_pe[[key]] = pi_star(P_pe[[key]])
    n_pe[[key]]      = nrow(sub)
    nrow_pe[[key]]   = rowSums(count_matrix(sub, "reltrad16_bp", "reltrad_bp",
                                            levels = states_bp))
  }
}

# ── SAMPLE-SIZE TABLE (group × cohort × origin state); flag N < 30 ────────────

size_tab = do.call(rbind, lapply(names(nrow_pe), function(key) {
  g   = sub("_\\d{4}$", "", key)
  mid = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  data.frame(group = grp_labels[g], cohort = cohort_lbl(mid),
             origin = names(nrow_pe[[key]]), n = as.integer(nrow_pe[[key]]),
             thin = nrow_pe[[key]] < 30, row.names = NULL)
}))
cat("\n── Origin-row N (group × cohort) ──\n")
print(size_tab |> pivot_wider(id_cols = c(group, cohort), names_from = origin,
                              values_from = n) |> as.data.frame(), row.names = FALSE)
if (any(size_tab$thin)) {
  cat("\nTHIN ORIGIN ROWS (N < 30):\n")
  print(size_tab[size_tab$thin, c("group", "cohort", "origin", "n")], row.names = FALSE)
}

# Console: matrices
for (key in names(P_pe)) {
  g   = sub("_\\d{4}$", "", key)
  mid = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  cat("\n──", grp_labels[g], "| Cohort", cohort_lbl(mid), "(N =", n_pe[[key]], ") ──\n")
  print(round(P_pe[[key]], 3))
}

dir.create("output/figures/parental-educ", recursive = TRUE, showWarnings = FALSE)

# ── FIGURE 1: π₀ vs π∞ BAR CHART ──────────────────────────────────────────────
# Same layout as pi_dist_20yr_nativity (13): bars colored by religion, striped = π₀,
# solid = π∞; rows = parental-education group, columns = cohort window.

pi_df = do.call(rbind, lapply(names(P_pe), function(key) {
  g   = sub("_\\d{4}$", "", key)
  mid = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  rbind(
    data.frame(cohort = cohort_lbl(mid), group = g, measure = "π₀",
               religion = names(pi0_pe[[key]]),    value = as.numeric(pi0_pe[[key]])),
    data.frame(cohort = cohort_lbl(mid), group = g, measure = "π∞",
               religion = names(pistar_pe[[key]]), value = as.numeric(pistar_pe[[key]]))
  )
}))
pi_df$religion = factor(pi_df$religion, levels = rel_order_bp)
pi_df$cohort   = factor(pi_df$cohort,   levels = sapply(mids20, cohort_lbl))
pi_df$measure  = factor(pi_df$measure,  levels = c("π₀", "π∞"))
pi_df$group    = factor(pi_df$group,    levels = grp_levels, labels = grp_labels[grp_levels])

p_pi = ggplot(pi_df, aes(x = religion, y = value)) +
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
  facet_grid(group ~ cohort) +
  labs(x = NULL, y = "Share") +
  theme_bc(base_size = 12, x_angle = 45) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("output/figures/parental-educ/pi_bar_parented.png", p_pi,
       width = 10, height = 6, dpi = 200)

# ── FIGURE 2: OVERALL MOBILITY ────────────────────────────────────────────────
# overall_mobility(P, π₀) = 1 − Σ π₀ᵢ Pᵢᵢ  (π₀-weighted probability of leaving origin)

om_df = do.call(rbind, lapply(names(P_pe), function(key) {
  g   = sub("_\\d{4}$", "", key)
  mid = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  data.frame(mid = mid, group = g,
             om = overall_mobility(P_pe[[key]], pi0_pe[[key]]))
}))
om_df$group = factor(om_df$group, levels = grp_levels)

# Styled like om_age_band (21): publication theme, legend inside at top.
p_om = ggplot(om_df, aes(x = mid, y = om, color = group, shape = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = grp_colors, labels = grp_labels, name = NULL) +
  scale_shape_manual(values = grp_shapes, labels = grp_labels, name = NULL) +
  scale_x_continuous(breaks = mids20, labels = sapply(mids20, cohort_lbl)) +
  scale_y_continuous(limits = c(0.2, 0.4)) +
  labs(x = "Birth cohorts", y = "Probability to move") +
  theme_bc(base_size = 14, x_angle = 0) +
  theme(legend.position = c(0.5, 0.92),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color = NA))

ggsave("output/figures/parental-educ/overall_mobility_parented.png", p_om,
       width = 7, height = 5, dpi = 200)

# ── FIGURE 3: MIXING TIME ─────────────────────────────────────────────────────
# One mixing time per group × cohort window (generations until all rows of P^t
# are within ε of each other in total variation).

mix_df = do.call(rbind, lapply(names(P_pe), function(key) {
  g   = sub("_\\d{4}$", "", key)
  mid = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  data.frame(mid = mid, group = g, t_mix = mixing_time(P_pe[[key]], eps = mix_eps))
}))
mix_df$group = factor(mix_df$group, levels = grp_levels)

cat("\n── Mixing time (max pairwise row TV <=", mix_eps, ") ──\n")
print(mix_df |> pivot_wider(id_cols = group, names_from = mid,
                            values_from = t_mix) |> as.data.frame(), row.names = FALSE)
if (anyNA(mix_df$t_mix))
  cat("NOTE: some matrices did not mix within 1000 generations (plotted as gaps).\n")

# Styled like tmix_age_band (21): publication theme, legend inside at top.
p_mix = ggplot(mix_df, aes(x = mid, y = t_mix, color = group, shape = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = grp_colors, labels = grp_labels, name = NULL) +
  scale_shape_manual(values = grp_shapes, labels = grp_labels, name = NULL) +
  scale_x_continuous(breaks = mids20, labels = sapply(mids20, cohort_lbl)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(0, 40)) +
  labs(x = "Birth cohorts", y = "Mixing Time (generations)") +
  theme_bc(base_size = 14, x_angle = 0) +
  theme(legend.position = c(0.5, 0.92),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color = NA))

ggsave("output/figures/parental-educ/mixing_time_parented.png", p_mix,
       width = 7, height = 5, dpi = 200)

# ── SAVE ESTIMATES ────────────────────────────────────────────────────────────
# .rds: everything, in list form. CSVs: long-format tables for easy sharing.

saveRDS(list(P = P_pe, pi0 = pi0_pe, pistar = pistar_pe, n = n_pe,
             n_row = nrow_pe, coverage = miss_tab, mixing = mix_df,
             overall_mobility = om_df),
        "data/derived/matrices_parented.rds")

dir.create("output/tables/parental-educ", recursive = TRUE, showWarnings = FALSE)

# key → group and cohort-window labels
key_group  = function(key) grp_labels[sub("_\\d{4}$", "", key)]
key_cohort = function(key) cohort_lbl(as.integer(sub(".*_(\\d{4})$", "\\1", key)))

# 1. Transition matrices: one row per group × cohort × origin × destination
csv_P = do.call(rbind, lapply(names(P_pe), function(key) {
  df = as.data.frame(as.table(P_pe[[key]]), stringsAsFactors = FALSE)
  names(df) = c("origin", "current", "p")
  df$n_origin = as.integer(nrow_pe[[key]][df$origin])
  data.frame(group = unname(key_group(key)), cohort = key_cohort(key), df)
}))
write.csv(csv_P, "output/tables/parental-educ/transition_matrices.csv", row.names = FALSE)

# 2. Origin (π₀) and stationary (π∞) distributions
csv_pi = do.call(rbind, lapply(names(P_pe), function(key) {
  data.frame(group = unname(key_group(key)), cohort = key_cohort(key),
             religion = names(pi0_pe[[key]]),
             pi0 = as.numeric(pi0_pe[[key]]),
             pistar = as.numeric(pistar_pe[[key]][names(pi0_pe[[key]])]))
}))
write.csv(csv_pi, "output/tables/parental-educ/pi_distributions.csv", row.names = FALSE)

# 3. Summary measures: N, overall mobility, mixing time per group × cohort
csv_sum = data.frame(
  group  = unname(sapply(names(P_pe), key_group)),
  cohort = sapply(names(P_pe), key_cohort),
  n      = unlist(n_pe[names(P_pe)]),
  overall_mobility = sapply(names(P_pe), function(k) overall_mobility(P_pe[[k]], pi0_pe[[k]])),
  mixing_time      = sapply(names(P_pe), function(k) mixing_time(P_pe[[k]], eps = mix_eps)),
  row.names = NULL
)
write.csv(csv_sum, "output/tables/parental-educ/summary_measures.csv", row.names = FALSE)

cat("Wrote output/figures/parental-educ/\n")
cat("Wrote output/tables/parental-educ/*.csv\n")
cat("Wrote data/derived/matrices_parented.rds\n")
