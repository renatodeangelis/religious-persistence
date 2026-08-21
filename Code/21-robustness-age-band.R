# ── 21 · AGE-BAND ROBUSTNESS (39–49) ─────────────────────────────────────────
# Restricts each birth cohort to respondents aged 39–49 at interview, then
# re-estimates transition matrices (6-state reltrad_bp) and computes π*, t_mix,
# and OM. Compares to full-sample estimates from matrices.rds.
#
# Input:  data/derived/gss_clean.rds, data/derived/matrices.rds
# Output: output/figures/age-band/

library(dplyr)
library(ggplot2)
library(ggpattern)
source("code/utils.R")

clean    = readRDS("data/derived/gss_clean.rds")
matrices = readRDS("data/derived/matrices.rds")

data      = clean$data
states_bp = clean$states_bp

cohorts_10 = c(1930, 1940, 1950, 1960, 1970, 1980)

cohort_labels = setNames(
  paste0(cohorts_10 - 5, "–", sprintf("%02d", (cohorts_10 + 4) %% 100)),
  as.character(cohorts_10)
)

dir.create("output/figures/age-band", recursive = TRUE, showWarnings = FALSE)

# ── CELL COUNT DIAGNOSTIC ─────────────────────────────────────────────────────

data_ab = data |>
  filter(
    age >= 39, age <= 49,
    !is.na(cohort_10), cohort_10 %in% cohorts_10,
    !is.na(reltrad16_bp), !is.na(reltrad_bp)
  )

cell_counts = data_ab |>
  group_by(cohort_10, reltrad16_bp) |>
  summarise(n = n(), .groups = "drop")

cat("\nCell counts (age 39–49, reltrad16_bp × cohort_10):\n")
print(as.data.frame(cell_counts))

thin_cells = filter(cell_counts, n < 30)
if (nrow(thin_cells) > 0) {
  cat("\nFlagged thin cells (N < 30):\n")
  print(thin_cells)
} else {
  cat("No thin cells detected.\n")
}

# ── ESTIMATE AGE-RESTRICTED MATRICES ─────────────────────────────────────────

P_ab      = list()
pi0_ab    = list()
pistar_ab = list()
n_ab      = list()

for (coh in cohorts_10) {
  sub = data_ab[data_ab$cohort_10 == coh, ]
  key = as.character(coh)
  if (nrow(sub) < 30) {
    cat(sprintf("Cohort %d: N = %d — skipped\n", coh, nrow(sub)))
    next
  }
  P_ab[[key]]      = p_matrix(sub, "reltrad16_bp", "reltrad_bp", levels = states_bp)
  pi0_ab[[key]]    = pi_0(sub, "reltrad16_bp")
  pistar_ab[[key]] = pi_star(P_ab[[key]])
  n_ab[[key]]      = nrow(sub)
}

# ── MIXING TIME ───────────────────────────────────────────────────────────────

t_mix_fn = function(P, pistar, max_t = 50, tol = 0.01) {
  for (t in seq_len(max_t)) {
    Pt    = P %^% t
    tv_mx = max(apply(Pt, 1, function(row) tv_norm(row, pistar)))
    if (tv_mx < tol) return(t)
  }
  cat(sprintf("  t_mix did not converge within %d steps\n", max_t))
  NA_integer_
}

# ── SUMMARY TABLE ─────────────────────────────────────────────────────────────

P_full      = matrices$nat10$P
pi0_full    = matrices$nat10$pi0
pistar_full = matrices$nat10$pistar
N_full      = matrices$nat10$N

summary_rows = lapply(cohorts_10, function(coh) {
  key = as.character(coh)

  om_full   = if (!is.null(P_full[[key]]))   overall_mobility(P_full[[key]],   pi0_full[[key]])   else NA_real_
  tmix_full = if (!is.null(P_full[[key]]))   t_mix_fn(P_full[[key]],   pistar_full[[key]])        else NA_integer_
  om_ab     = if (!is.null(P_ab[[key]]))     overall_mobility(P_ab[[key]],     pi0_ab[[key]])     else NA_real_
  tmix_ab   = if (!is.null(P_ab[[key]]))     t_mix_fn(P_ab[[key]],     pistar_ab[[key]])          else NA_integer_

  data.frame(
    cohort     = coh,
    cohort_lbl = cohort_labels[[key]],
    n_full     = if (!is.null(N_full[[key]])) sum(N_full[[key]]) else NA_integer_,
    n_ab       = if (!is.null(n_ab[[key]]))   n_ab[[key]]        else NA_integer_,
    om_full    = round(om_full,   4),
    om_ab      = round(om_ab,     4),
    tmix_full  = tmix_full,
    tmix_ab    = tmix_ab,
    stringsAsFactors = FALSE
  )
})

summary_df = do.call(rbind, summary_rows)

cat("\n── Summary: OM and t_mix by cohort ────────────────────────────────\n")
print(summary_df, row.names = FALSE)

# ── FIGURE A: π* PAIRED BAR ───────────────────────────────────────────────────

pistar_paired_df = do.call(rbind, lapply(cohorts_10, function(coh) {
  key      = as.character(coh)
  ps_full  = pistar_full[[key]]
  ps_ab    = pistar_ab[[key]]
  if (is.null(ps_full) || is.null(ps_ab)) return(NULL)
  lbl = cohort_labels[[key]]
  rbind(
    data.frame(cohort = lbl, religion = names(ps_full),
               value = as.numeric(ps_full), measure = "Full sample"),
    data.frame(cohort = lbl, religion = names(ps_ab),
               value = as.numeric(ps_ab),   measure = "Age 39–49")
  )
}))

pistar_paired_df$religion = factor(pistar_paired_df$religion, levels = rel_level_order)
pistar_paired_df$cohort   = factor(pistar_paired_df$cohort,
  levels = cohort_labels[as.character(cohorts_10)])
pistar_paired_df$measure  = factor(pistar_paired_df$measure,
  levels = c("Full sample", "Age 39–49"))

p_pistar = ggplot(pistar_paired_df, aes(x = religion, y = value)) +
  ggpattern::geom_col_pattern(
    aes(fill = religion, pattern = measure),
    position             = position_dodge(width = 0.8), width = 0.8,
    color                = "grey25", linewidth = 0.25,
    pattern_fill         = "white", pattern_colour = "white",
    pattern_angle        = 45,      pattern_density = 0.08,
    pattern_spacing      = 0.03,    pattern_key_scale_factor = 0.6
  ) +
  scale_fill_manual(values = reltrad_colors,
                    labels = reltrad_labels_tc, guide = "none") +
  ggpattern::scale_pattern_manual(
    values = c("Full sample" = "none", "Age 39–49" = "stripe"),
    labels = c(
      "Full sample"   = expression(Unadjusted~(pi[infinity])),
      "Age 39–49" = expression(Age~39-49~(pi[infinity]))
    )
  ) +
  guides(pattern = guide_legend(
    override.aes = list(fill = "grey55", color = "grey25")
  )) +
  scale_x_discrete(labels = reltrad_labels_tc) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ cohort, nrow = 2) +
  labs(x = NULL, y = "Share") +
  theme_bc(base_size = 12, x_angle = 45) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("output/figures/age-band/pistar_age_band.png", p_pistar,
       width = 10, height = 5, dpi = 200)
cat("Saved output/figures/age-band/pistar_age_band.png\n")

# ── FIGURE B: OM LINE PLOT ────────────────────────────────────────────────────

om_long = rbind(
  data.frame(cohort = summary_df$cohort, om = summary_df$om_full, sample = "Full sample"),
  data.frame(cohort = summary_df$cohort, om = summary_df$om_ab,   sample = "Age 39–49")
)
om_long = om_long[!is.na(om_long$om), ]
om_long$sample = factor(om_long$sample, levels = c("Full sample", "Age 39–49"))

p_om = ggplot(om_long, aes(x = cohort, y = om,
                            color = sample, shape = sample, group = sample)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Full sample" = "#0072B2", "Age 39–49" = "#D55E00"),
                     name = NULL) +
  scale_shape_manual(values = c("Full sample" = 16, "Age 39–49" = 17),
                     name = NULL) +
  scale_x_continuous(breaks = cohorts_10,
                     labels = cohort_labels[as.character(cohorts_10)]) +
  scale_y_continuous(limits = c(0.2, 0.4)) +
  labs(x = "Birth cohorts", y = "Probability to move") +
  theme_bc(base_size = 14, x_angle = 0) +
  theme(legend.position = c(0.5, 0.92),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color = NA))

ggsave("output/figures/age-band/om_age_band.png", p_om,
       width = 7, height = 5, dpi = 200)
cat("Saved output/figures/age-band/om_age_band.png\n")

# ── FIGURE C: t_mix LINE PLOT ─────────────────────────────────────────────────

tmix_long = rbind(
  data.frame(cohort = summary_df$cohort, tmix = summary_df$tmix_full, sample = "Full sample"),
  data.frame(cohort = summary_df$cohort, tmix = summary_df$tmix_ab,   sample = "Age 39–49")
)
tmix_long = tmix_long[!is.na(tmix_long$tmix), ]
tmix_long$sample = factor(tmix_long$sample, levels = c("Full sample", "Age 39–49"))

p_tmix = ggplot(tmix_long, aes(x = cohort, y = tmix,
                                color = sample, shape = sample, group = sample)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Full sample" = "#0072B2", "Age 39–49" = "#D55E00"),
                     name = NULL) +
  scale_shape_manual(values = c("Full sample" = 16, "Age 39–49" = 17),
                     name = NULL) +
  scale_x_continuous(breaks = cohorts_10,
                     labels = cohort_labels[as.character(cohorts_10)]) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), limits = c(0, 40)) +
  labs(x = "Birth cohorts",
       y = "Mixing Time (generations)") +
  theme_bc(base_size = 14, x_angle = 0) +
  theme(legend.position = c(0.5, 0.92),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color = NA))

ggsave("output/figures/age-band/tmix_age_band.png", p_tmix,
       width = 7, height = 5, dpi = 200)
cat("Saved output/figures/age-band/tmix_age_band.png\n")
