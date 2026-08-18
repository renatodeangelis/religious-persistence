# ── 06 · MOVEMENT MEASURES ────────────────────────────────────────────────────
# Implements mobility-as-movement measures from Wodtke et al. (2026, NBER w34800)
# for the observed one-step transition (t = 1) at each birth cohort window:
#
#   OM   = 1 − Σᵢ m_{i,0} · Pᵢᵢ           (overall mobility)
#   SM   = (1/2) Σᵢ |m_{i,0} − m_{i,1}|   (structural mobility)
#   EM   = OM − SM                          (exchange mobility)
#   SSM  = 1 − Σᵢ π*ᵢ · Pᵢᵢ               (steady-state mobility)
#   1−Pᵢᵢ per origin class                 (diagonal exit probability)
#
# overall_mobility(), sm(), and em() in utils.R implement these formulas directly;
# SSM reuses overall_mobility() with π* substituted for the initial distribution.
#
# Input:  data/derived/matrices.rds  (national 10-year matrices)
#         data/derived/gss_clean.rds  (5-year matrices built here)
# Output: output/figures/movement/*.png

source("code/utils.R")

matrices   = readRDS("data/derived/matrices.rds")
clean     = readRDS("data/derived/gss_clean.rds")
data      = clean$data
states_bp = clean$states_bp

# ── CONSTANTS ─────────────────────────────────────────────────────────────────

# Okabe-Ito palette — same keys and values as in 05-memory-measures.R
reltrad_colors = c(
  catholic           = "#0072B2",
  evangelical        = "#D55E00",
  `black protestant` = "#56B4E9",
  mainline           = "#009E73",
  other              = "#CC79A7",
  none               = "#999999"
)
reltrad_labels_tc = c(
  catholic           = "Catholic",
  evangelical        = "Evangelical",
  `black protestant` = "Black Protestant",
  mainline           = "Mainline",
  other              = "Other",
  none               = "None"
)

dir.create("output/figures/movement", recursive = TRUE, showWarnings = FALSE)

# ── HELPER: all four aggregate measures for one cohort ────────────────────────

movement_row = function(key, P, pi0, pistar) {
  data.frame(
    cohort = as.numeric(key),
    OM     = overall_mobility(P, pi0),
    SM     = sm(P, pi0, t = 0),
    EM     = em(P, pi0, t = 0),
    SSM    = overall_mobility(P, pistar),
    row.names = NULL
  )
}

# ── 10-YEAR COHORT COMPUTATION ────────────────────────────────────────────────

P_list_10      = matrices$nat10$P
pi0_list_10    = matrices$nat10$pi0
pistar_list_10 = matrices$nat10$pistar

# Aggregate movement measures across 10-year windows
mov_df_10 = do.call(rbind, lapply(names(P_list_10), function(key) {
  movement_row(key, P_list_10[[key]], pi0_list_10[[key]], pistar_list_10[[key]])
}))
mov_df_10 = mov_df_10[order(mov_df_10$cohort), ]

# Cohort axis labels: "1925–34", "1935–44", …  (follows 05 convention)
mov_df_10$label = paste0(mov_df_10$cohort - 5, "–",
                         sprintf("%02d", (mov_df_10$cohort + 4) %% 100))
mov_df_10$label = factor(mov_df_10$label, levels = mov_df_10$label)

# Diagonal exit probabilities (1 − Pᵢᵢ) per origin class
diag_df_10 = do.call(rbind, lapply(names(P_list_10), function(key) {
  d = 1 - diag(as.matrix(P_list_10[[key]]))
  data.frame(cohort = as.numeric(key), origin = names(d), exit_prob = d, row.names = NULL)
}))
diag_df_10 = diag_df_10[order(diag_df_10$cohort), ]
diag_df_10$origin = factor(diag_df_10$origin, levels = rel_level_order)
diag_df_10$label  = paste0(diag_df_10$cohort - 5, "–",
                            sprintf("%02d", (diag_df_10$cohort + 4) %% 100))
diag_df_10$label  = factor(diag_df_10$label, levels = unique(diag_df_10$label))

# ── 5-YEAR COHORT MATRICES (built here; not stored in matrices.rds) ───────────

# Bin start years: 1925, 1930, …, 1980 (12 five-year windows)
data$cohort_5 = floor((data$cohort - 1925) / 5) * 5 + 1925
cohorts_5 = sort(unique(data$cohort_5[!is.na(data$cohort_5) &
                                        data$cohort_5 >= 1925 &
                                        data$cohort_5 <= 1980]))

P_list_5      = list()
pi0_list_5    = list()
pistar_list_5 = list()

for (coh in cohorts_5) {
  sub = data[!is.na(data$cohort_5) & data$cohort_5 == coh &
               !is.na(data$reltrad16_bp) & !is.na(data$reltrad_bp), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)
  P_list_5[[key]]      = suppressWarnings(
    p_matrix(sub, "reltrad16_bp", "reltrad_bp", levels = states_bp)
  )
  pi0_list_5[[key]]    = pi_0(sub, "reltrad16_bp")
  pistar_list_5[[key]] = pi_star(P_list_5[[key]])
}

# Aggregate movement measures across 5-year windows
mov_df_5 = do.call(rbind, lapply(names(P_list_5), function(key) {
  movement_row(key, P_list_5[[key]], pi0_list_5[[key]], pistar_list_5[[key]])
}))
mov_df_5 = mov_df_5[order(mov_df_5$cohort), ]

# Cohort axis labels: "1925–29", "1930–34", …
mov_df_5$label = paste0(mov_df_5$cohort, "–",
                        sprintf("%02d", (mov_df_5$cohort + 4) %% 100))
mov_df_5$label = factor(mov_df_5$label, levels = mov_df_5$label)

# Diagonal exit probabilities (5yr)
diag_df_5 = do.call(rbind, lapply(names(P_list_5), function(key) {
  d = 1 - diag(as.matrix(P_list_5[[key]]))
  data.frame(cohort = as.numeric(key), origin = names(d), exit_prob = d, row.names = NULL)
}))
diag_df_5 = diag_df_5[order(diag_df_5$cohort), ]
diag_df_5$origin = factor(diag_df_5$origin, levels = rel_level_order)
diag_df_5$label  = paste0(diag_df_5$cohort, "–",
                           sprintf("%02d", (diag_df_5$cohort + 4) %% 100))
diag_df_5$label  = factor(diag_df_5$label, levels = unique(diag_df_5$label))

# ── FIGURE 1: OM TREND ────────────────────────────────────────────────────────

make_om_plot = function(df, title_str) {
  ggplot(df, aes(x = label, y = OM, group = 1)) +
    geom_line(linewidth = 0.9, color = "#0072B2") +
    geom_point(size = 2.5, color = "#0072B2") +
    scale_y_continuous(limits = c(0, 0.5)) +
    labs(x = NULL, y = "Overall mobility (OM)", title = title_str) +
    healy_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggsave("output/figures/movement/om_trend_10yr.png",
       make_om_plot(mov_df_10, "Overall Mobility by Birth Cohort (10-year bins)"),
       width = 7, height = 5, dpi = 200)

ggsave("output/figures/movement/om_trend_5yr.png",
       make_om_plot(mov_df_5, "Overall Mobility by Birth Cohort (5-year bins)"),
       width = 9, height = 5, dpi = 200)

# ── FIGURE 2: SM AND EM COMBINED ──────────────────────────────────────────────

make_sm_em_plot = function(df, title_str) {
  long = rbind(
    data.frame(cohort = df$cohort, label = df$label,
               measure = "Structural (SM)", value = df$SM),
    data.frame(cohort = df$cohort, label = df$label,
               measure = "Exchange (EM)",   value = df$EM)
  )
  long$measure = factor(long$measure, levels = c("Structural (SM)", "Exchange (EM)"))
  long$label   = factor(long$label, levels = levels(df$label))

  ggplot(long, aes(x = label, y = value, color = measure, group = measure)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("Structural (SM)" = "#D55E00",
                                  "Exchange (EM)"   = "#009E73"),
                       name = NULL) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(x = NULL, y = "Mobility", title = title_str) +
    healy_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggsave("output/figures/movement/sm_em_trend_10yr.png",
       make_sm_em_plot(mov_df_10,
                       "Structural and Exchange Mobility by Birth Cohort (10-year bins)"),
       width = 7, height = 5, dpi = 200)

ggsave("output/figures/movement/sm_em_trend_5yr.png",
       make_sm_em_plot(mov_df_5,
                       "Structural and Exchange Mobility by Birth Cohort (5-year bins)"),
       width = 9, height = 5, dpi = 200)

# ── FIGURE 3: SSM TREND ───────────────────────────────────────────────────────

make_ssm_plot = function(df, title_str) {
  ggplot(df, aes(x = label, y = SSM, group = 1)) +
    geom_line(linewidth = 0.9, color = "#CC79A7") +
    geom_point(size = 2.5, color = "#CC79A7") +
    scale_y_continuous(limits = c(0, 0.5)) +
    labs(x = NULL, y = "Steady-state mobility (SSM)", title = title_str) +
    healy_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggsave("output/figures/movement/ssm_trend_10yr.png",
       make_ssm_plot(mov_df_10, "Steady-State Mobility by Birth Cohort (10-year bins)"),
       width = 7, height = 5, dpi = 200)

ggsave("output/figures/movement/ssm_trend_5yr.png",
       make_ssm_plot(mov_df_5, "Steady-State Mobility by Birth Cohort (5-year bins)"),
       width = 9, height = 5, dpi = 200)

# ── FIGURE 4: DIAGONAL EXIT PROBABILITIES ─────────────────────────────────────

make_diag_exit_plot = function(df, title_str) {
  ggplot(df, aes(x = label, y = exit_prob, color = origin, group = origin)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc, name = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Exit probability (1 − Pᵢᵢ)", title = title_str) +
    healy_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggsave("output/figures/movement/diag_exit_10yr.png",
       make_diag_exit_plot(diag_df_10,
                           "Diagonal Exit Probability by Origin Class (10-year bins)"),
       width = 7, height = 5, dpi = 200)

ggsave("output/figures/movement/diag_exit_5yr.png",
       make_diag_exit_plot(diag_df_5,
                           "Diagonal Exit Probability by Origin Class (5-year bins)"),
       width = 9, height = 5, dpi = 200)

# ── CONSOLE SUMMARY ───────────────────────────────────────────────────────────

cat("\nMovement measures (10-year cohorts):\n")
print(mov_df_10[, c("label", "OM", "SM", "EM", "SSM")], row.names = FALSE, digits = 3)
cat("\nWrote 8 figures to output/figures/movement/\n")
