# ── 16 · ROBUSTNESS: CHILDHOOD ATTENDANCE-STRATIFIED MATRICES ────────────────
# Transition matrices stratified by attend12 (religious service attendance
# frequency at age 12). attend12 is available only in GSS waves 1991, 1998,
# 2008, and 2018, so the sample is restricted to those four waves and cohorts
# are pooled across them to preserve cell counts.
#
# Binary recode — nearly weekly or more vs. less than nearly weekly:
#   regular     — nearly every week / every week / several times a week  (values 6–8)
#   non_regular — never through 2-3 times a month                        (values 0–5)
# (Value 10 = not applicable / NA, treated as missing.)
#
# "Weekly" is the canonical threshold in the sociology of religion literature
# and maps directly onto the strict-church/boundary-maintenance prediction
# (Iannaccone 1994): high-demand households should show slower memory decay.
# A three-category version would require ~3× the per-cell N; not feasible
# given the four-wave restriction.
#
# Run manually — not sourced by 00-run-all.R.
#
# Input:  gss_all (gssr package)          — rebuilt here; not gss_clean.rds
# Output: console (cell counts, λ₂, diagonal persistence)
#         output/figures/attend12/*.png

library(dplyr)
library(ggplot2)
library(patchwork)
library(gssr)
source("code/utils.R")

reltrad_labels  = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")

reltrad_colors = c(
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

attend12_waves = c(1991, 1998, 2008, 2018)

# ── DATA ──────────────────────────────────────────────────────────────────────
# Applies the same sample restrictions as 01-prepare-data.R (age 30–75,
# cohort 1925–1994, same excluded years) but loads attend12 directly from
# gss_all because it is not included in gss_clean.rds.

data(gss_all)

data = gss_all |>
  filter(year %in% attend12_waves) |>
  select(year, cohort, reltrad, reltrad16, attend12) |>
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
  mutate(cohort = as.numeric(cohort),
         age    = year - cohort) |>
  filter(age >= 30, age <= 75, cohort >= 1925, cohort <= 1994) |>
  haven::zap_labels()

states_alt = sort(unique(c(data$reltrad_alt, data$reltrad16_alt)))
states_alt = states_alt[!is.na(states_alt)]

# ── ATTEND12 RECODE ───────────────────────────────────────────────────────────
# GSS numeric scale: 0=never, 1=<1/yr, 2=1-2/yr, 3=several/yr,
#                    4=once/mo, 5=2-3/mo, 6=nearly every wk,
#                    7=every wk, 8=several/wk; 10=NA.

data = data |>
  mutate(
    attend12_num = as.numeric(attend12),
    attend12_num = if_else(attend12_num == 10, NA_real_, attend12_num),
    attend12_bin = case_when(
      attend12_num %in% 6:8 ~ "regular",
      attend12_num %in% 0:5 ~ "non_regular",
      TRUE                  ~ NA_character_
    ),
    attend12_bin = factor(attend12_bin, levels = c("non_regular", "regular"))
  )

# ── CELL COUNT REPORT ─────────────────────────────────────────────────────────

cat("\n── Overall sample N by wave ──\n")
print(table(data$year))

cat("\n── Raw attend12 distribution (full scale, all waves) ──\n")
print(table(data$attend12_num, useNA = "ifany"))

cat("\n── Binary recode distribution ──\n")
print(table(data$attend12_bin, useNA = "ifany"))

cat("\n── Binary recode by wave (check for differential missingness) ──\n")
print(table(data$year, data$attend12_bin, useNA = "ifany"))

cat("\n── Origin × binary recode (cell count check before matrix estimation) ──\n")
print(table(data$reltrad16_alt, data$attend12_bin, useNA = "ifany"))

# ── BUILD MATRICES ────────────────────────────────────────────────────────────
# Pool across the four available waves. Cohort-stratified matrices are not
# estimated here because the four-wave restriction would produce very thin cells.

att_levels = c("non_regular", "regular")
att_labels = c(non_regular = "Non-regular (< nearly weekly)",
               regular     = "Regular (nearly weekly or more)")

P_att      = list()
pi0_att    = list()
pistar_att = list()
n_att      = list()

for (att in att_levels) {
  sub = data[!is.na(data$attend12_bin) & data$attend12_bin == att &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) {
    cat("Skipping stratum '", att, "' — N =", nrow(sub), "(too small)\n")
    next
  }
  P_att[[att]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_att[[att]]    = pi_0(sub, "reltrad16_alt")
  pistar_att[[att]] = pi_star(P_att[[att]])
  n_att[[att]]      = nrow(sub)
}

# ── MEMORY MEASURES (CONSOLE) ─────────────────────────────────────────────────

cat("\n── Transition matrices by attend12 stratum ──\n")
for (att in names(P_att)) {
  cat("\n── ", att_labels[att], " | N =", n_att[[att]], "──\n")
  print(round(P_att[[att]], 3))
}

# λ₂: second-largest eigenvalue by magnitude (the memory summary)
lambda2 = function(P) {
  eigs = sort(abs(Re(eigen(as.matrix(P))$values)), decreasing = TRUE)
  eigs[2]
}

cat("\n── λ₂ by childhood attendance stratum ──\n")
for (att in names(P_att)) {
  cat(sprintf("  %-10s  N = %4d  λ₂ = %.4f\n",
              att, n_att[[att]], lambda2(P_att[[att]])))
}

cat("\n── Diagonal persistence P[i→i] by stratum ──\n")
diag_tbl = do.call(rbind, lapply(names(P_att), function(att) {
  data.frame(
    stratum = att,
    origin  = names(diag(P_att[[att]])),
    persist = round(diag(P_att[[att]]), 3),
    n       = n_att[[att]],
    row.names = NULL
  )
}))
print(diag_tbl[order(diag_tbl$origin, diag_tbl$stratum), ])

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/attend12", recursive = TRUE, showWarnings = FALSE)

# Per-stratum transition heatmaps (P + π₀ + π*)
for (att in names(P_att)) {
  p = make_combined(
    P_att[[att]], pi0_att[[att]], pistar_att[[att]],
    levels    = rel_level_order,
    title_str = paste0(att_labels[att], "  (N = ", n_att[[att]], ")")
  )
  ggsave(paste0("output/figures/attend12/trans_", att, ".png"),
         p, width = 10, height = 7, dpi = 200)
}

# Diagonal persistence across attendance strata
diag_att = do.call(rbind, lapply(names(P_att), function(att) {
  data.frame(
    stratum = factor(att_labels[att], levels = att_labels),
    origin  = factor(rel_level_order, levels = rel_level_order),
    persist = diag(P_att[[att]])[rel_level_order],
    row.names = NULL
  )
}))

p_diag = ggplot(diag_att, aes(x = stratum, y = persist,
                               color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  labs(x = "Childhood attendance (age 12)",
       y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by Childhood Attendance Level") +
  healy_theme

ggsave("output/figures/attend12/diagonal_persistence_attend12.png",
       p_diag, width = 9, height = 5, dpi = 200)

# λ₂ bar chart across strata
lambda2_df = data.frame(
  stratum = factor(att_labels[names(P_att)], levels = att_labels),
  lambda2 = sapply(P_att, lambda2),
  n       = unlist(n_att)
)

p_lambda = ggplot(lambda2_df, aes(x = stratum, y = lambda2)) +
  geom_col(fill = "#0072B2", width = 0.4) +
  geom_text(aes(label = sprintf("%.3f\n(N=%d)", lambda2, n)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Childhood attendance (age 12)",
       y = "λ₂ (second eigenvalue)",
       title = "λ₂ by Childhood Attendance Level") +
  healy_theme

ggsave("output/figures/attend12/lambda2_attend12.png",
       p_lambda, width = 7, height = 5, dpi = 200)

cat("\nWrote output/figures/attend12/\n")
