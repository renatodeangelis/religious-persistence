# ── 03 · SAMPLE DIAGNOSTICS ────────────────────────────────────────────────────
# Cell-count and coverage diagnostics that justify the state space and the n < 30
# exclusion rule used in the estimation loops (02). Prints three tables to the
# console; no figures or saved artifacts.
#
# Input: data/derived/gss_clean.rds

library(dplyr)
source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

# ── COHORT SAMPLE SIZES (5/10-year windows) ──────────────────────────────────
# N per cohort window feeding the main RELIG16 -> RELIG transition matrices in
# 02 (same non-missing filter as the P_list_* loops). The n < 30 threshold
# matches the exclusion rule used in those loops.

valid_rows = !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt)

cohort_n_table = function(cohort_var, label) {
  coh = data[[cohort_var]]
  # upper bound widened to 1990 so binned-midpoint labels (10-yr up to 1985)
  # are not clipped; the 1-year window is unaffected.
  keep = valid_rows & !is.na(coh) & coh >= 1920 & coh <= 1995
  tab = table(coh[keep])
  data.frame(
    window   = label,
    cohort   = as.numeric(names(tab)),   # midpoint bins are non-integer (e.g. 1922.5)
    n        = as.integer(tab),
    below_30 = as.integer(tab) < 30,
    row.names = NULL
  )
}

cohort_n_df = do.call(rbind, list(
  cohort_n_table("cohort",    "1-year"),
  cohort_n_table("cohort_5",  "5-year"),
  cohort_n_table("cohort_10", "10-year")
))

print(cohort_n_df, row.names = FALSE)
