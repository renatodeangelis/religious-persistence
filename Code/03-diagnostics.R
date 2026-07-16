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

# ── COHORT SAMPLE SIZES (5/10/20-year windows) ───────────────────────────────
# N per cohort window feeding the main RELIG16 -> RELIG transition matrices in
# 02 (same non-missing filter as the P_list_* loops). The n < 30 threshold
# matches the exclusion rule used in those loops.

valid_rows = !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt)

cohort_n_table = function(cohort_var, label) {
  coh = data[[cohort_var]]
  keep = valid_rows & !is.na(coh) & coh >= 1920 & coh <= 1980
  tab = table(coh[keep])
  data.frame(
    window   = label,
    cohort   = as.integer(names(tab)),
    n        = as.integer(tab),
    below_30 = as.integer(tab) < 30,
    row.names = NULL
  )
}

cohort_n_df = do.call(rbind, list(
  cohort_n_table("cohort",    "1-year"),
  cohort_n_table("cohort_5",  "5-year"),
  cohort_n_table("cohort_10", "10-year"),
  cohort_n_table("cohort_20", "20-year")
))

print(cohort_n_df, row.names = FALSE)

# ── ATTITUDE COVERAGE SUMMARY ─────────────────────────────────────────────────

att_coverage = lapply(att_vars, function(v) {
  data |>
    filter(!is.na(reltrad_alt), !is.na(reltrad16_alt)) |>
    summarise(
      variable  = v,
      n_liberal = sum(.data[[v]] == 1L, na.rm = TRUE),
      n_conserv = sum(.data[[v]] == 0L, na.rm = TRUE),
      n_missing = sum(is.na(.data[[v]])),
      pct_cover = round((n_liberal + n_conserv) / n() * 100, 1)
    )
}) |> bind_rows()

print(att_coverage, n = Inf)

# ── CELL-COUNT FEASIBILITY (20-year cohort windows) ──────────────────────────

cell_rows = list()

for (v in att_vars) {
  for (grp in c(0L, 1L)) {
    for (coh in c(1920, 1940, 1960, 1980)) {
      sub = data[
        !is.na(data$cohort_20) & data$cohort_20 == coh &
        !is.na(data[[v]])      & data[[v]]       == grp, ]
      if (nrow(sub) == 0) next
      tab = table(
        factor(sub$reltrad16_alt, levels = states_alt),
        factor(sub$reltrad_alt,   levels = states_alt)
      )
      cell_rows[[length(cell_rows) + 1]] = data.frame(
        variable  = v,
        group     = if (grp == 1L) "liberal" else "conservative",
        cohort    = coh,
        n         = nrow(sub),
        min_cell  = min(tab),
        cells_lt5 = sum(tab < 5),
        cells_0   = sum(tab == 0),
        row.names = NULL
      )
    }
  }
}

cell_df = do.call(rbind, cell_rows)

print(cell_df, row.names = FALSE)
