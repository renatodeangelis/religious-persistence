# ── DIAGNOSTIC: Per-cell counts for stratified transition matrices ─────────────
# Checks whether all 6 cohort windows (1925–1984) support transition matrix
# estimation within each stratum × origin-state cell.
#
# Flags any cell with N < THRESHOLD (default 30) — a flagged row means that
# row of the transition matrix is unreliable.
#
# Run manually. No output files written.
#
# Stratifications checked: sex, nonblack, partyid (narrow/broad),
#                          polviews (narrow/broad)

library(dplyr)
library(tidyr)

clean  = readRDS("data/derived/gss_clean.rds")
data   = clean$data

THRESHOLD     = 30
mids_all      = c(1930, 1940, 1950, 1960, 1970, 1980)
cohort_labels = setNames(
  c("1925-34", "1935-44", "1945-54", "1955-64", "1965-74", "1975-84"),
  as.character(mids_all)
)

# ── HELPER ────────────────────────────────────────────────────────────────────
# df must have a 'stratum' column. Prints wide count table (cohort windows as
# columns) and a feasibility summary.

summarise_cut = function(df, label) {
  df = df |>
    filter(cohort_10 %in% mids_all,
           !is.na(stratum), !is.na(reltrad16_alt), !is.na(reltrad_alt)) |>
    group_by(cohort_10, stratum, reltrad16_alt) |>
    summarise(n = n(), .groups = "drop") |>
    mutate(
      cohort = cohort_labels[as.character(cohort_10)],
      cell   = ifelse(n < THRESHOLD, paste0(n, "*"), as.character(n))
    )

  cat("\n", strrep("=", 72), "\n")
  cat(" ", label, "\n")
  cat(strrep("=", 72), "\n\n")

  wide = df |>
    select(stratum, reltrad16_alt, cohort, cell) |>
    pivot_wider(names_from = cohort, values_from = cell, values_fill = "—") |>
    rename(origin = reltrad16_alt) |>
    arrange(stratum, origin)
  print(as.data.frame(wide), row.names = FALSE)

  # Window-level summary
  cat("\nFlagged cells per cohort window (N <", THRESHOLD, ", marked with *):\n")
  win_sum = df |>
    group_by(cohort) |>
    summarise(flagged = sum(n < THRESHOLD), total = n(), .groups = "drop")
  print(as.data.frame(win_sum), row.names = FALSE)

  # Drop-earliest / drop-latest advice
  n_total = sum(df$n < THRESHOLD)
  if (n_total == 0) {
    cat("\n  All 6 windows viable across all strata and origin states.\n")
  } else {
    n_keep_early = sum(df$n[df$cohort_10 != min(mids_all)] < THRESHOLD)
    n_keep_late  = sum(df$n[df$cohort_10 != max(mids_all)] < THRESHOLD)
    cat("\n  Total flagged cells:", n_total, "\n")
    cat("  Drop 1925-34: resolves", n_total - n_keep_early,
        "flagged cell(s);", n_keep_early, "remain.\n")
    cat("  Drop 1975-84: resolves", n_total - n_keep_late,
        "flagged cell(s);", n_keep_late, "remain.\n")
  }
}

# ── SEX ───────────────────────────────────────────────────────────────────────

data |>
  mutate(stratum = case_when(
    as.numeric(sex) == 1 ~ "male",
    as.numeric(sex) == 2 ~ "female"
  )) |>
  summarise_cut("Sex: male vs female")

# ── NON-BLACK ─────────────────────────────────────────────────────────────────

data |>
  mutate(stratum = ifelse(as.numeric(race) != 2, "non-Black", "Black")) |>
  summarise_cut("Race: non-Black vs Black")

# ── POLITICAL ─────────────────────────────────────────────────────────────────

pol_vars = c("partyid_narrow", "partyid_broad", "polviews_narrow", "polviews_broad")

for (vname in pol_vars) {
  data |>
    mutate(stratum = .data[[vname]]) |>
    summarise_cut(paste("Political:", gsub("_", " ", vname)))
}
