# ── 01 · DATA PREPARATION ─────────────────────────────────────────────────────
# Loads GSS from the gssr package, builds and cleans the master analysis frame
# with every derived column (reltrad recodes, cohort bins, party/polviews
# binaries), and the state space.
#
# Output: data/derived/gss_clean.rds — a list(data, states_alt) consumed by the
# estimation script (02) and every analysis script that recomputes matrices at
# 1-year resolution (03, 05, 06, 07).

library(dplyr)
library(tidyr)
library(gssr)
source("code/utils.R")

reltrad_labels = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

# ── DATA ──────────────────────────────────────────────────────────────────────

data(gss_all)
data = gss_all |>
  select(year, cohort, sex, reltrad, reltrad16, region, born,
         race, polviews, partyid) |>
  filter(!(year %in% c(1972, 1982, 1987, 2021, 2022, 2024))) |>
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
  # 6-state variant for the Black-Protestant robustness stage (12): collapse
  # only jewish → other, keeping "black protestant" as its own state.
  mutate(across(c(reltrad, reltrad16),
                ~ if_else(. == "jewish", "other", .),
                .names = "{.col}_bp")) |>
  # cohort arrives from gss_all as a haven_labelled vector; strip to plain
  # numeric so downstream median()/binning behave (median.haven_labelled errors)
  mutate(cohort = as.numeric(cohort)) |>
  mutate(age = year - cohort) |>
  filter(age >= 30, age <= 75, cohort >= 1925, cohort <= 1994) |>
  mutate(
    cohort_5   = (floor((cohort - 1900) / 5)  * 5  + 1900) + 2.5,
    cohort_10  = (floor((cohort - 1900) / 10) * 10 + 1900) + 5,
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

states_alt = sort(unique(c(data$reltrad_alt, data$reltrad16_alt)))
states_alt = states_alt[!is.na(states_alt)]

# 6-state space for the Black-Protestant robustness stage (12)
states_bp = sort(unique(c(data$reltrad_bp, data$reltrad16_bp)))
states_bp = states_bp[!is.na(states_bp)]

# ── SAVE ──────────────────────────────────────────────────────────────────────
# Strip haven value labels so the saved frame is plain numeric/character. The
# gssr columns arrive as haven_labelled, and as.numeric() on those only works
# while haven's S3 methods are attached (they are here, via library(gssr), but
# not in the downstream scripts that read this file).

data = haven::zap_labels(data)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = data, states_alt = states_alt, states_bp = states_bp),
        "data/derived/gss_clean.rds")
cat("Wrote data/derived/gss_clean.rds  (", nrow(data), "rows )\n")
