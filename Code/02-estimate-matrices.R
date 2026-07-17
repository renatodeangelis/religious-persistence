# ── 02 · MATRIX ESTIMATION ────────────────────────────────────────────────────
# Builds every count (N), probability (P), initial (pi0), and stationary (pistar)
# matrix used downstream, across all stratifications, and persists them in one
# nested list. Consumers (04–07) load this instead of recomputing.
#
# Raw count matrices (N) are the source of truth — P is just row-normalized N.
# They are saved ONLY where a consumer needs them: the national 5- and 10-year
# counts feed the Anderson-Goodman homogeneity tests in 04. Add N to other strata
# only if you later want homogeneity tests within them.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices.rds

source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

# ── ATTITUDE MATRICES: POOLED (no cohort stratification) ─────────────────────
# Mean birth year is tabulated to show cohort composition of each stratum.
# Note: item coverage varies by GSS year — inspect mean_cohort for drift.

P_list_att_pooled   = list()
pi0_list_att_pooled = list()
pooled_summary_rows = list()

for (v in att_vars) {
  for (grp in c(0L, 1L)) {
    sub = data[!is.na(data$reltrad_alt) & !is.na(data$reltrad16_alt) &
               !is.na(data[[v]])        & data[[v]] == grp, ]
    if (nrow(sub) < 30) next
    grp_lbl = if (grp == 1L) "liberal" else "conservative"
    key     = paste(v, grp_lbl, sep = "_")

    P_list_att_pooled[[key]]   = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0_list_att_pooled[[key]] = pi_0(sub, "reltrad16_alt")

    pooled_summary_rows[[key]] = data.frame(
      variable      = v,
      group         = grp_lbl,
      n             = nrow(sub),
      mean_cohort   = round(mean(sub$cohort,   na.rm = TRUE), 1),
      median_cohort = median(sub$cohort, na.rm = TRUE),
      row.names     = NULL
    )
  }
}

pooled_summary = do.call(rbind, pooled_summary_rows)

# ── NATIONAL COHORT MATRICES ─────────────────────────────────────────────────

# ── 5-year cohort loop ───────────────────────────────────────────────────────
cohorts_5 = sort(unique(data$cohort_5[!is.na(data$cohort_5) & data$cohort_5 >= 1922.5 & data$cohort_5 <= 1982.5]))

P_list_5      = list()
pi0_list_5    = list()
pistar_list_5 = list()
N_list_5      = list()   # raw counts feeding the rolling-window sensitivity (04)

for (coh in cohorts_5) {
  sub = data[!is.na(data$cohort_5) & data$cohort_5 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_5[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_5[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_5[[key]] = pi_star(P_list_5[[key]])
  N_list_5[[key]]      = count_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
}

# ── 10-year cohort loop ──────────────────────────────────────────────────────
cohorts_10 = sort(unique(data$cohort_10[!is.na(data$cohort_10) & data$cohort_10 >= 1925 & data$cohort_10 <= 1985]))

P_list_10      = list()
pi0_list_10    = list()
pistar_list_10 = list()
N_list_10      = list()   # raw count matrices feeding the homogeneity test (04)

for (coh in cohorts_10) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_10[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_10[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_10[[key]] = pi_star(P_list_10[[key]])
  N_list_10[[key]]      = count_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
}

# ── 20-year cohort loop ──────────────────────────────────────────────────────
cohorts_20_pooled = sort(unique(data$cohort_20[!is.na(data$cohort_20) & data$cohort_20 >= 1930 & data$cohort_20 <= 1990]))

P_list_20      = list()
pi0_list_20    = list()
pistar_list_20 = list()

for (coh in cohorts_20_pooled) {
  sub = data[!is.na(data$cohort_20) & data$cohort_20 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_20[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_20[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_20[[key]] = pi_star(P_list_20[[key]])
}

# ── REGIONAL COHORT MATRICES ─────────────────────────────────────────────────

cohorts_20    = c(1930, 1950, 1970)   # 20-year bin midpoints (edges 1920/1940/1960)
regions_broad = c("Midwest", "Northeast", "South", "West")

P_list_reg      = list()
pi0_list_reg    = list()
pistar_list_reg = list()

for (reg in regions_broad) {
  for (coh in cohorts_20) {
    sub = data[
      !is.na(data$cohort_20)     & data$cohort_20     == coh &
      !is.na(data$region_broad)  & data$region_broad  == reg &
      !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
    if (nrow(sub) < 30) next
    key = paste(reg, coh, sep = "_")

    P_list_reg[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0_list_reg[[key]]    = pi_0(sub, "reltrad16_alt")
    pistar_list_reg[[key]] = pi_star(P_list_reg[[key]])
  }
}

# ── BINARY (AFFILIATED / UNAFFILIATED) 2×2 MATRICES ──────────────────────────

states_2x2 = c("affiliated", "unaffiliated")

cohorts_10_2x2 = sort(unique(
  data$cohort_10[!is.na(data$cohort_10) & data$cohort_10 >= 1925 & data$cohort_10 <= 1985]
))

P_list_2x2      = list()
pi0_list_2x2    = list()
pistar_list_2x2 = list()
n_list_2x2      = list()

for (coh in cohorts_10_2x2) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == coh &
               !is.na(data$belief16) & !is.na(data$belief), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_2x2[[key]]      = p_matrix(sub, "belief16", "belief", levels = states_2x2)
  pi0_list_2x2[[key]]    = pi_0(sub, "belief16")
  pistar_list_2x2[[key]] = pi_star(P_list_2x2[[key]])
  n_list_2x2[[key]]      = nrow(sub)
}

# ── NATIVITY-SPLIT MATRICES (10-year cohorts: 1950, 1960, 1970) ──────────────

cohorts_nat    = c(1955, 1965, 1975)   # 10-year bin midpoints (edges 1950/1960/1970)
nativity_groups = c("Born in US", "Born abroad")

P_list_nat      = list()
pi0_list_nat    = list()
pistar_list_nat = list()
n_list_nat      = list()

for (nat in nativity_groups) {
  for (coh in cohorts_nat) {
    sub = data[
      !is.na(data$cohort_10)     & data$cohort_10 == coh &
      !is.na(data$nativity)      & data$nativity   == nat &
      !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
    if (nrow(sub) < 30) next
    key = paste(gsub(" ", "_", nat), coh, sep = "_")

    P_list_nat[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0_list_nat[[key]]    = pi_0(sub, "reltrad16_alt")
    pistar_list_nat[[key]] = pi_star(P_list_nat[[key]])
    n_list_nat[[key]]      = nrow(sub)
  }
}

# ── SEX-STRATIFIED DECADAL MATRICES (10-year cohorts, 1940–1980) ─────────────

cohorts_sex = c(1945, 1955, 1965, 1975, 1985)   # 10-year bin midpoints (edges 1940–1980)
sex_labels  = c("1" = "male", "2" = "female")

P_list_sex      = list()
pi0_list_sex    = list()
pistar_list_sex = list()
n_list_sex      = list()

for (sx in c(1, 2)) {
  for (coh in cohorts_sex) {
    sub = data[
      !is.na(data$cohort_10)     & data$cohort_10         == coh &
      !is.na(data$sex)           & as.numeric(data$sex)   == sx  &
      !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
    if (nrow(sub) < 30) next
    key = paste(sex_labels[as.character(sx)], coh, sep = "_")

    P_list_sex[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0_list_sex[[key]]    = pi_0(sub, "reltrad16_alt")
    pistar_list_sex[[key]] = pi_star(P_list_sex[[key]])
    n_list_sex[[key]]      = nrow(sub)
  }
}

# ── POLITICAL STRATIFICATION DECADAL MATRICES (10-year cohorts, 1940–1989) ───

pol_vars = list(
  partyid_narrow  = c("dem", "rep", "other"),
  partyid_broad   = c("dem", "rep", "other"),
  polviews_narrow = c("liberal", "moderate", "conservative"),
  polviews_broad  = c("liberal", "moderate", "conservative")
)

cohorts_pol = c(1945, 1955, 1965, 1975, 1985)   # 10-year bin midpoints (edges 1940–1980)

P_list_pol      = list()
pi0_list_pol    = list()
pistar_list_pol = list()
n_list_pol      = list()

for (vname in names(pol_vars)) {
  for (grp in pol_vars[[vname]]) {
    for (coh in cohorts_pol) {
      sub = data[
        !is.na(data$cohort_10)     & data$cohort_10 == coh &
        !is.na(data[[vname]])      & data[[vname]]  == grp &
        !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
      if (nrow(sub) < 30) next
      key = paste(vname, grp, coh, sep = "_")
      P_list_pol[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
      pi0_list_pol[[key]]    = pi_0(sub, "reltrad16_alt")
      pistar_list_pol[[key]] = pi_star(P_list_pol[[key]])
      n_list_pol[[key]]      = nrow(sub)
    }
  }
}

# ── SAVE ──────────────────────────────────────────────────────────────────────

matrices = list(
  att_pooled = list(P = P_list_att_pooled, pi0 = pi0_list_att_pooled,
                    summary = pooled_summary),
  nat5  = list(P = P_list_5,  pi0 = pi0_list_5,  pistar = pistar_list_5,  N = N_list_5),
  nat10 = list(P = P_list_10, pi0 = pi0_list_10, pistar = pistar_list_10, N = N_list_10),
  nat20 = list(P = P_list_20, pi0 = pi0_list_20, pistar = pistar_list_20),
  region    = list(P = P_list_reg, pi0 = pi0_list_reg, pistar = pistar_list_reg),
  binary    = list(P = P_list_2x2, pi0 = pi0_list_2x2, pistar = pistar_list_2x2,
                   n = n_list_2x2, states = states_2x2),
  nativity  = list(P = P_list_nat, pi0 = pi0_list_nat, pistar = pistar_list_nat,
                   n = n_list_nat),
  sex       = list(P = P_list_sex, pi0 = pi0_list_sex, pistar = pistar_list_sex,
                   n = n_list_sex),
  political = list(P = P_list_pol, pi0 = pi0_list_pol, pistar = pistar_list_pol,
                   n = n_list_pol)
)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
saveRDS(matrices, "data/derived/matrices.rds")
cat("Wrote data/derived/matrices.rds\n")
