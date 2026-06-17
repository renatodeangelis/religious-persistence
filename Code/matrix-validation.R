library(dplyr)
library(tidyr)
library(ggplot2)
source("Code/utils.R")

# ---------------------------------------------------------------------------
# Load Add Health data (Wave 1 and Wave 3 in-home surveys)
# ---------------------------------------------------------------------------

dir_data = "add-health"

addhealth_w1 = readRDS(file.path(dir_data, "w1inhome.rds"))
addhealth_w3 = readRDS(file.path(dir_data, "w3inhome.rds"))

glimpse(addhealth_w1)
glimpse(addhealth_w3)

# ---------------------------------------------------------------------------
# Extract religion variables and unique identifier
# W1: AID (respondent ID), H1RE1 (respondent religion), PA22 (parent religion)
# W3: AID, H3RE1 (respondent religion at Wave 3), H3RE26 (recalled childhood religion)
# ---------------------------------------------------------------------------

w1_relig = addhealth_w1 |>
  select(AID, H1RE1, PA22)

w3_relig = addhealth_w3 |>
  select(AID, H3RE1, H3RE2, H3RE6, H3RE26)

addhealth = inner_join(w1_relig, w3_relig, by = "AID")

# ---------------------------------------------------------------------------
# Recode PA22, H3RE1, and H3RE26 to a 4-category scheme:
#   "None", "Protestant", "Catholic", "Other"
# H1RE1 (respondent's own W1 religion) uses the same scheme for concordance checks.
# Jewish, Buddhist, Hindu, Islam etc. all collapse to "Other".
# Variables are ordered factors with labels like "(22) Catholic"; extract the
# leading numeric code into temporary columns, recode, then drop temporaries.
# ---------------------------------------------------------------------------

relig_lvls = c("none", "other", "protestant", "catholic")

addhealth = addhealth |>
  mutate(
    # Extract numeric codes from factor labels
    .h1re1  = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H1RE1))),
    .pa22   = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(PA22))),
    .h3re1  = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H3RE1))),
    .h3re2  = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H3RE2))),
    .h3re26 = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H3RE26))),

    # H1RE1: respondent's own religion at Wave 1
    # 0=None, 1–19=Protestant denominations, 22=Catholic, 27=Unitarian → Protestant
    # 20=Baha'i, 21=Buddhist, 23=E. Orthodox, 24=Hindu, 25=Islam, 26=Jewish, 28=Other → Other
    h1re1_relig = case_when(
      .h1re1 == 0                               ~ "none",
      .h1re1 %in% c(1:19, 27)                   ~ "protestant",
      .h1re1 == 22                               ~ "catholic",
      .h1re1 %in% c(20, 21, 23, 24, 25, 26, 28) ~ "other",
      TRUE                                       ~ NA_character_
    ),

    # PA22: parent's self-reported religion at Wave 1
    # 7=Catholic; 28=None; Protestant denominations → Protestant
    # 6=Buddhist, 11=E. Orthodox, 14=Hindu, 16=Islam, 18=Jewish, 24=Other → Other
    pa22_relig = case_when(
      .pa22 == 28                                            ~ "none",
      .pa22 %in% c(1:3, 5, 8:10, 12:13, 15, 20, 21,
                   23, 25:27)                               ~ "protestant",
      .pa22 == 7                                             ~ "catholic",
      .pa22 %in% c(6, 11, 14, 16:19, 24)                   ~ "other",
      TRUE                                                   ~ NA_character_
    ),

    # H3RE1: respondent's current religion at Wave 3
    h3re1_relig = case_when(
      .h3re1 == 0 | (.h3re1 == 7 & .h3re2 %in% c(0, 11)) ~ "none",
      .h3re1 %in% c(1, 8) | (.h3re1 == 7 & .h3re2 == 1) ~ "protestant",
      .h3re1 == 2 | (.h3re1 == 7 & .h3re2 == 2) ~ "catholic",
      .h3re1 %in% c(3, 4, 5, 6, 7)  ~ "other",
      TRUE ~ NA_character_
    ),

    # H3RE26: retrospective recall of childhood religion at Wave 3
    h3re26_relig = case_when(
      .h3re26 == 0                    ~ "none",
      .h3re26 %in% c(1, 8)           ~ "protestant",
      .h3re26 == 2                    ~ "catholic",
      .h3re26 %in% c(3, 4, 5, 6, 7)  ~ "other",
      TRUE                            ~ NA_character_
    )
  ) |>
  select(-.h1re1, -.pa22, -.h3re1, -.h3re2, -.h3re26)

# ---------------------------------------------------------------------------
# PA22 missingness diagnostic
# Cross-tabulate child's own W1 religion (always observed) by whether
# parental religion (PA22) is missing — to assess whether missingness is
# random with respect to religious background
# ---------------------------------------------------------------------------

miss_tab = addhealth |>
  mutate(
    pa22_missing  = is.na(pa22_relig),
    child_relig   = factor(h1re1_relig, levels = relig_lvls)
  ) |>
  count(child_relig, pa22_missing) |>
  group_by(pa22_missing) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup() |>
  mutate(pa22_missing = if_else(pa22_missing, "PA22 missing", "PA22 observed"))

print(miss_tab |> tidyr::pivot_wider(names_from = pa22_missing,
                                     values_from = c(n, pct),
                                     names_glue = "{pa22_missing}: {.value}"))

# ---------------------------------------------------------------------------
# Recall bias cross-tabulation
# Rows: actual parental religion (PA22, ground truth)
# Cols: respondent's retrospective recall of childhood religion (H3RE26)
# Complete cases only (both PA22 and H3RE26 observed)
# Off-diagonal cells = recall mismatch
# ---------------------------------------------------------------------------

recall_dat = addhealth |>
  filter(!is.na(pa22_relig), !is.na(h3re26_relig)) |>
  mutate(
    parent  = factor(pa22_relig,    levels = relig_lvls),
    recall  = factor(h3re26_relig,  levels = relig_lvls)
  )

# Raw counts
recall_tab = table(parent = recall_dat$parent, recall = recall_dat$recall)
cat("\nRecall bias table — counts (rows = actual PA22, cols = W3 recall):\n")
print(recall_tab)

# Row-proportions: conditional on actual parental religion, what did they recall?
cat("\nRecall bias table — row proportions:\n")
print(round(prop.table(recall_tab, margin = 1), 3))

# Overall and category-level mismatch rates
recall_dat = recall_dat |>
  mutate(mismatch = parent != recall)

cat("\nOverall mismatch rate:", round(mean(recall_dat$mismatch), 3), "\n")

cat("\nMismatch rate by actual parental religion:\n")
recall_dat |>
  group_by(parent) |>
  summarise(n = n(), mismatch_rate = round(mean(mismatch), 3)) |>
  print()

# ---------------------------------------------------------------------------
# Recall bias table — concordant-at-W1 subsample
# Restricted to respondents whose own Wave 1 religion (relig_child_w1)
# matched their parent's Wave 1 religion (relig_parent_w1). In these
# households there is no ambiguity about the childhood religious environment,
# so off-diagonal cells are more clearly attributable to recall error rather
# than genuine household heterogeneity.
# ---------------------------------------------------------------------------

recall_concordant = addhealth |>
  filter(
    !is.na(pa22_relig), !is.na(h3re26_relig),
    !is.na(h1re1_relig),
    h1re1_relig == pa22_relig
  ) |>
  mutate(
    parent  = factor(pa22_relig,   levels = relig_lvls),
    recall  = factor(h3re26_relig, levels = relig_lvls)
  )

cat("\nConcordant-at-W1 subsample: n =", nrow(recall_concordant), "\n")
cat("(respondents whose own W1 religion matched parent's W1 religion)\n")

recall_tab_conc = table(parent = recall_concordant$parent,
                        recall = recall_concordant$recall)

cat("\nRecall bias table (concordant-at-W1) — counts:\n")
print(recall_tab_conc)

cat("\nRecall bias table (concordant-at-W1) — row proportions:\n")
print(round(prop.table(recall_tab_conc, margin = 1), 3))

recall_concordant = recall_concordant |>
  mutate(mismatch = parent != recall)

cat("\nOverall mismatch rate (concordant-at-W1):",
    round(mean(recall_concordant$mismatch), 3), "\n")

cat("\nMismatch rate by actual parental religion (concordant-at-W1):\n")
recall_concordant |>
  group_by(parent) |>
  summarise(n = n(), mismatch_rate = round(mean(mismatch), 3)) |>
  print()

# ---------------------------------------------------------------------------
# Child–parent religious discordance at Wave 1
# How many respondents already differed from their parent's religion at age ~16?
# relig_child_w1 (H1RE1) vs relig_parent_w1 (PA22)
# This matters because PA22 is only a clean ground truth if the household
# shared a single religion — prior divergence inflates apparent recall mismatch
# ---------------------------------------------------------------------------

discord_dat = addhealth |>
  filter(!is.na(h1re1_relig), !is.na(pa22_relig)) |>
  mutate(
    child   = factor(h1re1_relig, levels = relig_lvls),
    parent  = factor(pa22_relig,  levels = relig_lvls),
    discord = child != parent
  )

cat("\nChild–parent discordance at Wave 1:\n")
cat("N (complete cases):", nrow(discord_dat), "\n")
cat("Overall discordance rate:", round(mean(discord_dat$discord), 3), "\n")

cat("\nDiscordance rate by parental religion:\n")
discord_dat |>
  group_by(parent) |>
  summarise(n = n(), discord_rate = round(mean(discord), 3)) |>
  print()

cat("\nChild–parent discordance table (rows = parent, cols = child at W1):\n")
print(table(parent = discord_dat$parent, child = discord_dat$child))

# ---------------------------------------------------------------------------
# Adult religion and direction of recall mismatch (endogeneity test)
# For each current W3 religion, what share of respondents recalled:
#   (a) accurately (recall == parent)
#   (b) more secular than parent (recall = None, parent != None)
#   (c) more religious than parent (parent = None, recall != None)
#   (d) a different religion (both non-None but differ)
# The Hayward pattern predicts (b) is concentrated among current Nones
# ---------------------------------------------------------------------------

endogeneity_dat = addhealth |>
  filter(!is.na(pa22_relig), !is.na(h3re26_relig),
         !is.na(h3re1_relig)) |>
  mutate(
    adult  = factor(h3re1_relig,  levels = relig_lvls),
    parent = factor(pa22_relig,   levels = relig_lvls),
    recall = factor(h3re26_relig, levels = relig_lvls),
    mismatch_type = case_when(
      recall == parent                          ~ "Concordant",
      recall == "none" & parent != "none"       ~ "Recalled more secular",
      parent == "none" & recall != "none"       ~ "Recalled more religious",
      TRUE                                      ~ "Recalled different religion"
    )
  )

cat("\nMismatch type by current adult religion (W3) — row percentages:\n")
endogeneity_dat |>
  count(adult, mismatch_type) |>
  group_by(adult) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup() |>
  select(adult, mismatch_type, pct) |>
  tidyr::pivot_wider(names_from = mismatch_type, values_from = pct) |>
  print()

# ---------------------------------------------------------------------------
# Logistic regression: does current adult religion predict recalling
# a more secular upbringing than parents reported?
#
# Outcome: binary — recalled None when parent was non-None ("more secular")
# Sample: complete cases where parent != None (outcome is only possible here)
# Reference category: Protestant for both adult and parental religion
# Controls: parental religion dummies absorb origin-state differences in
#           base rates of misclassification
# ---------------------------------------------------------------------------

logreg_dat = endogeneity_dat |>
  filter(parent != "none") |>
  mutate(
    more_secular = as.integer(mismatch_type == "Recalled more secular"),
    adult_f  = relevel(droplevels(adult),  ref = "protestant"),
    parent_f = relevel(droplevels(parent), ref = "protestant")
  )

mod = glm(more_secular ~ adult_f + parent_f,
          data   = logreg_dat,
          family = binomial(link = "logit"))

coef_tab = coef(summary(mod))
or_tab   = data.frame(
  term  = rownames(coef_tab),
  OR    = round(exp(coef_tab[, "Estimate"]), 3),
  CI_lo = round(exp(coef_tab[, "Estimate"] - 1.96 * coef_tab[, "Std. Error"]), 3),
  CI_hi = round(exp(coef_tab[, "Estimate"] + 1.96 * coef_tab[, "Std. Error"]), 3),
  p     = round(coef_tab[, "Pr(>|z|)"], 4)
)

cat("\nLogistic regression — outcome: recalled more secular than parent\n")
cat("Reference: Protestant (adult and parent) | n =", nrow(logreg_dat), "\n\n")
print(or_tab, row.names = FALSE)

# ---------------------------------------------------------------------------
# Model 2: add respondent's own religion at age ~16 (relig_child_w1)
# Tests how much of the None effect is explained by having already left
# religion by Wave 1, vs. adult identity independently shaping recall
# ---------------------------------------------------------------------------

logreg_dat2 = logreg_dat |>
  filter(!is.na(h1re1_relig)) |>
  mutate(child_f = relevel(
    factor(h1re1_relig, levels = relig_lvls),
    ref = "protestant"
  ))

# Refit baseline on same restricted sample for a clean comparison
mod_base2 = glm(more_secular ~ adult_f + parent_f,
                data   = logreg_dat2,
                family = binomial(link = "logit"))

mod_child  = glm(more_secular ~ adult_f + parent_f + child_f,
                 data   = logreg_dat2,
                 family = binomial(link = "logit"))

fmt_or = function(mod) {
  ct  = coef(summary(mod))
  data.frame(
    term  = rownames(ct),
    OR    = round(exp(ct[, "Estimate"]), 3),
    CI_lo = round(exp(ct[, "Estimate"] - 1.96 * ct[, "Std. Error"]), 3),
    CI_hi = round(exp(ct[, "Estimate"] + 1.96 * ct[, "Std. Error"]), 3),
    p     = round(ct[, "Pr(>|z|)"], 4)
  )
}

cat("\nModel 1 (no child W1 religion) | n =", nrow(logreg_dat2), "\n\n")
print(fmt_or(mod_base2), row.names = FALSE)

cat("\nModel 2 (+ respondent religion at age ~16) | n =", nrow(logreg_dat2), "\n\n")
print(fmt_or(mod_child), row.names = FALSE)

cat("\nAttenuation in None OR:",
    round(exp(coef(mod_base2)["adult_fnone"]), 3), "->",
    round(exp(coef(mod_child)["adult_fnone"]), 3), "\n")

# ---------------------------------------------------------------------------
# Transition matrices: PA22 → H3RE1 and H3RE26 → H3RE1
# Compare "true" parent-to-child transmission (PA22) against the
# recall-based origin measure (H3RE26) as a diagnostic for recall bias
# ---------------------------------------------------------------------------

dir.create("output/figures/validation", recursive = TRUE, showWarnings = FALSE)

# Matrix 1: actual parental religion (PA22) → respondent's adult religion (W3)
pmat_dat_pa22 = addhealth |>
  filter(!is.na(pa22_relig), !is.na(h3re1_relig))

P_pa22   = p_matrix(pmat_dat_pa22, "pa22_relig",  "h3re1_relig", levels = relig_lvls)
pi0_pa22 = pi_0(pmat_dat_pa22, "pa22_relig")
pis_pa22 = pi_star(P_pa22)

# Matrix 2: recalled childhood religion (H3RE26) → respondent's adult religion (W3)
pmat_dat_h3re26 = addhealth |>
  filter(!is.na(h3re26_relig), !is.na(h3re1_relig))

P_h3re26   = p_matrix(pmat_dat_h3re26, "h3re26_relig", "h3re1_relig", levels = relig_lvls)
pi0_h3re26 = pi_0(pmat_dat_h3re26, "h3re26_relig")
pis_h3re26 = pi_star(P_h3re26)

p_pa22 = make_combined(
  P_pa22, pi0_pa22, pis_pa22,
  levels    = relig_lvls,
  title_str = "PA22 → H3RE1\n(Actual parent religion W1 → Adult religion W3)"
)

p_h3re26 = make_combined(
  P_h3re26, pi0_h3re26, pis_h3re26,
  levels    = relig_lvls,
  title_str = "H3RE26 → H3RE1\n(Recalled childhood religion W3 → Adult religion W3)"
)

ggsave("output/figures/validation/pmat_pa22_h3re1.png",
       p_pa22,   width = 10, height = 7, dpi = 200)
ggsave("output/figures/validation/pmat_h3re26_h3re1.png",
       p_h3re26, width = 10, height = 7, dpi = 200)