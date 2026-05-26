library(dplyr)
library(tidyr)
library(ggplot2)
library(haven)
library(janitor)

# ---------------------------------------------------------------------------
# Load Add Health data (Wave 1 and Wave 3 in-home surveys, RData format)
# Object name in both files is "x"; load into isolated environments to rename
# ---------------------------------------------------------------------------

e1 = new.env()
load("~/Downloads/w1inhome_dvn.RData", envir = e1)
addhealth_w1 = e1$x

e3 = new.env()
load("~/Downloads/w3inhome_dvn.RData", envir = e3)
addhealth_w3 = e3$x

rm(e1, e3)

glimpse(addhealth_w1)
glimpse(addhealth_w3)

# ---------------------------------------------------------------------------
# Extract religion variables and unique identifier
# W1: AID (respondent ID), H1RE1 (respondent religion), PA22 (parent religion)
# W3: aid (lowercase in this wave), h3re1 (respondent religion at Wave 3)
# ---------------------------------------------------------------------------

w1_relig = addhealth_w1 |>
  select(AID, H1RE1, PA22)

w3_relig = addhealth_w3 |>
  select(AID = aid, H3RE1 = h3re1, H3RE26 = h3re26)

glimpse(w1_relig)
glimpse(w3_relig)

# ---------------------------------------------------------------------------
# Join waves on AID (inner join — keep only respondents present in both waves)
# W1: 6,504 | W3: 4,882 | Joined: 4,882
# All W3 AIDs matched to W1; 1,622 W1 respondents dropped due to Wave 3 attrition
# ---------------------------------------------------------------------------

addhealth = inner_join(w1_relig, w3_relig, by = "AID")

cat("W1 rows:    ", nrow(w1_relig), "\n")
cat("W3 rows:    ", nrow(w3_relig), "\n")
cat("Joined rows:", nrow(addhealth), "\n")
cat("Dropped (W1 attrition):", nrow(w1_relig) - nrow(addhealth), "\n")

# ---------------------------------------------------------------------------
# Recode religion variables to a common 7-category scheme
#
# Categories:
#   0 = None / atheist / agnostic
#   1 = Protestant  (incl. JW, LDS, Christian Science, Unitarian, "just Christian")
#   2 = Catholic
#   3 = Jewish
#   4 = Buddhist
#   5 = Hindu
#   6 = Moslem
#   7 = Other       (Eastern Orthodox, Baha'i, Other religion)
#   NA = Refused / Don't know / Not applicable
# ---------------------------------------------------------------------------

# Variables are ordered factors with labels like "(22) Catholic"; extract the
# leading numeric code into temporary columns, recode, then drop temporaries

addhealth = addhealth |>
  mutate(
    # Extract numeric codes from factor labels
    .h1re1  = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H1RE1))),
    .pa22   = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(PA22))),
    .h3re1  = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H3RE1))),
    .h3re26 = as.integer(sub("^\\(([0-9]+)\\).*", "\\1", as.character(H3RE26))),

    # H1RE1: respondent's current religion at Wave 1
    # 0=None, 1–19=Protestant denominations, 20=Baha'i, 21=Buddhist, 22=Catholic,
    # 23=Eastern Orthodox, 24=Hindu, 25=Islam, 26=Jewish, 27=Unitarian, 28=Other religion
    relig_child_w1 = case_when(
      .h1re1 == 0                ~ 0L,  # None
      .h1re1 %in% c(1:19, 27)   ~ 1L,  # Protestant
      .h1re1 == 22               ~ 2L,  # Catholic
      .h1re1 == 26               ~ 3L,  # Jewish
      .h1re1 %in% c(20, 21, 23, 24, 25, 28) ~ 4L,  # Other
      TRUE                       ~ NA_integer_
    ),

    # PA22: parent's self-reported religion at Wave 1
    # 1–3=Protestant, 5=Baptist, 6=Buddhist, 7=Catholic, 8–10=Protestant,
    # 11=Eastern Orthodox, 12–13=Protestant, 14=Hindu, 15=Protestant, 16=Islam,
    # 17=JW, 18=Jewish, 19=LDS, 20–21=Protestant, 23=Other Protestant,
    # 24=Other religion, 25–27=Protestant, 28=None
    relig_parent_w1 = case_when(
      .pa22 == 28                                       ~ 0L,  # None
      .pa22 %in% c(1:3, 5, 8:10, 12:13, 15, 17, 19:21,
                   23, 25:27)                           ~ 1L,  # Protestant
      .pa22 == 7                                        ~ 2L,  # Catholic
      .pa22 == 18                                       ~ 3L,  # Jewish
      .pa22 %in% c(6, 11, 14, 16, 24)                   ~ 4L,  # Other
      TRUE                                              ~ NA_integer_
    ),

    # H3RE1: respondent's current religion at Wave 3
    # "just Christian" (8) collapsed into Protestant (1)
    relig_adult_w3 = case_when(
      .h3re1 == 0           ~ 0L,  # None
      .h3re1 %in% c(1, 8)  ~ 1L,  # Protestant
      .h3re1 == 2           ~ 2L,  # Catholic
      .h3re1 == 3           ~ 3L,  # Jewish
      .h3re1 %in% c(4, 5, 6, 7) ~ 4L,  # Other
      TRUE                  ~ NA_integer_
    ),

    # H3RE26: retrospective recall of childhood religion at Wave 3
    # Same scheme as H3RE1
    relig_recall_w3 = case_when(
      .h3re26 == 0          ~ 0L,  # None
      .h3re26 %in% c(1, 8)       ~ 1L,  # Protestant
      .h3re26 == 2                ~ 2L,  # Catholic
      .h3re26 == 3                ~ 3L,  # Jewish
      .h3re26 %in% c(4, 5, 6, 7) ~ 4L,  # Other
      TRUE                  ~ NA_integer_
    )
  ) |>
  select(-.h1re1, -.pa22, -.h3re1, -.h3re26)

relig_labs = c("None", "Protestant", "Catholic", "Jewish", "Other")

# ---------------------------------------------------------------------------
# PA22 missingness diagnostic
# Cross-tabulate child's own W1 religion (always observed) by whether
# parental religion (PA22) is missing — to assess whether missingness is
# random with respect to religious background
# ---------------------------------------------------------------------------

miss_tab = addhealth |>
  mutate(
    pa22_missing  = is.na(relig_parent_w1),
    child_relig   = factor(relig_child_w1, levels = 0:4, labels = relig_labs)
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
  filter(!is.na(relig_parent_w1), !is.na(relig_recall_w3)) |>
  mutate(
    parent  = factor(relig_parent_w1, levels = 0:4, labels = relig_labs),
    recall  = factor(relig_recall_w3, levels = 0:4, labels = relig_labs)
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
# Child–parent religious discordance at Wave 1
# How many respondents already differed from their parent's religion at age ~16?
# relig_child_w1 (H1RE1) vs relig_parent_w1 (PA22)
# This matters because PA22 is only a clean ground truth if the household
# shared a single religion — prior divergence inflates apparent recall mismatch
# ---------------------------------------------------------------------------

discord_dat = addhealth |>
  filter(!is.na(relig_child_w1), !is.na(relig_parent_w1)) |>
  mutate(
    child  = factor(relig_child_w1, levels = 0:4, labels = relig_labs),
    parent = factor(relig_parent_w1, levels = 0:4, labels = relig_labs),
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
  filter(!is.na(relig_parent_w1), !is.na(relig_recall_w3),
         !is.na(relig_adult_w3)) |>
  mutate(
    adult  = factor(relig_adult_w3,  levels = 0:4, labels = relig_labs),
    parent = factor(relig_parent_w1, levels = 0:4, labels = relig_labs),
    recall = factor(relig_recall_w3, levels = 0:4, labels = relig_labs),
    mismatch_type = case_when(
      recall == parent                          ~ "Concordant",
      recall == "None" & parent != "None"       ~ "Recalled more secular",
      parent == "None" & recall != "None"       ~ "Recalled more religious",
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
  filter(parent != "None") |>
  mutate(
    more_secular = as.integer(mismatch_type == "Recalled more secular"),
    adult_f  = relevel(droplevels(adult),  ref = "Protestant"),
    parent_f = relevel(droplevels(parent), ref = "Protestant")
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
  filter(!is.na(relig_child_w1)) |>
  mutate(child_f = relevel(
    factor(relig_child_w1, levels = 0:4, labels = relig_labs),
    ref = "Protestant"
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
    round(exp(coef(mod_base2)["adult_fNone"]), 3), "->",
    round(exp(coef(mod_child)["adult_fNone"]), 3), "\n")