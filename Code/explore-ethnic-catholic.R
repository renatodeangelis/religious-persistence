# ── EXPLORE: ETHNIC CATHOLICISM AND NOMINAL ORIGIN CLAIMING ───────────────────
# Throwaway script. Tests whether the quasi-ethnic Catholic hypothesis can
# explain low Catholic retention via a measurement artifact: Italian/Irish/Polish-
# identified respondents claim a Catholic upbringing at lower levels of actual
# childhood religiosity, inflating the Catholic origin denominator with nominal
# Catholics who were never going to retain.
#
# Data: attend12 waves only (1991, 1998, 2008, 2018) — same restriction as 16.
#       ethnic variable added.
#
# Three-step test:
#   Step 1 — exposure → mediator:  do ethnic identifiers report lower attend12?
#   Step 2 — total effect:         do ethnic identifiers retain less as adults?
#   Step 3 — mediation:            does the ethnic coefficient attenuate after
#                                  conditioning on attend12?

library(dplyr)
library(ggplot2)
library(gssr)
source("code/utils.R")

reltrad_labels = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

attend12_waves = c(1991, 1998, 2008, 2018)

# ── DATA ──────────────────────────────────────────────────────────────────────

data(gss_all)

d = gss_all |>
  filter(year %in% attend12_waves) |>
  select(year, cohort, reltrad, reltrad16, attend12, ethnic) |>
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
  haven::zap_labels() |>
  mutate(
    attend12_num = as.numeric(attend12),
    attend12_num = if_else(attend12_num == 10, NA_real_, attend12_num),
    attend12_bin = case_when(
      attend12_num %in% 6:8 ~ "regular",
      attend12_num %in% 0:5 ~ "non_regular",
      TRUE                  ~ NA_character_
    ),
    ethnic_num = as.numeric(ethnic),
    # Ireland=14, Italy=15, Poland=21
    ital_irish_polish = as.integer(ethnic_num %in% c(14L, 15L, 21L)),
    american_only     = as.integer(ethnic_num == 97L)
  )

# ── RESTRICT TO CATHOLIC ORIGIN ───────────────────────────────────────────────

cath = d |>
  filter(reltrad16_alt == "catholic", !is.na(ethnic_num)) |>
  mutate(retained = as.integer(reltrad_alt == "catholic"))

cat("── Catholic-origin sample (attend12 waves, with ethnic) ──\n")
cat("N total:", nrow(cath), "\n")
cat("N with attend12:", sum(!is.na(cath$attend12_num)), "\n")
cat("\nEthnic group breakdown:\n")
print(table(
  group = ifelse(cath$ital_irish_polish == 1, "Italian/Irish/Polish",
          ifelse(cath$american_only == 1,     "American only",
                                              "Other")),
  useNA = "ifany"
))

# ── STEP 1: EXPOSURE → MEDIATOR ───────────────────────────────────────────────
# Do ethnic identifiers report lower childhood attendance?

cat("\n\n── STEP 1: attend12 distribution by ethnic group ──\n")

attend_by_ethnic = cath |>
  filter(!is.na(attend12_num)) |>
  mutate(group = ifelse(ital_irish_polish == 1, "Italian/Irish/Polish", "Other Catholic")) |>
  group_by(group) |>
  summarise(
    n         = n(),
    mean_att  = round(mean(attend12_num), 3),
    pct_regular = round(mean(attend12_bin == "regular", na.rm = TRUE), 3),
    .groups = "drop"
  )
print(attend_by_ethnic)

# Cross-tab: attend12 raw scale × ethnic group
cat("\nRaw attend12 scale by group (row %)\n")
ctab = table(
  group   = ifelse(cath$ital_irish_polish == 1, "Ital/Irish/Pol", "Other Cath"),
  attend12 = cath$attend12_num
)
print(round(prop.table(ctab, margin = 1), 3))

# OLS: attend12 ~ ethnic + cohort (just to get a coefficient with SE)
m_att = lm(attend12_num ~ ital_irish_polish + cohort,
           data = filter(cath, !is.na(attend12_num)))
cat("\n── OLS: attend12_num ~ ital_irish_polish + cohort ──\n")
print(summary(m_att)$coefficients)

# ── STEP 2: TOTAL EFFECT ON RETENTION ─────────────────────────────────────────
# Do ethnic identifiers retain at lower rates (unconditional)?

cat("\n\n── STEP 2: retention by ethnic group (no attend12 control) ──\n")

ret_by_ethnic = cath |>
  mutate(group = ifelse(ital_irish_polish == 1, "Italian/Irish/Polish", "Other Catholic")) |>
  group_by(group) |>
  summarise(n = n(), pct_retained = round(mean(retained), 3), .groups = "drop")
print(ret_by_ethnic)

m_ret1 = glm(retained ~ ital_irish_polish + cohort,
             data = cath, family = binomial())
cat("\n── Logit: retained ~ ital_irish_polish + cohort ──\n")
print(round(summary(m_ret1)$coefficients, 4))

# ── STEP 3: MEDIATION — ADD ATTEND12 ─────────────────────────────────────────
# Does the ethnic coefficient shrink after conditioning on attend12?

m_ret2 = glm(retained ~ ital_irish_polish + attend12_num + cohort,
             data = filter(cath, !is.na(attend12_num)), family = binomial())
cat("\n\n── STEP 3: Logit: retained ~ ital_irish_polish + attend12_num + cohort ──\n")
print(round(summary(m_ret2)$coefficients, 4))

cat("\n── Coefficient on ital_irish_polish: without vs. with attend12 ──\n")
# Refit step 2 model on the same attend12-non-missing sample for comparability
m_ret1b = glm(retained ~ ital_irish_polish + cohort,
              data = filter(cath, !is.na(attend12_num)), family = binomial())
coef_without = coef(m_ret1b)["ital_irish_polish"]
coef_with    = coef(m_ret2)["ital_irish_polish"]
cat(sprintf("  Without attend12:  %.4f\n  With attend12:     %.4f\n  Attenuation:       %.1f%%\n",
            coef_without, coef_with,
            100 * (coef_without - coef_with) / coef_without))

# ── FIGURE: attend12 distribution by group ────────────────────────────────────

dir.create("output/figures/explore", recursive = TRUE, showWarnings = FALSE)

att_plot_df = cath |>
  filter(!is.na(attend12_num)) |>
  mutate(group = factor(
    ifelse(ital_irish_polish == 1, "Italian/Irish/Polish", "Other Catholic"),
    levels = c("Italian/Irish/Polish", "Other Catholic")
  )) |>
  group_by(group, attend12_num) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(group) |>
  mutate(pct = n / sum(n))

p_att = ggplot(att_plot_df, aes(x = factor(attend12_num), y = pct, fill = group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Italian/Irish/Polish" = "#5C4A8F",
                               "Other Catholic"       = "#0072B2"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(labels = c("0"="Never","1"="<1/yr","2"="1-2/yr","3"="Sev/yr",
                               "4"="1/mo","5"="2-3/mo","6"="~Wkly","7"="Wkly","8"="Sev/wk")) +
  labs(x = "Attendance at age 12", y = "Share of group",
       title = "Childhood Attendance Distribution by Ethnic Group",
       subtitle = "Catholic-origin respondents only (attend12 waves: 1991, 1998, 2008, 2018)") +
  healy_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("output/figures/explore/ethnic_attend12_dist.png",
       p_att, width = 9, height = 5, dpi = 200)

cat("\nDone. Figure: output/figures/explore/ethnic_attend12_dist.png\n")
