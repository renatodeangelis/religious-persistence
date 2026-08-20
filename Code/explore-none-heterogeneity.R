# ── EXPLORE · INTRA-"NONE" HETEROGENEITY ──────────────────────────────────────
# Exploratory analysis of heterogeneity within the religiously unaffiliated.
# Distinguishes cradle nones (raised without religion) from nonverts (adults
# who left a tradition), and examines how belief and practice vary by pathway
# and — for nonverts — by tradition of origin.
#
# Not destined for the paper as currently drafted; for coauthor circulation.
# Run manually — not sourced by 00-run-all.R.
#
# Input:  gss_all (gssr package)   — rebuilt here; not gss_clean.rds
# Output: console (cell counts, item-coverage rates)
#         output/figures/none-heterogeneity/*.png

if (FALSE) {

library(dplyr)
library(tidyr)
library(ggplot2)
library(gssr)
source("code/utils.R")

# ── LOCAL CONSTANTS ───────────────────────────────────────────────────────────
# Copied from script 16 for consistent tradition colors across the project.

reltrad_labels = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none"
)


pathway_colors = c("cradle none" = "#999999", "nonvert" = "#E69F00")

group_labels = c(
  "cradle none" = "Cradle none",
  "catholic"    = "Catholic",
  "evangelical" = "Evangelical",
  "mainline"    = "Mainline",
  "other"       = "Other"
)

# ── DATA ──────────────────────────────────────────────────────────────────────
# Rebuilt from gss_all rather than gss_clean.rds because the belief and
# practice variables (attend, god, pray) are not carried in the derived file.
# Applies the same sample restrictions as 01-prepare-data.R.

data(gss_all)

data = gss_all |>
  select(year, cohort, reltrad, reltrad16, attend, god, pray) |>
  filter(!(year %in% c(1972, 1982, 1987, 2021))) |>
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
  mutate(cohort   = as.numeric(cohort),
         age      = year - cohort,
         cohort_5 = (floor((cohort - 1925) / 5) * 5 + 1925) + 2.5) |>
  filter(age >= 30, age <= 75, cohort >= 1925, cohort <= 1984) |>
  haven::zap_labels()

# ── RECODES ───────────────────────────────────────────────────────────────────

data = data |>
  mutate(
    # Pathway into non-affiliation; defined only where current state is "none".
    none_type = case_when(
      reltrad_alt == "none" & reltrad16_alt == "none" ~ "cradle none",
      reltrad_alt == "none" & reltrad16_alt != "none" ~ "nonvert",
      TRUE                                            ~ NA_character_
    ),
    # Origin tradition among nonverts, for figures 3 and 4.
    nonvert_origin = if_else(none_type == "nonvert", reltrad16_alt, NA_character_),

    # attend: 0=never … 8=several times/wk (higher = more frequent)
    attends_yearly  = as.numeric(attend) >= 2,
    attends_monthly = as.numeric(attend) >= 4,

    # pray: 1=several times/day … 6=never (LOWER values = MORE frequent)
    prays_ever   = as.numeric(pray) <= 5,
    prays_weekly = as.numeric(pray) <= 4,

    # god: 1=don't believe … 6=know God exists
    god_cat = case_when(
      as.numeric(god) %in% 1:2 ~ "atheist/agnostic",
      as.numeric(god) == 3     ~ "higher power",
      as.numeric(god) %in% 4:6 ~ "believer",
      TRUE                      ~ NA_character_
    )
  )

# ── CONSOLE: CELL COUNTS AND ITEM COVERAGE ───────────────────────────────────

cat("\n── Overall N by survey year ──\n")
print(table(data$year))

cat("\n── N by none_type and survey year ──\n")
print(table(data$year, data$none_type, useNA = "ifany"))

cat("\n── Cradle-none N per year (thin-cell diagnostic) ──\n")
print(
  data |>
    filter(none_type == "cradle none") |>
    count(year),
  n = Inf
)

cat("\n── Non-missing rate for 'god' by year ──\n")
print(
  data |>
    group_by(year) |>
    summarise(n_total = n(),
              n_nonmiss = sum(!is.na(god_cat)),
              pct_nonmiss = round(100 * n_nonmiss / n_total, 1)),
  n = Inf
)

cat("\n── Non-missing rate for 'pray' by year ──\n")
print(
  data |>
    group_by(year) |>
    summarise(n_total = n(),
              n_nonmiss = sum(!is.na(prays_ever)),
              pct_nonmiss = round(100 * n_nonmiss / n_total, 1)),
  n = Inf
)

# ── OUTPUT DIRECTORY ─────────────────────────────────────────────────────────

dir.create("output/figures/none-heterogeneity", recursive = TRUE,
           showWarnings = FALSE)

# ── FIGURE 1: CRADLE-NONE VS. NONVERT SHARE BY BIRTH COHORT ─────────────────
# Among currently unaffiliated respondents, proportion who are cradle nones vs.
# nonverts, by 5-year birth cohort bin.

fig1_data = data |>
  filter(!is.na(none_type)) |>
  group_by(cohort_5, none_type) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(cohort_5) |>
  mutate(share = n / sum(n)) |>
  ungroup() |>
  mutate(se   = sqrt(share * (1 - share) / n),
         ymin = pmax(0, share - 1.96 * se),
         ymax = pmin(1, share + 1.96 * se))

cat("\n── Figure 1: none_type share by birth cohort ──\n")
print(fig1_data, n = Inf)

fig1 = ggplot(fig1_data,
              aes(x = cohort_5, y = share, color = none_type, group = none_type)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = none_type),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = pathway_colors,
                     labels = c("cradle none" = "Cradle none",
                                "nonvert"     = "Nonvert")) +
  scale_fill_manual(values = pathway_colors, guide = "none") +
  labs(x = "Birth cohort (5-year midpoint)",
       y = "Share of all nones",
       color = NULL,
       title = "Composition of the Unaffiliated by Pathway into Non-Affiliation",
       caption = "Ages 30–75, cohorts 1925–1984. Unweighted.") +
  healy_theme

ggsave("output/figures/none-heterogeneity/fig1_cradle_share.png",
       fig1, width = 9, height = 5, dpi = 200)

# ── FIGURE 2: ATHEIST/AGNOSTIC SHARE BY PATHWAY, 1988+ ───────────────────────
# Share of nones in the "atheist/agnostic" god_cat category, by survey year
# and pathway. Truncated at 1988 (first god wave). Pre-2006 points are open
# with no connecting line to flag sporadic item coverage.

fig2_data = data |>
  filter(!is.na(none_type), !is.na(god_cat), year >= 1988) |>
  group_by(year, none_type) |>
  summarise(n = n(),
            share_ath = mean(god_cat == "atheist/agnostic"),
            .groups = "drop") |>
  mutate(era  = if_else(year < 2006, "pre2006", "post2006"),
         se   = sqrt(share_ath * (1 - share_ath) / n),
         ymin = pmax(0, share_ath - 1.96 * se),
         ymax = pmin(1, share_ath + 1.96 * se))

cat("\n── Figure 2: atheist/agnostic share among nones by year and pathway ──\n")
print(fig2_data, n = Inf)

fig2 = ggplot(fig2_data,
              aes(x = year, y = share_ath, color = none_type)) +
  geom_ribbon(data = filter(fig2_data, era == "post2006"),
              aes(ymin = ymin, ymax = ymax, fill = none_type, group = none_type),
              alpha = 0.2, color = NA) +
  geom_line(data = filter(fig2_data, era == "post2006"),
            aes(group = none_type), linewidth = 0.8) +
  geom_point(data = filter(fig2_data, era == "post2006"),
             size = 2.5, shape = 16) +
  geom_point(data = filter(fig2_data, era == "pre2006"),
             size = 2.5, shape = 1) +
  geom_vline(xintercept = 2006, linetype = "dashed", color = "grey55",
             linewidth = 0.5) +
  annotate("text", x = 2007, y = 0.04,
           label = "2006: consistent\ncoverage begins",
           hjust = 0, size = 3, color = "grey45") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = pathway_colors,
                     labels = c("cradle none" = "Cradle none",
                                "nonvert"     = "Nonvert")) +
  scale_fill_manual(values = pathway_colors, guide = "none") +
  labs(x = "Survey year",
       y = "Share atheist or agnostic (god)",
       color = NULL,
       title = "Atheist/Agnostic Belief among the Unaffiliated, by Pathway",
       caption = "Open points (pre-2006): sporadic god item coverage. Ages 30–75, cohorts 1925–1984. Unweighted.") +
  healy_theme

ggsave("output/figures/none-heterogeneity/fig2_god_by_pathway.png",
       fig2, width = 9, height = 5, dpi = 200)

# ── FIGURE 3: NONVERT ORIGIN COMPOSITION BY BIRTH COHORT ────────────────────
# Among nonverts only, composition by origin tradition across 5-year birth
# cohort bins. Reveals which traditions fed non-affiliation for each cohort.

origin_order = c("catholic", "evangelical", "mainline", "other")

fig3_data = data |>
  filter(none_type == "nonvert", !is.na(nonvert_origin)) |>
  mutate(nonvert_origin = factor(nonvert_origin, levels = origin_order)) |>
  group_by(cohort_5, nonvert_origin) |>
  summarise(n = n(), .groups = "drop")

cat("\n── Figure 3: nonvert origin counts by birth cohort ──\n")
print(fig3_data, n = Inf)

fig3 = ggplot(fig3_data,
              aes(x = cohort_5, y = n, fill = nonvert_origin, group = nonvert_origin)) +
  geom_area(position = "fill", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = reltrad_colors[origin_order],
                    labels = reltrad_labels_tc[origin_order]) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Birth cohort (5-year midpoint)",
       y = "Share of nonverts",
       fill = "Origin tradition",
       title = "Origin Tradition among Nonverts, by Birth Cohort",
       caption = "Ages 30–75, cohorts 1925–1984. Unweighted.") +
  healy_theme

ggsave("output/figures/none-heterogeneity/fig3_nonvert_origins.png",
       fig3, width = 9, height = 5, dpi = 200)

# ── FIGURE 4: RELIGIOUS RESIDUE GRADIENT ─────────────────────────────────────
# Pooled across all years. Two-facet grouped bar chart: share attends_yearly
# and share prays_ever, by group (cradle none + four nonvert origins). Groups
# sorted by attends_yearly share (descending) on x-axis.

fig4_summary = data |>
  filter(!is.na(none_type)) |>
  mutate(group = if_else(none_type == "cradle none", "cradle none",
                         nonvert_origin)) |>
  filter(!is.na(group)) |>
  group_by(group) |>
  summarise(
    att_share  = mean(attends_yearly, na.rm = TRUE),
    att_n      = sum(!is.na(attends_yearly)),
    pray_share = mean(prays_ever,     na.rm = TRUE),
    pray_n     = sum(!is.na(prays_ever)),
    n          = n(),
    .groups = "drop"
  )

group_order = fig4_summary |>
  arrange(desc(att_share)) |>
  pull(group)

fig4_long = fig4_summary |>
  pivot_longer(cols = c(att_share, pray_share),
               names_to = "measure", values_to = "share") |>
  mutate(
    n_measure = if_else(measure == "att_share", att_n, pray_n),
    se        = sqrt(share * (1 - share) / n_measure),
    ymin      = pmax(0, share - 1.96 * se),
    ymax      = pmin(1, share + 1.96 * se),
    group     = factor(group, levels = group_order),
    measure   = factor(measure,
                       levels = c("att_share", "pray_share"),
                       labels = c("Attends at least yearly",
                                  "Prays at least occasionally"))
  )

cat("\n── Figure 4: pooled practice rates by group ──\n")
print(select(fig4_summary, group, att_share, att_n, pray_share, pray_n, n))

fig4 = ggplot(fig4_long, aes(x = group, y = share, fill = group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey30") +
  geom_text(aes(y = ymax + 0.03, label = sprintf("n=%d", n_measure)),
            vjust = 0, size = 2.8) +
  facet_wrap(~ measure) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(reltrad_colors, "cradle none" = "#999999"),
                    guide = "none") +
  scale_x_discrete(labels = group_labels) +
  labs(x = NULL,
       y = "Share",
       title = "Religious Residue by Pathway and Origin Tradition",
       caption = "Error bars: 95% CI. Groups sorted by attendance share (descending). Pooled across all years. Ages 30–75, cohorts 1925–1984. Unweighted.") +
  healy_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/none-heterogeneity/fig4_residue_gradient.png",
       fig4, width = 10, height = 5, dpi = 200)

# ── FIGURE 5: ATTENDANCE BY AGE AND PATHWAY ──────────────────────────────────
# Share reporting attends_yearly by 5-year age bin, one line each for cradle
# nones and nonverts. Cross-sectional — cohort and period effects are
# confounded.

fig5_data = data |>
  filter(!is.na(none_type), !is.na(attends_yearly)) |>
  mutate(age_bin = cut(age, breaks = seq(30, 75, by = 5), right = FALSE,
                       labels = seq(30, 70, by = 5))) |>
  filter(!is.na(age_bin)) |>
  group_by(age_bin, none_type) |>
  summarise(n = n(),
            share = mean(attends_yearly),
            .groups = "drop") |>
  mutate(age_bin = as.numeric(as.character(age_bin)),
         se      = sqrt(share * (1 - share) / n),
         ymin    = pmax(0, share - 1.96 * se),
         ymax    = pmin(1, share + 1.96 * se))

cat("\n── Figure 5: attends_yearly by age bin and pathway ──\n")
print(fig5_data, n = Inf)

fig5 = ggplot(fig5_data,
              aes(x = age_bin, y = share,
                  color = none_type, group = none_type)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = none_type),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(30, 70, by = 5)) +
  scale_color_manual(values = pathway_colors,
                     labels = c("cradle none" = "Cradle none",
                                "nonvert"     = "Nonvert")) +
  scale_fill_manual(values = pathway_colors, guide = "none") +
  labs(x = "Age at interview (5-year bins)",
       y = "Share attends at least yearly",
       color = NULL,
       title = "Attendance by Age and Pathway into Non-Affiliation",
       caption = "Cross-sectional — cohort and period effects are confounded. Ages 30–75, cohorts 1925–1984. Unweighted.") +
  healy_theme

ggsave("output/figures/none-heterogeneity/fig5_age_profile.png",
       fig5, width = 9, height = 5, dpi = 200)

# ── NOTES ─────────────────────────────────────────────────────────────────────
# 1. Cradle nones are very rare among early cohorts. The left edge of figures
#    1 and 3 (cohorts born ~1925–1944) is unstable; see cradle-none N per
#    cohort_5 bin in the fig1_data console output above.
# 2. `god` coverage is sporadic before 2006 and concentrated in later waves.
#    Any pooled comparison involving god_cat weights later years heavily (a
#    compositional artifact). Figure 4 uses attend and pray instead.
# 3. `god`, `pray`, and `postlife` were on rotating ballots in some years;
#    see item-coverage rates in console output above.
# 4. No origin-side intensity measure exists in the GSS. These figures
#    characterize the destination none state only; they cannot speak to whether
#    the origin none state has changed in character. (NORC Methodological
#    Report #146 notes the same constraint in explaining why RELTRAD16 is
#    constructed without an attendance component while RELTRAD includes one.)

cat("\nWrote output/figures/none-heterogeneity/\n")

} # end if (FALSE)
