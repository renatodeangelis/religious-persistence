library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gssr)
source("Code/utils.R")

data(gss_all)
data = gss_all |>
  select(year, cohort, relig, relig16, denom, denom16, other, oth16,
         attend, reliten, relitenv, relitennv, sprel, spden, sprel16, spden16, educ,
         race, hispanic, born, parborn) |>
  filter(!(year %in% c(1972, 2021)))

data_filtered = data |>
  mutate(age = year - cohort,
         year   = as.numeric(year),
         cohort = as.numeric(cohort),
         relig   = as.numeric(relig),
         relig16 = as.numeric(relig16),
         reliten    = as.numeric(reliten),
         relitenv   = as.numeric(relitenv),
         relitennv  = as.numeric(relitennv),
         reliten_comb = coalesce(reliten, relitenv, relitennv),
         cohort_5 = floor((cohort - 1900) / 5) * 5 + 1900,
         educ_group = case_when(
           educ < 12       ~ "Less than HS",
           educ %in% 12:15 ~ "HS or some college",
           educ >= 16      ~ "College or more"
         ))

dir.create("output/figures/hout",     recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/catholic", recursive = TRUE, showWarnings = FALSE)

## HOUT FIGURE 1

fig1_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16)) |>
  mutate(currently_catholic = relig == 2,
         raised_catholic    = relig16 == 2) |>
  group_by(year) |>
  summarise(
    pct_current = mean(currently_catholic, na.rm = TRUE) * 100,
    pct_raised = mean(raised_catholic, na.rm = TRUE) * 100) |>
  pivot_longer(cols = c(pct_current, pct_raised),
               names_to = "series", values_to = "pct") |>
  mutate(series = recode(series,
                         pct_current = "Currently Catholic",
                         pct_raised = "Raised Catholic")) 


ggplot(fig1_data, aes(x = year, y = pct, color = series)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 10)) +
  scale_color_manual(values = c("Currently Catholic" = "gray50",
                                "Raised Catholic" = "black")) +
  labs(x = "Year", y = "Catholic (%)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/figures/hout/fig1.png", width = 7, height = 5)

## HOUT FIGURE 2
fig2_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16),
         year %in% 1974:2024) |>
  mutate(group = case_when(
           relig16 != 2 & relig == 2 ~ "Converted to Catholic",
           relig16 == 2 & relig == 2 & reliten_comb == 1 ~ "Strong Catholic",
           relig16 == 2 & relig == 2 & (reliten_comb %in% 2:3 | is.na(reliten_comb)) ~ "Not strong Catholic",
           relig16 == 2 & relig == 4 ~ "No religion",
           relig16 == 2 & relig != 2 & reliten_comb == 1 ~ "Strong in new religion",
           relig16 == 2 & relig != 2 & (reliten_comb %in% 2:3 | is.na(reliten_comb)) ~ "Not strong in new religion")) |>
  count(year, group) |>
  group_by(year) |>
  mutate(pct_raw = n / sum(n) * 100,
         pct_raw = if_else(group %in% c("No religion", "Strong in new religion",
                                        "Not strong in new religion"), -pct_raw, pct_raw)) |>
  ungroup() |>
  filter(!is.na(group)) |>
  reframe(year = year,
          pct_smooth = predict(loess(pct_raw ~ year, span = 0.4)),
          .by = group)
    
group_levels = c(
  "Converted to Catholic",
  "Strong Catholic",
  "Not strong Catholic",
  "No religion",
  "Not strong in new religion",
  "Strong in new religion")

group_colors = c(
  "Converted to Catholic"       = "#FFFFFF",  # white
  "Strong Catholic"             = "#0072B2",  # Okabe-Ito blue
  "Not strong Catholic"         = "#56B4E9",  # Okabe-Ito sky
  "Strong in new religion"      = "#D55E00",  # Okabe-Ito vermillion
  "Not strong in new religion"  = "#E69F00",  # Okabe-Ito orange
  "No religion"                 = "#555555")  # dark gray

band_labels = tibble::tribble(
  ~group,                      ~x,     ~y,    ~color,
  "Converted to Catholic",      1999,   23, "black",
  "Strong Catholic",            1999,   18, "white",
  "Not strong Catholic",        1999,    6.5, "black",
  "Strong in new religion",     1999,   -1, "white",
  "Not strong in new religion", 1999,   -3, "black",
  "No religion",                1999,   -5.5, "white")

fig2_data |>
  mutate(group = factor(group, levels = group_levels)) |>
  ggplot(aes(x = year, y = pct_smooth, fill = group)) +
  geom_area(position = "stack", color = "black", linewidth = 0.3, outline.type = "full") +
  geom_hline(yintercept = 0, linewidth = 0) +
  geom_text(data = band_labels, aes(x = x, y = y, label = group, color = color),
            inherit.aes = FALSE, size = 3) +
  scale_color_identity() +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(breaks = seq(-10, 30, 10),
                     sec.axis = dup_axis(name = NULL)) +
  annotate("text", x = 1971, y = 15,  label = "Currently →",
           angle = 90, size = 3) +
  annotate("text", x = 1971, y = -5, label = "← Formerly",
           angle = 90, size = 3) +
  labs(x = "Year", y = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y.right = element_text(color = "gray50"))

ggsave("output/figures/hout/fig2.png", width = 7, height = 5)

## HOUT FIGURE 3
fig3_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(cohort),
         year %in% 1974:2024,
         cohort %in% 1905:1995) |>
  mutate(group = case_when(
    relig16 != 2 & relig == 2 ~ "Converted to Catholic",
    relig16 == 2 & relig == 2 & reliten_comb == 1 ~ "Strong Catholic",
    relig16 == 2 & relig == 2 & (reliten_comb %in% 2:3 | is.na(reliten_comb)) ~ "Not strong Catholic",
    relig16 == 2 & relig == 4 ~ "No religion",
    relig16 == 2 & relig != 2 & reliten_comb == 1 ~ "Strong in new religion",
    relig16 == 2 & relig != 2 & (reliten_comb %in% 2:3 | is.na(reliten_comb)) ~ "Not strong in new religion")) |>
  count(cohort, group) |>
  group_by(cohort) |>
  mutate(pct_raw = n / sum(n) * 100,
         pct_raw = if_else(group %in% c("No religion", "Strong in new religion",
                                        "Not strong in new religion"), -pct_raw, pct_raw)) |>
  ungroup() |>
  filter(!is.na(group)) |>
  reframe(cohort = cohort,
          pct_smooth = predict(loess(pct_raw ~ cohort, span = 0.4)),
          .by = group)

band_labels = tibble::tribble(
  ~group,                      ~x,     ~y,    ~color,
  "Converted to Catholic",      1960,   24, "black",
  "Strong Catholic",            1950,   18, "white",
  "Not strong Catholic",        1960,    6.5, "black",
  "Strong in new religion",     1955,   -1, "white",
  "Not strong in new religion", 1970,   -3.5, "black",
  "No religion",                1975,   -7, "white")

fig3_data |>
  mutate(group = factor(group, levels = group_levels)) |>
  ggplot(aes(x = cohort, y = pct_smooth, fill = group)) +
  geom_area(position = "stack", color = "black", linewidth = 0.3, outline.type = "full") +
  geom_hline(yintercept = 0, linewidth = 0) +
  geom_text(data = band_labels, aes(x = x, y = y, label = group, color = color),
            inherit.aes = FALSE, size = 3) +
  scale_color_identity() +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(breaks = seq(-10, 30, 10),
                     sec.axis = dup_axis(name = NULL)) +
  annotate("text", x = 1900, y = 15,  label = "Currently →",
           angle = 90, size = 3) +
  annotate("text", x = 1900, y = -5, label = "← Formerly",
           angle = 90, size = 3) +
  labs(x = "Birth Year", y = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y.right = element_text(color = "gray50"))

ggsave("output/figures/hout/fig3.png", width = 7, height = 5)

## HOUT FIGURE 4

cons_other = c(
  68, 64, 35, 31, 58, 12, 36, 77, 65, 56, 45, 32, 69, 71, 16, 51, 83, 63, 24,
  97, 52, 47, 9, 53, 18, 23, 3, 43, 66, 76, 67, 41, 13, 102, 117, 42, 138,
  103, 155, 27, 107, 181, 91, 28, 84, 100, 111, 113, 125, 22, 57, 92, 109,
  112, 2, 6, 26, 115, 154, 182, 120, 121, 118, 176, 208, 197)

ml_other = c(
  40, 81, 80, 44, 70, 20, 54, 46, 48, 25, 72, 73, 89, 99, 96, 119, 150, 153,
  105, 188)

classify_relig = function(relig, denom, other) {
  relig = as.numeric(relig)
  denom = as.numeric(denom)
  other = as.numeric(other)
  case_when(
    relig == 2 ~ "Catholic",
    relig == 13 | denom %in% c(2, 3, 4, 5, 7) | other %in% ml_other ~ "Mainline",
    relig == 11 | denom == 1 | other %in% cons_other ~ "Conservative",
    relig == 4 ~ "No religion",
    relig == 1 ~ "Mainline",
    TRUE ~ "Other"
  )
}

fig4_base = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), cohort %in% 1905:1995,
         year %in% 1974:2024) |>
  mutate(relig_macro   = classify_relig(relig,   denom,   other),
         relig16_macro = classify_relig(relig16, denom16, oth16))

fig4_cohort_n = count(fig4_base, cohort, name = "total_n")

# Above zero: panel = current religion, includes converts in
fig4_above = fig4_base |>
  filter(relig_macro %in% c("Catholic", "Mainline", "Conservative")) |>
  mutate(panel = relig_macro,
         group = case_when(
           relig16_macro != relig_macro ~ "Converted",
           reliten_comb == 1            ~ "Strong",
           TRUE                         ~ "Not strong"
         ))

# Below zero: panel = origin religion, switchers only
fig4_below = fig4_base |>
  filter(relig16_macro %in% c("Catholic", "Mainline", "Conservative"),
         relig_macro != relig16_macro) |>
  mutate(panel = relig16_macro,
         group = case_when(
           relig_macro == "No religion" ~ "No religion",
           reliten_comb == 1            ~ "Strong in new religion",
           TRUE                         ~ "Not strong in new religion"
         ))

fig4_data = bind_rows(fig4_above, fig4_below) |>
  count(cohort, panel, group) |>
  left_join(fig4_cohort_n, by = "cohort") |>
  mutate(pct_raw = n / total_n * 100,
         pct_raw = if_else(group %in% c("No religion", "Strong in new religion",
                                        "Not strong in new religion"),
                           -pct_raw, pct_raw)) |>
  reframe(cohort     = cohort,
          panel      = panel,
          pct_smooth = predict(loess(pct_raw ~ cohort, span = 0.4)),
          .by = c(panel, group))

fig4_group_levels = c("Converted", "Strong", "Not strong",
                      "No religion", "Not strong in new religion", "Strong in new religion")

fig4_group_colors = c(
  "Converted"                  = "#FFFFFF",  # white            — converts in
  "Strong"                     = "#0072B2",  # Okabe-Ito blue   — persistent, strong
  "Not strong"                 = "#56B4E9",  # Okabe-Ito sky    — persistent, weak
  "No religion"                = "#555555",  # dark gray        — left for none
  "Strong in new religion"     = "#D55E00",  # Okabe-Ito verm.  — switcher, strong
  "Not strong in new religion" = "#E69F00")  # Okabe-Ito orange — switcher, weak

fig4_panel_labels = c(
  Catholic     = "Catholic",
  Conservative = "Conservative Protestant",
  Mainline     = "Mainline Protestant")

fig4_data |>
  mutate(group = factor(group, levels = fig4_group_levels),
         panel = factor(panel, levels = c("Catholic", "Conservative", "Mainline"))) |>
  ggplot(aes(x = cohort, y = pct_smooth, fill = group)) +
  geom_area(position = "stack", color = "black", linewidth = 0.3, outline.type = "full") +
  geom_hline(yintercept = 0, linewidth = 0) +
  scale_fill_manual(values = fig4_group_colors) +
  scale_y_continuous(breaks = seq(-15, 45, 15),
                     sec.axis = dup_axis(name = NULL)) +
  facet_wrap(~ panel, nrow = 1,
             labeller = labeller(panel = fig4_panel_labels)) +
  labs(x = "Year of Birth", y = "Share of total population (%)", fill = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y.right = element_text(color = "gray50"),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("output/figures/hout/fig4.png", width = 10, height = 5)

## FIGURE 5: PERSISTENCE BY EDUCATION

fig5_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(educ_group), !is.na(cohort),
         cohort %in% 1905:1995,
         relig16 != 4) |>
  mutate(stayed = relig == relig16) |>
  group_by(cohort, educ_group) |>
  summarise(pct_stayed = mean(stayed, na.rm = TRUE) * 100, .groups = "drop")

fig5_data |>
  mutate(educ_group = factor(educ_group,
                             levels = c("Less than HS", "HS or some college", "College or more"))) |>
  ggplot(aes(x = cohort, y = pct_stayed)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE, color = "black") +
  facet_wrap(~ educ_group, nrow = 1) +
  labs(x = "Birth Year", y = "Stayed in origin religion (%)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/hout/fig5.png", width = 10, height = 5)

## FIGURE 6: CURRENTLY VS RAISED CATHOLIC BY EDUCATION AND YEAR

fig6_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(educ_group),
         year %in% 1974:2024) |>
  mutate(currently_catholic = relig   == 2,
         raised_catholic    = relig16 == 2) |>
  group_by(year, educ_group) |>
  summarise(pct_current = mean(currently_catholic, na.rm = TRUE) * 100,
            pct_raised  = mean(raised_catholic,    na.rm = TRUE) * 100,
            .groups = "drop") |>
  pivot_longer(cols = c(pct_current, pct_raised),
               names_to = "series", values_to = "pct") |>
  mutate(series = recode(series,
                         pct_current = "Currently Catholic",
                         pct_raised  = "Raised Catholic"))

fig6_data |>
  mutate(educ_group = factor(educ_group,
                             levels = c("Less than HS", "HS or some college", "College or more"))) |>
  ggplot(aes(x = year, y = pct, color = series)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE) +
  scale_color_manual(values = c("Currently Catholic" = "gray50",
                                "Raised Catholic"    = "black")) +
  facet_wrap(~ educ_group, nrow = 1) +
  labs(x = "Year", y = "Catholic (%)", color = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("output/figures/hout/fig6.png", width = 10, height = 5)

## FIGURE 7: RAISED VS CURRENTLY CATHOLIC BY NATIVITY (5-YEAR BIRTH COHORT)

fig7_base = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(cohort_5),
         year >= 1977,
         cohort >= 1930, cohort <= 1985,
         as.numeric(born) %in% c(1, 2)) |>
  mutate(
    raised_catholic  = relig16 == 2,
    current_catholic = relig   == 2,
    nativity = factor(
      if_else(as.numeric(born) == 1, "Born in US", "Born abroad"),
      levels = c("Born in US", "Born abroad")
    )
  )

fig7_catholic = fig7_base |>
  group_by(cohort_5, nativity) |>
  summarise(
    n_obs       = n(),
    pct_raised  = mean(raised_catholic,  na.rm = TRUE),
    pct_current = mean(current_catholic, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  pivot_longer(cols = c(pct_raised, pct_current),
               names_to = "series", values_to = "pct") |>
  mutate(
    se     = sqrt(pct * (1 - pct) / n_obs),
    ymin   = pct - 1.96 * se,
    ymax   = pct + 1.96 * se,
    series = recode(series, pct_raised = "Raised Catholic", pct_current = "Currently Catholic"),
    series = factor(series, levels = c("Raised Catholic", "Currently Catholic"))
  )

# Nativity composition strip: % born abroad by 5-year cohort
fig7_abroad = fig7_base |>
  group_by(cohort_5) |>
  summarise(n_obs = n(), pct = mean(nativity == "Born abroad"), .groups = "drop") |>
  mutate(
    se   = sqrt(pct * (1 - pct) / n_obs),
    ymin = pct - 1.96 * se,
    ymax = pct + 1.96 * se
  )

# Gray = raised (baseline); blue = currently (outcome) — both Okabe-Ito compliant
fig7_colors = c("Raised Catholic" = "#999999", "Currently Catholic" = "#0072B2")

p_fig7_main = ggplot() +
  geom_ribbon(data = fig7_catholic,
              aes(x = cohort_5, ymin = ymin, ymax = ymax, fill = series, group = series),
              alpha = 0.2, color = NA) +
  geom_line(data  = fig7_catholic,
            aes(x = cohort_5, y = pct, color = series, group = series),
            linewidth = 0.8) +
  geom_point(data = fig7_catholic,
             aes(x = cohort_5, y = pct, color = series, group = series),
             shape = 16, size = 2) +
  facet_wrap(~ nativity) +
  scale_x_continuous(breaks = seq(1930, 1980, by = 10),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(name = "Proportion Catholic") +
  scale_color_manual(values = fig7_colors, name = NULL) +
  scale_fill_manual(values  = fig7_colors, name = NULL, guide = "none") +
  labs(x = NULL) +
  healy_theme

# Strip spans full width below both panels — no spacers, so referent is unambiguous
p_fig7_strip = ggplot(fig7_abroad, aes(x = cohort_5, y = pct)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, fill = "#555555", color = NA) +
  geom_line(color = "#555555", linewidth = 0.8) +
  geom_point(color = "#555555", shape = 16, size = 2) +
  scale_x_continuous(breaks = seq(1930, 1980, by = 10),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(limits = c(0, 0.22),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Birth cohort (5-year bin)", y = "Share\nborn abroad") +
  healy_theme

patchwork::wrap_plots(p_fig7_main, p_fig7_strip, ncol = 1, heights = c(5, 1.5)) +
  patchwork::plot_annotation(
    caption = "Source: GSS, 1977–present. 5-year birth cohort bins, 1930–1985. Ribbons are 95% CIs.",
    theme   = theme(plot.caption = element_text(size = 9, hjust = 0))
  )

ggsave("output/figures/catholic/fig7.png", width = 9, height = 6.5, dpi = 200)

## 5-STATE TRANSITION MATRICES BY BIRTH COHORT AND NATIVITY

state_levels = c("Catholic", "Mainline", "Conservative", "No religion", "Other")

mat_base = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(cohort),
         as.numeric(born) %in% c(1, 2),
         cohort >= 1950, cohort < 1980) |>
  mutate(
    origin  = factor(classify_relig(relig16, denom16, oth16), levels = state_levels),
    dest    = factor(classify_relig(relig,   denom,   other),  levels = state_levels),
    cohort_10 = case_when(
      cohort < 1960 ~ "1950-1959",
      cohort < 1970 ~ "1960-1969",
      cohort < 1980 ~ "1970-1979"
    ),
    nativity = if_else(as.numeric(born) == 1, "Born in US", "Born abroad")
  ) |>
  filter(!is.na(origin), !is.na(dest))

cohort_windows = c("1950-1959", "1960-1969", "1970-1979")
nativity_groups = c("Born in US", "Born abroad")

for (coh in cohort_windows) {
  for (nat in nativity_groups) {
    sub = filter(mat_base, cohort_10 == coh, nativity == nat)
    cat("\n=== Cohort:", coh, "|", nat, "(N =", nrow(sub), ") ===\n")
    tbl = table(Origin = sub$origin, Destination = sub$dest)
    print(round(prop.table(tbl, margin = 1), 3))
  }
}
