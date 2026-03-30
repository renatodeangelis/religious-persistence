library(dplyr)
library(tidyr)
library(ggplot2)
library(haven)

data = read_dta(
  file = "https://www.dropbox.com/scl/fi/4tvctfek9wqskk7v34k57/gss7224_r3.dta?rlkey=dm9c0smlgbqneutvordpnr8yr&st=76655say&dl=1",
  col_select = c(year, cohort, relig, relig16, denom, denom16, other, oth16,
                 attend, reliten, relitenv, relitennv, sprel, spden, sprel16, spden16, educ))

data_filtered = data |>
  mutate(age = year - cohort,
         year = as.numeric(year),
         reliten = as.numeric(reliten),
         relitenv = as.numeric(relitenv),
         relitennv = as.numeric(relitennv),
         reliten_comb = coalesce(reliten, relitenv, relitennv))

## HOUT FIGURE 1

fig1_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16)) |>
  mutate(currently_catholic = as.numeric(relig) == 2,
         raised_catholic = as.numeric(relig16) == 2) |>
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
  
## HOUT FIGURE 2
fig2_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16),
         year %in% 1974:2024) |>
  mutate(group = case_when(
           relig16 != 2 & relig == 2 ~ "Converted to Catholic",
           relig16 == 2 & relig == 2 & reliten_comb == 1 ~ "Strong Catholic",
           relig16 == 2 & relig == 2 & reliten_comb %in% 2:3 ~ "Not strong Catholic",
           relig16 == 2 & relig == 4 ~ "No religion",
           relig16 == 2 & relig != 2 & reliten_comb == 1 ~ "Strong in new religion",
           relig16 == 2 & relig != 2 & reliten_comb %in% 2:3 ~ "Not strong in new religion")) |>
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
  "Converted to Catholic"       = "white",
  "Strong Catholic"             = "gray30",
  "Not strong Catholic"         = "gray65",
  "Strong in new religion"      = "gray40",
  "Not strong in new religion"  = "gray55",
  "No religion"                 = "gray10")

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
  mutate(cohort        = as.numeric(cohort),
         relig_macro   = classify_relig(relig,   denom,   other),
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

## FIGURE 5: PERSISTENCE BY EDUCATION

fig5_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(educ), !is.na(cohort),
         cohort %in% 1905:1995) |>
  mutate(cohort = as.numeric(cohort),
         educ_group = case_when(
           educ < 12       ~ "Less than HS",
           educ %in% 12:15 ~ "HS or some college",
           educ >= 16      ~ "College or more"
         ),
         stayed = as.numeric(relig) == as.numeric(relig16)) |>
  filter(!is.na(educ_group),
         as.numeric(relig16) != 4) |>
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

## FIGURE 6: CURRENTLY VS RAISED CATHOLIC BY EDUCATION AND YEAR

fig6_data = data_filtered |>
  filter(!is.na(relig), !is.na(relig16), !is.na(educ),
         year %in% 1974:2024) |>
  mutate(educ_group = case_when(
           educ < 12       ~ "Less than HS",
           educ %in% 12:15 ~ "HS or some college",
           educ >= 16      ~ "College or more"
         ),
         currently_catholic = as.numeric(relig)   == 2,
         raised_catholic    = as.numeric(relig16) == 2) |>
  filter(!is.na(educ_group)) |>
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
