library(dplyr)
library(tidyr)
library(ggplot2)
library(survey)
library(haven)

data = read_dta(
  file = "https://www.dropbox.com/scl/fi/4tvctfek9wqskk7v34k57/gss7224_r3.dta?rlkey=dm9c0smlgbqneutvordpnr8yr&st=76655say&dl=1",
  col_select = c(year, cohort, relig, relig16, denom, denom16, other, oth16,
                 attend, reliten, sprel, spden, sprel16, spden16, educ))

data_filtered = data |>
  filter(!(year%in% 1972:1973)) |>
  mutate(wtssps = as.numeric(wtssps))

## HOUT FIGURE 1

fig1_data = data_filtered |>
  filter(!is.na(cohort), !is.na(relig), !is.na(relig16)) |>
  mutate(year = as.numeric(year),
         currently_catholic = as.numeric(relig) == 2,
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
  filter(!is.na(year), !is.na(relig16)) |>
  mutate(
    year = as.numeric(year),
    raised_catholic = as.numeric(relig16) == 2,
    currently_catholic = as.numeric(relig) == 2,
    strong = as.numeric(reliten) == 1,
    group = case_when(
      !raised_catholic ~ NA_character_,
      currently_catholic & strong ~ "Strong Catholic",
      currently_catholic & !strong ~ "Not strong Catholic",
      !currently_catholic ~ "Switched"
    )
  ) |>
  filter(!is.na(group), raised_catholic) |>
  group_by(year, group) |>
  summarise(n = sum(wtssps, na.rm = TRUE), .groups = "drop") |>
  group_by(year) |>
  mutate(pct = n / sum(n) * 100)


svy = svydesign(ids = ~1, weights = ~wtssps, data = data)
svymean(~as_factor(relig), design = svy, na.rm = TRUE)
