library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gssr)
source("Code/utils.R")

reltrad_labels  = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")

reltrad_colors = c(
  catholic    = "#0072B2",
  evangelical = "#D55E00",
  mainline    = "#009E73",
  other       = "#CC79A7",
  none        = "#999999"
)
reltrad_labels_tc = c(
  catholic    = "Catholic",
  evangelical = "Evangelical",
  mainline    = "Mainline",
  other       = "Other",
  none        = "None"
)

# GSS decade windows (survey year, not birth cohort)
# Excluded years: 1972 and 2021 (standard); 1982, 1987, 2022, 2024 (oversamples)
oversample_years = c(1972, 1982, 1987, 2021, 2022, 2024)

decade_windows = c(
  "1973-1981", "1983-1992", "1993-2002", "2003-2012", "2013-2020"
)

# ── DATA ──────────────────────────────────────────────────────────────────────

data(gss_all)
data_gd = gss_all |>
  select(year, cohort, reltrad, reltrad16) |>
  filter(!(year %in% oversample_years)) |>
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
  mutate(age = year - cohort) |>
  filter(age > 25, cohort >= 1900) |>
  mutate(
    gss_decade = case_when(
      year >= 1973 & year <= 1981 ~ "1973-1981",
      year >= 1983 & year <= 1992 ~ "1983-1992",
      year >= 1993 & year <= 2002 ~ "1993-2002",
      year >= 2003 & year <= 2012 ~ "2003-2012",
      year >= 2013 & year <= 2020 ~ "2013-2020",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(gss_decade))

states_alt = sort(unique(c(data_gd$reltrad_alt, data_gd$reltrad16_alt)))
states_alt = states_alt[!is.na(states_alt)]

# Confirm which survey years are present in each window (sanity check)
cat("\n── Survey years by GSS decade window ──\n")
year_check = data_gd |>
  distinct(year, gss_decade) |>
  arrange(gss_decade, year) |>
  group_by(gss_decade) |>
  summarise(years = paste(year, collapse = ", "), .groups = "drop")
print(year_check, n = Inf)

# ── CELL-COUNT DIAGNOSTIC ─────────────────────────────────────────────────────

cell_diag = lapply(decade_windows, function(dw) {
  sub = data_gd[!is.na(data_gd$gss_decade) & data_gd$gss_decade == dw &
                !is.na(data_gd$reltrad16_alt) & !is.na(data_gd$reltrad_alt), ]
  if (nrow(sub) == 0) return(NULL)
  tab = table(
    factor(sub$reltrad16_alt, levels = states_alt),
    factor(sub$reltrad_alt,   levels = states_alt)
  )
  data.frame(
    gss_decade = dw,
    n          = nrow(sub),
    min_cell   = min(tab),
    cells_lt5  = sum(tab < 5),
    cells_0    = sum(tab == 0),
    row.names  = NULL
  )
}) |> do.call(what = rbind)

cat("\n── Cell-count diagnostic (5-state, GSS decade windows) ──\n")
print(cell_diag, row.names = FALSE)

# ── TRANSITION MATRICES BY GSS DECADE ────────────────────────────────────────

P_list_gd      = list()
pi0_list_gd    = list()
pistar_list_gd = list()
n_list_gd      = list()

for (dw in decade_windows) {
  sub = data_gd[!is.na(data_gd$gss_decade) & data_gd$gss_decade == dw &
                !is.na(data_gd$reltrad16_alt) & !is.na(data_gd$reltrad_alt), ]
  if (nrow(sub) < 30) next

  P_list_gd[[dw]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_gd[[dw]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_gd[[dw]] = pi_star(P_list_gd[[dw]])
  n_list_gd[[dw]]      = nrow(sub)
}

# Print matrices
for (dw in names(P_list_gd)) {
  cat("\n── GSS window", dw, "  (N =", n_list_gd[[dw]], ") ──\n")
  print(round(P_list_gd[[dw]], 3))
}

# ── IM MEMORY CURVES (t = 0:6) ───────────────────────────────────────────────

im_rows_gd = vector("list", length(P_list_gd))
names(im_rows_gd) = names(P_list_gd)

for (dw in names(P_list_gd)) {
  rows = lapply(0:6, function(t) {
    vals = im_from_P(P_list_gd[[dw]], t = t)
    data.frame(gss_decade = dw, t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_gd[[dw]] = do.call(rbind, rows)
}

im_df_gd = do.call(rbind, im_rows_gd)
im_df_gd$origin     = factor(im_df_gd$origin,     levels = rel_level_order)
im_df_gd$gss_decade = factor(im_df_gd$gss_decade, levels = decade_windows)

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/gss-decade", recursive = TRUE, showWarnings = FALSE)

# Heatmaps with π₀ and π*
for (dw in names(P_list_gd)) {
  safe_dw = gsub("-", "_", dw)
  p = make_combined(
    P_list_gd[[dw]], pi0_list_gd[[dw]], pistar_list_gd[[dw]],
    levels    = rel_level_order,
    title_str = paste0("GSS ", dw, "  (N = ", n_list_gd[[dw]], ")")
  )
  ggsave(paste0("output/figures/gss-decade/trans_", safe_dw, ".png"),
         p, width = 10, height = 7, dpi = 200)
}

# IM curves: faceted by GSS window, colored by origin
p_im_gd = ggplot(im_df_gd, aes(x = t, y = im, color = origin, group = origin)) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_hline(yintercept = log(0.01), linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ gss_decade, nrow = 1) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = NULL,
       title = "Individual Memory by GSS Survey Period (5-state, t = 0–6)") +
  healy_theme

ggsave("output/figures/gss-decade/im_memory_gss_decade.png",
       p_im_gd, width = 14, height = 5, dpi = 200)

# Diagonal persistence across GSS windows
diag_df_gd = do.call(rbind, lapply(names(P_list_gd), function(dw) {
  P = P_list_gd[[dw]]
  data.frame(gss_decade = dw, origin = rownames(P),
             persistence = diag(P), row.names = NULL)
}))
diag_df_gd$origin     = factor(diag_df_gd$origin,     levels = rel_level_order)
diag_df_gd$gss_decade = factor(diag_df_gd$gss_decade, levels = decade_windows)

p_diag_gd = ggplot(diag_df_gd,
    aes(x = gss_decade, y = persistence, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "GSS survey period", y = "Diagonal persistence P[i → i]",
       color = NULL,
       title = "Diagonal Persistence by GSS Survey Period (5-state)") +
  healy_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/gss-decade/diagonal_persistence_gss_decade.png",
       p_diag_gd, width = 8, height = 5, dpi = 200)

# ── YEAR-BY-YEAR DIAGONAL PERSISTENCE ────────────────────────────────────────
# One transition matrix per individual survey year (existing exclusions already
# applied in data_gd). We keep only the diagonal P[i -> i] = share of origin i
# who still report i as adults. Points are suppressed where an origin has fewer
# than `min_origin_n` respondents in that year, to avoid spiky lines from thin
# cells (small-N years such as biennial gaps).

min_origin_n = 25

diag_year = do.call(rbind, lapply(sort(unique(data_gd$year)), function(yr) {
  sub = data_gd[data_gd$year == yr &
                !is.na(data_gd$reltrad16_alt) & !is.na(data_gd$reltrad_alt), ]
  if (nrow(sub) < 30) return(NULL)

  tab      = table(factor(sub$reltrad16_alt, levels = states_alt),
                   factor(sub$reltrad_alt,   levels = states_alt))
  row_tot  = rowSums(tab)
  stay     = diag(tab)
  pii      = ifelse(row_tot > 0, stay / row_tot, NA_real_)

  data.frame(
    year        = as.numeric(yr),
    origin      = states_alt,
    persistence = as.numeric(pii),
    origin_n    = as.numeric(row_tot),
    row.names   = NULL
  )
}))

# Suppress thin-cell points
diag_year$persistence[diag_year$origin_n < min_origin_n] = NA_real_
diag_year$origin = factor(diag_year$origin, levels = rel_level_order)

p_diag_year = ggplot(diag_year,
    aes(x = year, y = persistence, color = origin, group = origin)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(1975, 2020, 5)) +
  labs(x = "GSS survey year", y = "Percent staying  P[i → i]",
       color = NULL,
       title = "Year-by-Year Diagonal Persistence by Religion (5-state)") +
  healy_theme

ggsave("output/figures/gss-decade/diagonal_persistence_by_year.png",
       p_diag_year, width = 10, height = 5, dpi = 200)

# ── YEAR-BY-YEAR STATIONARY DISTRIBUTION (π*) ────────────────────────────────
# One transition matrix per individual survey year, then its implied stationary
# distribution π* (left eigenvector for eigenvalue 1). Unlike diagonal
# persistence, π* is a property of the whole matrix, so a year is dropped
# entirely if ANY origin has fewer than `min_origin_n` respondents (a thin row
# destabilizes the eigen-decomposition). Same exclusions as data_gd throughout.

pistar_year = do.call(rbind, lapply(sort(unique(data_gd$year)), function(yr) {
  sub = data_gd[data_gd$year == yr &
                !is.na(data_gd$reltrad16_alt) & !is.na(data_gd$reltrad_alt), ]
  if (nrow(sub) < 30) return(NULL)

  tab = table(factor(sub$reltrad16_alt, levels = states_alt),
              factor(sub$reltrad_alt,   levels = states_alt))
  if (any(rowSums(tab) < min_origin_n)) return(NULL)   # thin row → skip year

  P  = tab / rowSums(tab)
  class(P) = "matrix"
  ps = pi_star(P)

  data.frame(
    year   = as.numeric(yr),
    origin = names(ps),
    pistar = as.numeric(ps),
    row.names = NULL
  )
}))

pistar_year$origin = factor(pistar_year$origin, levels = rel_level_order)

p_pistar_year = ggplot(pistar_year,
    aes(x = year, y = pistar, color = origin, group = origin)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(1975, 2020, 5)) +
  labs(x = "GSS survey year", y = "Stationary share (π*)",
       color = NULL,
       title = "Year-by-Year Implied Stationary Distribution (π*) by Religion (5-state)") +
  healy_theme

ggsave("output/figures/gss-decade/pistar_by_year.png",
       p_pistar_year, width = 10, height = 5, dpi = 200)

# π* over GSS windows — tracks long-run composition implied by each period's matrix
pistar_df_gd = do.call(rbind, lapply(names(P_list_gd), function(dw) {
  data.frame(gss_decade = dw, origin = names(pistar_list_gd[[dw]]),
             pistar = as.numeric(pistar_list_gd[[dw]]), row.names = NULL)
}))
pistar_df_gd$origin     = factor(pistar_df_gd$origin,     levels = rel_level_order)
pistar_df_gd$gss_decade = factor(pistar_df_gd$gss_decade, levels = decade_windows)

p_pistar_gd = ggplot(pistar_df_gd,
    aes(x = gss_decade, y = pistar, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "GSS survey period", y = "Stationary share (π*)",
       color = NULL,
       title = "Implied Stationary Distribution (π*) by GSS Survey Period") +
  healy_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/gss-decade/pistar_gss_decade.png",
       p_pistar_gd, width = 8, height = 5, dpi = 200)
