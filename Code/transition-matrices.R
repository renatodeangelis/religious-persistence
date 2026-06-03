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

# ── DATA ──────────────────────────────────────────────────────────────────────

data(gss_all)
data = gss_all |>
  select(year, cohort, reltrad, reltrad16, region) |>
  filter(!(year %in% c(1972, 2021))) |>
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
    cohort_5   = floor((cohort - 1900) / 5)  * 5  + 1900,
    cohort_10  = floor((cohort - 1900) / 10) * 10 + 1900,
    cohort_20  = floor((cohort - 1900) / 20) * 20 + 1900,
    region_broad = case_when(
      as.numeric(region) == 1 ~ "Northeast",
      as.numeric(region) == 2 ~ "Midwest",
      as.numeric(region) == 3 ~ "South",
      as.numeric(region) == 4 ~ "West",
      TRUE ~ NA_character_
    )
  )

# ── STATE SPACE ───────────────────────────────────────────────────────────────

states_alt = sort(unique(c(data$reltrad_alt, data$reltrad16_alt)))
states_alt = states_alt[!is.na(states_alt)]

# ── NATIONAL COHORT MATRICES ─────────────────────────────────────────────────

# ── 5-year cohort loop ───────────────────────────────────────────────────────
cohorts_5 = sort(unique(data$cohort_5[!is.na(data$cohort_5) & data$cohort_5 >= 1920 & data$cohort_5 <= 1980]))

P_list_5      = list()
pi0_list_5    = list()
pistar_list_5 = list()

for (coh in cohorts_5) {
  sub = data[!is.na(data$cohort_5) & data$cohort_5 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_5[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_5[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_5[[key]] = pi_star(P_list_5[[key]])
}

# ── 10-year cohort loop ──────────────────────────────────────────────────────
cohorts_10 = sort(unique(data$cohort_10[!is.na(data$cohort_10) & data$cohort_10 >= 1920 & data$cohort_10 <= 1980]))

P_list_10      = list()
pi0_list_10    = list()
pistar_list_10 = list()

for (coh in cohorts_10) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) next
  key = as.character(coh)

  P_list_10[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0_list_10[[key]]    = pi_0(sub, "reltrad16_alt")
  pistar_list_10[[key]] = pi_star(P_list_10[[key]])
}

# ── 20-year cohort loop ──────────────────────────────────────────────────────
cohorts_20_pooled = sort(unique(data$cohort_20[!is.na(data$cohort_20) & data$cohort_20 >= 1920 & data$cohort_20 <= 1980]))

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

# ── NATIONAL IM COMPUTATION ──────────────────────────────────────────────────

# ── IM LOOP (10-year cohorts, t = 0:4) ──────────────────────────────────────
im_rows_10 = vector("list", length(P_list_10))
names(im_rows_10) = names(P_list_10)

for (key in names(P_list_10)) {
  rows = lapply(0:4, function(t) {
    vals = im_from_P(P_list_10[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_10[[key]] = do.call(rbind, rows)
}

im_df_10 = do.call(rbind, im_rows_10)

# ── IM LOOP (20-year cohorts, t = 0:4) ──────────────────────────────────────
im_rows_20 = vector("list", length(P_list_20))
names(im_rows_20) = names(P_list_20)

for (key in names(P_list_20)) {
  rows = lapply(0:4, function(t) {
    vals = im_from_P(P_list_20[[key]], t = t)
    data.frame(cohort = as.integer(key), t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_20[[key]] = do.call(rbind, rows)
}

im_df_20 = do.call(rbind, im_rows_20)

# ── NATIONAL MOBILITY ────────────────────────────────────────────────────────

mob_rows = lapply(1920:1980, function(coh) {
  sub = data[!is.na(data$cohort) & data$cohort == coh &
               !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
  if (nrow(sub) < 30) return(NULL)
  P   = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  pi0 = pi_0(sub, "reltrad16_alt")
  data.frame(cohort = coh, mobility = overall_mobility(P, pi0))
})
mob_df = do.call(rbind, Filter(Negate(is.null), mob_rows))

# ── EXCHANGE/STRUCTURAL MOBILITY (10-year cohorts, t = 0:4) ──────────────────

em_sm_rows_10 = lapply(names(P_list_10), function(key) {
  P   = P_list_10[[key]]
  pi0 = pi0_list_10[[key]]
  lapply(0:4, function(t) {
    om_v = overall_mobility(P, mu_t(pi0, P, t))
    sm_v = sm(P, pi0, t)
    data.frame(cohort = as.integer(key), t = t, EM = om_v - sm_v, SM = sm_v, OM = om_v,
               row.names = NULL)
  })
})
em_sm_df_10 = do.call(rbind, unlist(em_sm_rows_10, recursive = FALSE))

# ── NATIONAL FIGURES ──────────────────────────────────────────────────────────

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_5)) {
  p = make_combined(P_list_5[[key]], pi0_list_5[[key]], pistar_list_5[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", key, "–", as.integer(key) + 4))
  ggsave(paste0("output/figures/trans_", key, "_5yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

for (key in names(P_list_10)) {
  p = make_combined(P_list_10[[key]], pi0_list_10[[key]], pistar_list_10[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", key, "–", as.integer(key) + 9))
  ggsave(paste0("output/figures/trans_", key, "_10yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

for (key in names(P_list_20)) {
  p = make_combined(P_list_20[[key]], pi0_list_20[[key]], pistar_list_20[[key]],
                    levels = rel_level_order,
                    title_str = paste0("Cohort ", key, "–", as.integer(key) + 19))
  ggsave(paste0("output/figures/trans_", key, "_20yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

im_df_10$origin = factor(im_df_10$origin, levels = rel_level_order)

p_im = ggplot(im_df_10, aes(x = t, y = im, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = "Origin", title = "Individual Memory by Cohort (10-year bins, t = 0–4)") +
  theme_minimal() +
  theme(
    strip.text      = element_text(size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

ggsave("output/figures/im_memory_10yr.png", p_im, width = 12, height = 7, dpi = 200)

im_df_20$origin = factor(im_df_20$origin, levels = rel_level_order)

p_im_20 = ggplot(im_df_20, aes(x = t, y = im, color = origin, group = origin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Step (t)", y = "log(TV distance from π*)",
       color = "Origin", title = "Individual Memory by Cohort (20-year bins, t = 0–4)") +
  theme_minimal() +
  theme(
    strip.text      = element_text(size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

ggsave("output/figures/im_memory_20yr.png", p_im_20, width = 10, height = 5, dpi = 200)

p_mob = ggplot(mob_df, aes(x = cohort, y = mobility)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, span = 0.4) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Birth cohort", y = "Overall mobility (1 − weighted diagonal)",
       title = "Religious Mobility by Birth Cohort (pooled, 1-year bins)") +
  theme_minimal()

ggsave("output/figures/mobility_pooled.png", p_mob, width = 9, height = 5, dpi = 200)

em_sm_long = pivot_longer(em_sm_df_10, cols = c(EM, OM),
                           names_to = "measure", values_to = "value")

p_em_sm = ggplot(em_sm_long, aes(x = t, y = value, color = measure,
                                  linetype = measure, shape = measure, group = measure)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ cohort, nrow = 2) +
  scale_color_manual(name   = "Measure",
                     values = c(EM = "#D55E00", OM = "#0072B2"),
                     labels = c(EM = "Exchange mobility", OM = "Overall mobility")) +
  scale_linetype_manual(name   = "Measure",
                        values = c(EM = "dashed", OM = "solid"),
                        labels = c(EM = "Exchange mobility", OM = "Overall mobility")) +
  scale_shape_manual(name   = "Measure",
                     values = c(EM = 17, OM = 16),
                     labels = c(EM = "Exchange mobility", OM = "Overall mobility")) +
  scale_x_continuous(breaks = 0:4) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Step (t)", y = "Probability to move",
       title = "Exchange and Structural Mobility by Cohort (10-year bins, t = 0–4)") +
  theme_minimal() +
  theme(strip.text      = element_text(size = 12),
        legend.position = "bottom",
        legend.title    = element_text(size = 12),
        legend.text     = element_text(size = 11))

ggsave("output/figures/em_sm_10yr.png", p_em_sm, width = 12, height = 7, dpi = 200)

# ── REGIONAL COHORT MATRICES ─────────────────────────────────────────────────

cohorts_20    = c(1920, 1940, 1960)
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

# ── REGIONAL IM AND MOBILITY COMPUTATION ─────────────────────────────────────

im_rows_reg = vector("list", length(P_list_reg))
names(im_rows_reg) = names(P_list_reg)

for (key in names(P_list_reg)) {
  parts = strsplit(key, "_")[[1]]
  reg   = paste(parts[-length(parts)], collapse = "_")
  coh   = as.integer(parts[length(parts)])
  rows  = lapply(0:4, function(t) {
    vals = im_from_P(P_list_reg[[key]], t = t)
    data.frame(region = reg, cohort = coh, t = t, origin = names(vals), im = vals,
               row.names = NULL)
  })
  im_rows_reg[[key]] = do.call(rbind, rows)
}

im_df_reg = do.call(rbind, im_rows_reg)

mob_reg_rows = lapply(regions_broad, function(reg) {
  lapply(1920:1980, function(coh) {
    sub = data[!is.na(data$cohort)        & data$cohort        == coh &
                 !is.na(data$region_broad) & data$region_broad  == reg &
                 !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
    if (nrow(sub) < 30) return(NULL)
    P   = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
    pi0 = pi_0(sub, "reltrad16_alt")
    data.frame(cohort = coh, region = reg, mobility = overall_mobility(P, pi0))
  })
})
mob_reg_df = do.call(rbind,
  Filter(Negate(is.null), unlist(mob_reg_rows, recursive = FALSE)))

# ── REGIONAL FIGURES ──────────────────────────────────────────────────────────

dir.create("output/figures/region", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_reg)) {
  coh_yr  = as.integer(sub(".*_(\\d{4})$", "\\1", key))
  reg_lbl = gsub("_", " ", sub("_\\d{4}$", "", key))
  p = make_combined(P_list_reg[[key]], pi0_list_reg[[key]], pistar_list_reg[[key]],
                    levels = rel_level_order,
                    title_str = paste0(reg_lbl, " – ", coh_yr, "–", coh_yr + 19))
  ggsave(paste0("output/figures/region/trans_", key, "_20yr.png"), p,
         width = 10, height = 7, dpi = 200)
}

for (reg in regions_broad) {
  reg_plots = lapply(cohorts_20, function(coh) {
    key = paste(reg, coh, sep = "_")
    if (!is.null(P_list_reg[[key]])) {
      make_combined(
        P_list_reg[[key]], pi0_list_reg[[key]], pistar_list_reg[[key]],
        levels    = rel_level_order,
        title_str = paste0(reg, "\n", coh, "–", coh + 19)
      )
    } else {
      patchwork::plot_spacer()
    }
  })
  p_reg = patchwork::wrap_plots(reg_plots, nrow = 1, ncol = length(cohorts_20))
  fname = paste0("output/figures/region/trans_grid_", gsub(" ", "_", reg), "_20yr.png")
  ggsave(fname, p_reg, width = 30, height = 7, dpi = 200)
}

mob_reg_df$region = factor(mob_reg_df$region, levels = regions_broad)

p_mob_reg = ggplot(mob_reg_df, aes(x = cohort, y = mobility)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, span = 0.5) +
  facet_wrap(~ region, nrow = 2) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Birth cohort", y = "Overall mobility (1 − weighted diagonal)",
       title = "Religious Mobility by Birth Cohort and Region (1-year bins)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12))

ggsave("output/figures/region/mobility_region.png", p_mob_reg,
       width = 11, height = 7, dpi = 200)
