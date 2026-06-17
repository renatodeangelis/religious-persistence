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
  select(year, cohort, reltrad, reltrad16, region, evolved, abany, homosex, premarsx, pornlaw, cappun, cappun2) |>
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
    cohort_3   = case_when(
      cohort < 1940                  ~ "1900-1939",
      cohort >= 1940 & cohort < 1960 ~ "1940-1959",
      cohort >= 1960                 ~ paste0("1960-", max(cohort, na.rm = TRUE))
    ),
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

# ── ATTITUDE BINARY RECODES ───────────────────────────────────────────────────
# GSS codings:
#   evolved:   1 = True, 2 = False
#   abany:     1 = Yes (any reason), 2 = No
#   homosex:   1 = Always wrong … 4 = Not wrong at all  → 1–2 conservative, 3–4 liberal
#   premarsx:  1 = Always wrong … 4 = Not wrong at all  → 1–2 conservative, 3–4 liberal
#   pornlaw:   1 = Illegal to all, 2 = Illegal under 18, 3 = Legal to all → 1 conservative, 2–3 liberal

data = data |>
  mutate(
    evolved_bin   = case_when(evolved  == 1 ~ 1L, evolved  == 2 ~ 0L),
    abany_bin     = case_when(abany    == 1 ~ 1L, abany    == 2 ~ 0L),
    homosex_bin   = case_when(homosex  %in% 3:4 ~ 1L, homosex  %in% 1:2 ~ 0L),
    premarsx_bin  = case_when(premarsx %in% 3:4 ~ 1L, premarsx %in% 1:2 ~ 0L),
    pornlaw_bin   = case_when(pornlaw  %in% 2:3 ~ 1L, pornlaw  == 1      ~ 0L),
    # cappun2 covers 1972-73; cappun covers 1974-present; same 1/2 coding
    cappun_merged = coalesce(as.numeric(cappun), as.numeric(cappun2)),
    cappun_bin    = case_when(cappun_merged == 2 ~ 1L, cappun_merged == 1 ~ 0L)
  )

# ── COVERAGE SUMMARY ─────────────────────────────────────────────────────────

att_vars = c(
  "evolved_bin", "abany_bin",
  "homosex_bin", "premarsx_bin", "pornlaw_bin", "cappun_bin"
)

att_coverage = lapply(att_vars, function(v) {
  data |>
    filter(!is.na(reltrad_alt), !is.na(reltrad16_alt)) |>
    summarise(
      variable  = v,
      n_liberal = sum(.data[[v]] == 1L, na.rm = TRUE),
      n_conserv = sum(.data[[v]] == 0L, na.rm = TRUE),
      n_missing = sum(is.na(.data[[v]])),
      pct_cover = round((n_liberal + n_conserv) / n() * 100, 1)
    )
}) |> bind_rows()

print(att_coverage, n = Inf)

# ── CELL-COUNT FEASIBILITY (20-year cohort windows) ──────────────────────────

cell_rows = list()

for (v in att_vars) {
  for (grp in c(0L, 1L)) {
    for (coh in c(1920, 1940, 1960, 1980)) {
      sub = data[
        !is.na(data$cohort_20) & data$cohort_20 == coh &
        !is.na(data[[v]])      & data[[v]]       == grp, ]
      if (nrow(sub) == 0) next
      tab = table(
        factor(sub$reltrad16_alt, levels = states_alt),
        factor(sub$reltrad_alt,   levels = states_alt)
      )
      cell_rows[[length(cell_rows) + 1]] = data.frame(
        variable  = v,
        group     = if (grp == 1L) "liberal" else "conservative",
        cohort    = coh,
        n         = nrow(sub),
        min_cell  = min(tab),
        cells_lt5 = sum(tab < 5),
        cells_0   = sum(tab == 0),
        row.names = NULL
      )
    }
  }
}

cell_df = do.call(rbind, cell_rows)

print(cell_df, row.names = FALSE)

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
print(pooled_summary, row.names = FALSE)

# ── ATTITUDE MATRICES: THREE-COHORT STRATIFICATION ───────────────────────────
# Groups: pre-1940 / 1940-1959 / 1960+ (open-ended)

max_cohort_yr   = max(data$cohort, na.rm = TRUE)
cohort_3_levels = c("1900-1939", "1940-1959", paste0("1960-", max_cohort_yr))

P_list_att_3coh   = list()
pi0_list_att_3coh = list()
coh3_summary_rows = list()

for (v in att_vars) {
  for (grp in c(0L, 1L)) {
    for (cg in cohort_3_levels) {
      sub = data[!is.na(data$reltrad_alt) & !is.na(data$reltrad16_alt) &
                 !is.na(data[[v]])        & data[[v]]     == grp &
                 !is.na(data$cohort_3)   & data$cohort_3 == cg, ]
      if (nrow(sub) < 30) next
      grp_lbl = if (grp == 1L) "liberal" else "conservative"
      key     = paste(v, grp_lbl, cg, sep = "_")

      P_list_att_3coh[[key]]   = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
      pi0_list_att_3coh[[key]] = pi_0(sub, "reltrad16_alt")

      coh3_summary_rows[[key]] = data.frame(
        variable    = v,
        group       = grp_lbl,
        cohort_grp  = cg,
        n           = nrow(sub),
        mean_cohort = round(mean(sub$cohort, na.rm = TRUE), 1),
        row.names   = NULL
      )
    }
  }
}

coh3_summary = do.call(rbind, coh3_summary_rows)
print(coh3_summary, row.names = FALSE)

# ── ATTITUDE FIGURES ──────────────────────────────────────────────────────────

att_labels = c(
  evolved_bin  = "Evolution",
  abany_bin    = "Abortion (any reason)",
  homosex_bin  = "Homosexuality",
  premarsx_bin = "Premarital sex",
  pornlaw_bin  = "Pornography law",
  cappun_bin   = "Capital punishment"
)

dir.create("output/figures/attitude/pooled", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/attitude/3coh",   recursive = TRUE, showWarnings = FALSE)

for (v in att_vars) {
  for (grp_lbl in c("conservative", "liberal")) {
    key = paste(v, grp_lbl, sep = "_")
    if (is.null(P_list_att_pooled[[key]])) next
    P       = P_list_att_pooled[[key]]
    pi0     = pi0_list_att_pooled[[key]]
    pistar  = pi_star(P)
    med_coh = pooled_summary$median_cohort[pooled_summary$variable == v &
                                             pooled_summary$group    == grp_lbl]
    ttl     = paste0(att_labels[v], " — ", tools::toTitleCase(grp_lbl), " (pooled)")
    p       = make_combined(P, pi0, pistar, levels = rel_level_order, title_str = ttl) +
                patchwork::plot_annotation(subtitle = paste0("Median birth cohort: ", med_coh))
    ggsave(paste0("output/figures/attitude/pooled/", key, ".png"),
           p, width = 10, height = 7, dpi = 200)
  }
}

for (v in att_vars) {
  for (grp_lbl in c("conservative", "liberal")) {
    for (cg in cohort_3_levels) {
      key = paste(v, grp_lbl, cg, sep = "_")
      if (is.null(P_list_att_3coh[[key]])) next
      P      = P_list_att_3coh[[key]]
      pi0    = pi0_list_att_3coh[[key]]
      pistar = pi_star(P)
      ttl    = paste0(att_labels[v], " — ", tools::toTitleCase(grp_lbl), " — ", cg)
      p      = make_combined(P, pi0, pistar, levels = rel_level_order, title_str = ttl)
      safe_cg = gsub("[^a-zA-Z0-9]", "_", cg)
      ggsave(paste0("output/figures/attitude/3coh/", v, "_", grp_lbl, "_", safe_cg, ".png"),
             p, width = 10, height = 7, dpi = 200)
    }
  }
}

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

# ── BINARY (AFFILIATED / UNAFFILIATED) 2×2 MATRICES ──────────────────────────
# Affiliated = any reltrad_alt other than "none"; unaffiliated = "none"

data = data |>
  mutate(
    belief   = if_else(reltrad_alt   == "none", "unaffiliated", "affiliated"),
    belief16 = if_else(reltrad16_alt == "none", "unaffiliated", "affiliated")
  )

states_2x2 = c("affiliated", "unaffiliated")

cohorts_10_2x2 = sort(unique(
  data$cohort_10[!is.na(data$cohort_10) & data$cohort_10 >= 1920 & data$cohort_10 <= 1980]
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

# Print matrices with N
for (key in names(P_list_2x2)) {
  cat("\n── Cohort", key, "–", as.integer(key) + 9, "  (N =", n_list_2x2[[key]], ") ──\n")
  print(round(P_list_2x2[[key]], 3))
}

# ── FIGURES: individual make_combined heatmaps per cohort ────────────────────

dir.create("output/figures/binary", recursive = TRUE, showWarnings = FALSE)

for (key in names(P_list_2x2)) {
  p = make_combined(
    P_list_2x2[[key]], pi0_list_2x2[[key]], pistar_list_2x2[[key]],
    levels    = states_2x2,
    title_str = paste0("Cohort ", key, "–", as.integer(key) + 9,
                       "  (N = ", n_list_2x2[[key]], ")")
  )
  ggsave(paste0("output/figures/binary/trans_", key, "_10yr_2x2.png"),
         p, width = 7, height = 5, dpi = 200)
}

# ── FIGURE: all 4 cells over cohorts ─────────────────────────────────────────

cells_2x2 = do.call(rbind, lapply(names(P_list_2x2), function(key) {
  P   = P_list_2x2[[key]]
  df  = as.data.frame(as.table(P), stringsAsFactors = FALSE)
  names(df) = c("origin", "current", "prob")
  df$cohort = as.integer(key)
  df
}))

p_cells_2x2 = ggplot(cells_2x2,
    aes(x = cohort, y = prob, color = current, group = current)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ origin, nrow = 1,
             labeller = labeller(origin = c(
               affiliated   = "Origin: Affiliated",
               unaffiliated = "Origin: Unaffiliated"))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    values = c(affiliated = "#0072B2", unaffiliated = "#D55E00"),
    labels = c(affiliated = "→ Affiliated", unaffiliated = "→ Unaffiliated")) +
  labs(x = "Birth cohort (10-year bins)", y = "Transition probability",
       color = NULL,
       title = "2×2 Transition Probabilities by Birth Cohort") +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text = element_text(size = 11))

ggsave("output/figures/binary/cells_2x2_10yr.png", p_cells_2x2,
       width = 10, height = 5, dpi = 200)

# ── FIGURE: diagonal persistence over cohorts ─────────────────────────────────

persistence_2x2 = do.call(rbind, lapply(names(P_list_2x2), function(key) {
  P = P_list_2x2[[key]]
  data.frame(cohort = as.integer(key), state = rownames(P),
             persistence = diag(P), row.names = NULL)
}))

p_persist_2x2 = ggplot(persistence_2x2,
    aes(x = cohort, y = persistence, color = state, group = state)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_color_manual(
    values = c(affiliated = "#0072B2", unaffiliated = "#D55E00")) +
  labs(x = "Birth cohort (10-year bins)",
       y = "Diagonal persistence P[i → i]",
       color = "Origin state",
       title = "Diagonal Persistence: Affiliated vs. Unaffiliated (10-year cohorts)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/figures/binary/persistence_2x2_10yr.png", p_persist_2x2,
       width = 8, height = 5, dpi = 200)

# ── ATTITUDE TRENDS BY RELTRAD AND BIRTH COHORT ──────────────────────────────
# Conservative positions: evolved_bin == 0, abany_bin == 0, homosex_bin == 0

dir.create("output/figures/attitude/cohort_trends", recursive = TRUE, showWarnings = FALSE)

reltrad_colors_att = c(
  "Catholic"    = "#0072B2",
  "Evangelical" = "#D55E00",
  "Mainline"    = "#009E73",
  "Other"       = "#CC79A7",
  "None"        = "#E69F00"
)

att_cohort_df = data |>
  filter(!is.na(reltrad_alt), cohort_5 >= 1940, cohort_5 <= 1980) |>
  pivot_longer(
    cols      = c(evolved_bin, abany_bin, homosex_bin, cappun_bin),
    names_to  = "attitude",
    values_to = "liberal"
  ) |>
  filter(!is.na(liberal)) |>
  group_by(attitude, reltrad_alt, cohort_5) |>
  summarise(
    n                = n(),
    pct_conservative = mean(liberal == 0L) * 100,
    se               = sqrt(pct_conservative / 100 * (1 - pct_conservative / 100) / n) * 100,
    .groups          = "drop"
  ) |>
  mutate(
    attitude = factor(attitude,
      levels = c("evolved_bin", "abany_bin", "homosex_bin", "cappun_bin"),
      labels = c("Evolution (% deny)", "Abortion (% oppose)", "Homosexuality (% morally wrong)",
                 "Capital punishment (% favor)")
    ),
    reltrad_alt = factor(reltrad_alt,
      levels = c("catholic", "evangelical", "mainline", "other", "none"),
      labels = c("Catholic", "Evangelical", "Mainline", "Other", "None")
    )
  )

make_att_cohort_plot = function(att_label) {
  att_cohort_df |>
    filter(attitude == att_label) |>
    ggplot(aes(x = cohort_5, y = pct_conservative, color = reltrad_alt, fill = reltrad_alt, group = reltrad_alt)) +
    geom_ribbon(aes(ymin = pct_conservative - 1.96 * se, ymax = pct_conservative + 1.96 * se),
                alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    scale_color_manual(values = reltrad_colors_att, name = NULL) +
    scale_fill_manual(values = reltrad_colors_att, name = NULL) +
    scale_x_continuous(breaks = seq(1940, 1980, by = 10)) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    labs(title = att_label, x = NULL, y = sub(".*\\((.+)\\)", "\\1", att_label)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(size = 12, face = "bold"))
}

p_evolved = make_att_cohort_plot("Evolution (% deny)")
p_abany   = make_att_cohort_plot("Abortion (% oppose)")
p_homosex = make_att_cohort_plot("Homosexuality (% morally wrong)")
p_cappun  = make_att_cohort_plot("Capital punishment (% favor)")

caption_str = "Source: GSS. 5-year birth cohort bins, 1940–1980. Ribbons show 95% CIs. Current religious affiliation (reltrad_alt)."

p_evolved = p_evolved + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)
p_abany   = p_abany   + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)
p_homosex = p_homosex + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)
p_cappun  = p_cappun  + labs(x = "Birth cohort (5-year bin)") +
  plot_annotation(caption = caption_str)

ggsave("output/figures/attitude/cohort_trends/att_trends_evolution.png",
       p_evolved, width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends/att_trends_abortion.png",
       p_abany,   width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends/att_trends_homosexuality.png",
       p_homosex, width = 7, height = 5, dpi = 200)
ggsave("output/figures/attitude/cohort_trends/att_trends_cappun.png",
       p_cappun,  width = 7, height = 5, dpi = 200)
